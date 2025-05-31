library(tidyverse)
library(seqinr)
library(viridis)
library(ape)
library(vegan)

set.seed(123)

#########################################
########### INITIAL PROCESS: ZOOX #####
#########################################
#read in blast assignments
blastout <- read.table("data/PoritesASVsBlastOut.tab")
colnames(blastout) <- c("query", "hit", "eval", "match_length", "percent_iden_match")
length(unique(blastout$query))
blastout<- blastout %>% filter(eval < 1e-100) %>% filter(percent_iden_match > 95) %>% separate(hit, into = c("Clade", "Subclade", "substrain"), sep="(?=[a-z])", remove = F) %>% dplyr::select(!c(Subclade, substrain)) %>% separate(Clade, into =c("Clade", "Subclade"), sep="(?<=[A-Z])")
length(unique(blastout$query)) #only 244 great hits. Let's roll with those
blast_best <- blastout %>% group_by(query) %>% filter(percent_iden_match==max(percent_iden_match))
length(unique(blast_best$query)) #more than one best hit per read. Let's get consensus.
table(blastout$hit)

good_hits_by_query <- blastout %>% group_by(query) %>% summarise(totalhits=n())
blast_consensus_all <- blastout %>% group_by(Clade, Subclade, query) %>% summarise(count_sub=n()) %>% ungroup()
blast_consensus_all <- left_join(blast_consensus_all, good_hits_by_query) %>% mutate(perc_hit =count_sub/totalhits) %>% ungroup() %>% group_by(query) %>% filter(perc_hit==max(perc_hit)) 
hist(blast_consensus_all$perc_hit)

best_hits_by_query <- blast_best %>% group_by(query) %>% summarise(totalhits=n())
blast_consensus_best <- blast_best %>% group_by(Clade, Subclade, query) %>% summarise(count_sub=n()) %>% ungroup()
blast_consensus_best <- left_join(blast_consensus_best, best_hits_by_query) %>% mutate(perc_hit =count_sub/totalhits) %>% ungroup() %>%
  group_by(query) %>% filter(perc_hit==max(perc_hit)) 
hist(blast_consensus_best$perc_hit) 
keep_asv_ass <- blast_consensus_best %>% filter(perc_hit > .9) %>% separate(query, into = c("ASV", "SEQ"))

#grab the ITS2 OTU table
itsin<-readRDS('data/UpdatedSeqTabITS2.Rdata')
its2.counts <- data.frame(itsin)
new_ids <- c()
for(i in 1:length(its2.counts$SampleID)){
  x=ifelse(its2.counts$SampleID[i] %in% new_ids, paste(its2.counts$SampleID[i], "_1", sep = ""), its2.counts$SampleID[i])
  x=ifelse(x %in% new_ids, gsub("_1", "_2", x), x)
  new_ids <- c(new_ids, x)
}
rownames(its2.counts) <- new_ids
its2.counts$SampleID <- NULL

#################get co-occuring ASVS##############

its2.pa <- its2.counts
its2.pa[its2.pa > 0] <- 1
library(pheatmap)
its2.cor <- cor(its2.pa)
png("figures/SuppFig_ZooxCorrCutoff.png", width = 6.5, height= 4, units = "in", res= 500)
hist(its2.cor[its2.cor>0], main="Cross-Correlation of Symbiont ASVs", xlab="Correlation Value") #set a threshhold at .8
abline(v=.8, col = 'red')
dev.off()

pheatmap(cor(its2.pa), show_rownames = F, show_colnames = F)
plot(hclust(as.dist(cor(its2.pa))), labels = F)
its2.df.cor <- data.frame(its2.cor)
tomerge <- list()
i = 1
cols_remaining = ncol(its2.df.cor)
while(i < cols_remaining){
  grouped_seqs <-data.frame(SEQ=c(colnames(its2.df.cor)[i], rownames(its2.df.cor)[which(its2.df.cor[,i]>.8)]))
  #assignments <- keep_asv_ass %>% filter(SEQ%in%grouped_seqs) %>% dplyr::select(Clade, Subclade, SEQ)
  its2.df.cor <- its2.df.cor[setdiff(rownames(its2.df.cor), grouped_seqs$SEQ), setdiff(colnames(its2.df.cor), grouped_seqs$SEQ)]
  cols_remaining <- ncol(its2.df.cor)
  i = i + 1
  tomerge[[i]] <- left_join(grouped_seqs, keep_asv_ass) %>% dplyr::select(SEQ, Clade, Subclade)
}
length(tomerge)
#how many groups?
for(entry in 1:length(tomerge)){
  print(length(tomerge[[entry]]))
}

ASV_out <- bind_rows(tomerge, .id="ASV_group")
ASV_consensus <- ASV_out %>% group_by(ASV_group, Clade, Subclade) %>% summarise(hits=n()) %>% ungroup() 
group_totals <- ASV_out %>% group_by(ASV_group) %>% summarise(total_entries=n())
consensus_groups <- left_join(ASV_consensus, group_totals) %>% mutate(perc_covd=hits/total_entries) %>% filter(perc_covd==max(perc_covd)) %>% mutate(consensus_call=ifelse(perc_covd > .8, paste(Clade, Subclade, sep = "_"), Clade))
#do groups have consensus?
table(consensus_groups$consensus_call)

######### get a new abundance table by group
#its2.counts
outdf <- data.frame()
for(i in 1:length(unique(ASV_out$ASV_group))){
  group_asvs <- ASV_out %>% filter(ASV_group == unique(ASV_out$ASV_group)[i])
  ncol(its2.counts[,group_asvs$SEQ]) == length(group_asvs$SEQ)
  newname=paste("ASV_group_", i, sep ="")
  tmp <- t(data.frame(rowSums(its2.counts[,group_asvs$SEQ])))
  rownames(tmp) <- newname
  outdf <- rbind(outdf, tmp)
}

out.pa <- outdf
out.pa[out.pa>0]=1
hist(log10(rowSums(out.pa)))
keeprows <- rownames(outdf)[rowSums(out.pa) >= 3] #must be present in at least 3 samps
keepcols <- colnames(outdf)[colSums(outdf) >=100] #trash samples with fewer than 100 total ASV counts
keepcols <- intersect(keepcols, colnames(outdf)[colSums(out.pa) > 1])
keeprows <- intersect(keepcols, colnames(outdf)[rowSums(outdf) > 100])

its2.filtered <- outdf[keeprows,keepcols] #got rid of 50 entries

meta = read.csv("data/Updated_Porites_Meta.csv", header = T)
meta <- meta %>% mutate(Site=ifelse(Site=="sOrpheus2", "sOrpheus", Site))
rownames(meta) <- meta$ID

#############################################
### combine across replicates
its2.combinereps <- data.frame(t(outdf))
its2.combinereps$ID <- rownames(its2.combinereps)
its2.combinereps <- its2.combinereps %>% pivot_longer(starts_with("ASV"), names_to = "ASV_group", values_to="readcount")
its2.combinereps <- its2.combinereps %>% left_join(meta) %>% group_by(ColonyID, ASV_group) %>% summarise(summed_count=sum(readcount)) %>%
  pivot_wider(names_from=ASV_group, values_from=summed_count) %>% drop_na()
its2.combinereps <- data.frame(its2.combinereps)
rownames(its2.combinereps) <- its2.combinereps$ColonyID
its2.combinereps$ColonyID <- NULL
#filter by presense absense in colonies....
out.pa <- its2.combinereps
out.pa[out.pa>0]=1
hist(log10(rowSums(out.pa)))
keeprows <- rownames(its2.combinereps)[rowSums(out.pa) >= 3] #must be present in at least 3 COLONIES

keepcols <- colnames(its2.combinereps)[colSums(its2.combinereps) >=100] #trash samples with fewer than 100 total ASV counts
keepcols <- intersect(keepcols, colnames(its2.combinereps)[colSums(out.pa) > 1])
keeprows <- intersect(keeprows, rownames(its2.combinereps)[rowSums(its2.combinereps) > 1000])

rowSums(outdf) 
its2.filtered.colony <- its2.combinereps[keeprows,keepcols]

bycolony.meta <- meta %>% dplyr::select(ColonyID, Site, Age) %>% mutate(Site=gsub("2", "", Site)) %>% unique()


#############################################

its2.total <- data.frame(decostand(data.frame(its2.filtered.colony), method = "total"))
its2.total$ID <- rownames(its2.total)
its2.total %>% pivot_longer(starts_with("ASV"), names_to="ASV_group", values_to="perc_sample") %>%
  separate(ASV_group, into = c("asv", "group", "number")) %>% 
  left_join(consensus_groups, by=c("number"="ASV_group")) %>%
  ggplot(aes(x=ID, y=perc_sample, fill=number)) + geom_bar(stat="identity", width =1) + theme_bw() + scale_fill_viridis_d()


#######################################################################
####### HOST GENETICS & ADMIXTURE  ####################################
#######################################################################

dev.off()
ibs <- read.table("data/clone_qual_filteredPorites.ibsMat")
nam <- read.table("data/clone_qual_filteredPorites.names")
nam$V1 <- gsub("E75-1", "E75", nam$V1)
rownames(ibs) <- nam$V1
colnames(ibs) <- nam$V1
ibs <- as.dist(ibs)
plot(hclust(ibs))
abline(h=.1)

ad5 <- read.table("data/clone_qual_filteredPorites.admix.5.Q")
ad5$SampleID <- nam$V1
ad5 <- ad5 %>% separate(SampleID, into =c("ID"), sep = "\\.") %>% separate(ID, into=c("smallID", "techrep"), remove = F)
admix_prop <- test %>% pivot_longer(!c(ID, smallID, techrep), names_to = "AdmixGroup", values_to = "Prop")
admix_prop <- admix_prop %>% left_join(meta, by=c("smallID"="ID")) %>% filter(is.na(Site)==F)

admix_prop <- admix_prop %>% mutate(short_site=case_when(Site=="nePelorus"~"neP", 
                                           Site=="pioneerBay" ~"PB",
                                           Site=="sePelorus"~"seP",
                                           Site=="sOrpheus" ~ "sO",
                                           Site =="sOrpheus2"~"sO")) %>%
                      mutate(short_age=case_when(Age=="adult"~"A",
                                          Age=="juvenile"~"J")) %>% mutate(short_id=paste(short_site, short_age, sep="-"))

max_ad <- admix_prop %>% group_by(ID) %>% filter(Prop==max(Prop))

plot_levels <- max_ad %>% arrange(AdmixGroup,Site, Age, smallID)
id_levels <- unique(plot_levels$ID)
admix_prop$ID <- factor(admix_prop$ID, levels=id_levels)

uq_samps <- admix_prop %>% dplyr::select(Site, Age, ID, smallID, techrep) %>% unique()
add_site_verts <- cumsum(table(uq_samps$Site)) - 1
add_site_verts <- add_site_verts[1:(length(add_site_verts))]+.5
add_age_verts <- cumsum(table(uq_samps$Age, uq_samps$Site)) - 1
add_age_verts <- add_age_verts[1:(length(add_age_verts))]+.5
admix_prop %>% ggplot(aes(x=ID,y=Prop, fill=AdmixGroup)) + geom_bar(stat="identity", width = .95) + scale_fill_viridis_d(option="turbo") + theme_classic() +
  theme(axis.text.x=element_blank(), legend.position = "none", axis.text =element_text(size = 8),
        axis.title=element_text(size = 10)) + scale_x_discrete(labels=admix_prop$short_id) + 
  xlab("Sample") + ylab("Prop. ADMIX Group")# +

#ggsave("results/ADMIX5.MEREVIEW.png", width = 3, height = 1.75, units = "in", dpi = 500)

#make little wheels for the map to show site differences. 
groups <- unique(max_ad$AdmixGroup)
for(i in 1:length(groups)){
  sub <- max_ad %>% filter(AdmixGroup == groups[i])
  tmp <- gsub("E75", "E75-1", paste(sub$ID, ".trim.master.local.coral.sorted.bam",sep=""))
  write.table(tmp, file=paste("data/bamlist.", groups[i], ".admix", sep = ""), row.names = F, col.names = F, quote = F)
}


#####Unrooted Tree of Admixture FST ####
fst <- read.delim("data/named.fst.admix5", header = F)
FST_matrix <- fst %>% separate(V1, into = c("ad", "popA", "popB", "fststat")) %>%
  dplyr::select(popA, popB, V3) %>% drop_na() %>%
  mutate(V3 = ifelse(popA == popB, 0, V3)) %>%
  mutate(popA = paste("Admix_", popA, sep =""),
         popB = paste("Admix_", popB, sep = "")) %>%
  pivot_wider(names_from = popB, values_from = V3)

FST_matrix <- data.frame(FST_matrix)
rownames(FST_matrix) = FST_matrix$popA
FST_matrix$popA = NULL
FST_matrix = as.dist(FST_matrix)
plot(hclust(FST_matrix))
f.hc = hclust(FST_matrix, method = "ave")

tr <- as.phylo(f.hc)
par(omi = c(0.1, 0.1, 0.1, 0.1), 
    mai= c(0, 0, 0, 0))
pdf("figures/unrootedFSTtree.pdf", width = 5, height = 3.5)
plot(tr, type="unrooted", rotate.tree = 100, edge.width = 2, cex = 1)
add.scale.bar(x = 0.25, y = -.01)
dev.off()

png(filename = "figures/Supp_PairedFST_5clust.png", width = 5, height = 4, units = "in", res = 500)
pheatmap(as.matrix(FST_matrix),show_rownames = T,
         show_colnames = T, display_numbers = T,
         cluster_rows = T, 
         cluster_cols = T, 
         na_col = "white", fontsize_number = 10)
dev.off()

###### PCoA of Genetic Distances ###########

ibs.cap <- capscale(ibs~1)

ibs.scores <- scores(ibs.cap, tidy = T, choices=c(1:3), scaling = "species") %>% filter(score=="sites")

ibs.meta <- ibs.scores %>% left_join(max_ad, by=c("label"="ID"))

ggplot(ibs.meta, aes(x=MDS1, y=MDS2, col=AdmixGroup)) +
  geom_vline(xintercept = 0, lty='dashed', col = 'grey')+
  geom_hline(yintercept = 0, lty='dashed', col = 'grey')+
  geom_point(size =4, alpha =.7) + theme_classic() + coord_equal() + scale_color_viridis_d(option="turbo")


adonis2(ibs~ibs.meta$Site*ibs.meta$Age, method = "marginal") #yeah, there's sig structure by site, age. Mostly Site.

#add a label instead of a legend to plot
cent_label = ibs.meta %>% dplyr::select(MDS1, MDS2, AdmixGroup) %>%
  group_by(AdmixGroup) %>% summarise(centx=mean(MDS1), centy=mean(MDS2)) %>% mutate(AdmixGroup=gsub("V", "Admix_", AdmixGroup))
library(ggrepel)
ggplot() +
  geom_vline(xintercept = 0, lty='dashed', col = 'grey')+
  geom_hline(yintercept = 0, lty='dashed', col = 'grey')+
  geom_point(data=ibs.meta, aes(x=MDS1, y=MDS2, col=AdmixGroup), size =4, alpha =.4) + theme_classic() + 
  scale_color_viridis_d(option = "turbo") +
  geom_text_repel(data=cent_label, aes(x=centx, y=centy, label=AdmixGroup), col = 'gray4') + 
  theme(legend.position = "none", text=element_text(size = 10)) + #coord_fixed() +
  ylim(c(-1.5, .5)) + xlim(c(-.7,.7))

ggsave("figures/ADMIX5.PCA.PLOS.pdf", height=4, width = 4, units = "in")


virdpal <- viridis::viridis(n=5, option="turbo")
names(virdpal) <- c("V1", "V2", "V3", "V4", "V5")
virdpal = data.frame(virdpal)
virdpal$AdmixGroup = rownames(virdpal)

labcol = ibs.meta %>% dplyr::select(AdmixGroup, label) %>% left_join(virdpal)
#ibs.meta.col <- ibs.meta %>% mutate(mclust=as.numeric(mclust)) %>% left_join(ibs.meta, tipdf, by=c("mclust"="group"))
library(dendextend)
dend <- as.dendrogram(hclust(as.dist(ibs), method="ave"))

pdf("figures/PLOS_unrooted_samps.pdf", width = 6, height = 6 )
plot(as.phylo(dend), type = "unrooted", tip.color = labcol[labcol$label == labels(as.phylo(dend)),]$virdpal, cex = .4)
dev.off()


#######################################
###### get pie charts by loc ##########
#######################################
max_ad <- max_ad %>% mutate(Site=ifelse(Site=="sOrpheus2", "sOrpheus", Site))
ad_pie <- max_ad %>% group_by(Site, AdmixGroup) %>% summarise(count_gs=n()) %>% ungroup() %>% group_by(Site) %>% mutate(pr_gs=count_gs/sum(count_gs))
ggplot(ad_pie, aes(x="", y=pr_gs, fill=AdmixGroup)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + facet_wrap(~Site) + scale_fill_viridis_d(option="turbo") + theme_void() + 
  theme(legend.position = "none", text=element_blank())
ggsave("figures/MEREVIEW.piesadmix.png",width = 2, height = 2, units = "in", dpi=500, bg='transparent')



########################################
####### INTEGRATE ZOOX WITH AD GROUPS ###
########################################

its2.hell <- decostand(data.frame(its2.filtered.colony), method = "hellinger")
its2.hell.dist <- vegdist(its2.hell, method = "bray")
its2.hell.cap <- capscale(its2.hell.dist~1)
plot(its2.hell.cap, choices=c(1,2))
pheatmap(1-its2.hell.dist)


its2.total <- data.frame(decostand(data.frame(its2.filtered.colony), method = "total"))
its2.total$ID <- rownames(its2.total)
its2.total <- its2.total %>% pivot_longer(starts_with("ASV"), names_to = "ASV", values_to = "prop") %>% 
  separate(ASV, into=c("asv", "group", "ASV_group")) %>%
  left_join(consensus_groups)
its2.total <- its2.total %>% filter(prop > 0)
its2.total <- its2.total %>% left_join(meta, by=c("ID"="ColonyID")) %>%
  dplyr::select(ID, ASV_group, prop, consensus_call, Site) %>% unique()
its2.total <- left_join(its2.total, ibs.meta, by=c("ID"="ColonyID"))# %>% left_join(meta, by=c("ID"="ColonyID"))
order_by_site <- its2.total %>% arrange(Site.x) 
group2consensus <- its2.total %>% dplyr::select(ASV_group, consensus_call) %>% unique() %>% arrange(ASV_group)
table(group2consensus$consensus_call)
library(colorspace)
C15_cols <-  sequential_hcl(n=27, palette = "Mint")
c15.tmp <- group2consensus %>% filter(consensus_call == "C_15") %>% mutate(color=C15_cols)
D4_cols <- "#c13c9b"
d4.tmp <- group2consensus %>% filter(consensus_call == "D_4") %>% mutate(color=D4_cols)
unknown_cols <- sequential_hcl(n=14, palette = "Lajolla")
unk.tmp <- group2consensus %>% filter(consensus_call == "NA_NA" | is.na(consensus_call == T)) %>%
  mutate(color=unknown_cols[1:10])

getasv_ord <- its2.total %>% mutate(ID=factor(ID, levels=unique(order_by_site$ID))) %>%
  mutate(ASV_group=as.factor(ASV_group)) %>%
  mutate(ASV_group=reorder(ASV_group, prop))
  
#make a stacked barplot of zoox profiles according to admixture group
col_order <- rbind(c15.tmp, d4.tmp, unk.tmp) 
col_order$ASV_group <- factor(col_order$ASV_group, levels=levels(getasv_ord$ASV_group))
col_order <- col_order %>% arrange(ASV_group)
#doesn't matter how it's ordered
its2.total %>% mutate(ID=factor(ID, levels=unique(order_by_site$ID)), Site.y=ifelse(Site.y=="sOrpheus2", "sOrpheus", Site.y)) %>%
  mutate(ASV_group=as.factor(ASV_group)) %>%
  mutate(ASV_group=reorder(ASV_group, prop)) %>%
  mutate(AdmixGroup=ifelse(is.na(AdmixGroup), "Unknown", AdmixGroup)) %>%
  mutate(AdmixGroup=gsub("V", "Admix_", AdmixGroup)) %>%
  mutate(AdmixGroup=factor(AdmixGroup, levels=c("Admix_1", "Admix_2", "Admix_3", "Admix_4", "Admix_5", "Unknown"))) %>%
  ggplot(aes(x=ID, y=prop, fill=ASV_group)) + geom_bar(stat="identity", width = 1, col = 'white', lwd = .075) + theme_classic() +
  scale_fill_manual(values=col_order$color) +
  ylab("Proportion") + xlab("Sample") + 
  theme(legend.position = "none", text=element_text(size = 8), axis.text.x=element_blank()) + facet_wrap(~AdmixGroup, scales = "free_x", nrow = 1) 
  
#ggsave("figures/ITS2_site_admix_sep.png", width = 6.5, height = 2, units = "in", dpi = 500)

png("figures/ITS2.legend.png", res = 500)
par(mfrow=c(3,1), mar=c(1,1,1,1))
image(1:nrow(unk.tmp), 1, as.matrix(1:nrow(unk.tmp)), 
             col=unknown_cols,
             xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
image(1:nrow(c15.tmp), 1, as.matrix(1:nrow(c15.tmp)), 
      col=C15_cols,
      xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
image(1:nrow(d4.tmp), 1, as.matrix(1:nrow(d4.tmp)), 
      col=D4_cols,
      xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
dev.off()
#keep_meta <- meta_big %>% dplyr::select(maxgroup, Site.x, Age.x, ID.x) %>% drop_na()


########################################
########################################
#### NOW DO ZOOX ANALYSIS ##############

bycolony.meta = left_join(bycolony.meta, ibs.meta[,c("ColonyID", "AdmixGroup")]) %>%
  mutate(AdmixGroup = ifelse(is.na(AdmixGroup), NA, AdmixGroup))

#normalize differently
its2.hell <- decostand(data.frame(its2.filtered.colony), method = "hellinger")
its2.hell.dist <- vegdist(its2.hell, method = "bray")
its2.hell.cap <- capscale(its2.hell.dist~1)
plot(its2.hell.cap, choices=c(1,2))
pheatmap(1-its2.hell.dist)

logfold_its2abun = log10(t(data.frame(its2.filtered.colony))+1)
its2.group = decostand(data.frame(t(its2.filtered.colony)), method = "hellinger")
its2.group.dist = vegdist(its2.group, method = "bray")

#get better ordering
library(dendsort)
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
mat_cluster_rows <- sort_hclust(hclust(its2.group.dist))
group2name = group2consensus %>% mutate(longname=paste(consensus_call, " (AG_", ASV_group, ")",  sep = ""))
group2name[is.na(group2name)] = "No Consensus"
newnames = group2name$longname
names(newnames) = paste("ASV_group_", group2name$ASV_group, sep = "")

rownames(logfold_its2abun) = newnames[rownames(logfold_its2abun)]
getpheatinfo = bycolony.meta
rownames(getpheatinfo) = getpheatinfo$ColonyID
getpheatinfo = getpheatinfo[colnames(logfold_its2abun),]

getpheatinfo = getpheatinfo %>% arrange(AdmixGroup)
bks = getpheatinfo %>% group_by(AdmixGroup) %>% summarise(tot=n())

anot = data.frame(AdmixGroup = getpheatinfo[,c("AdmixGroup")])
rownames(anot) = getpheatinfo[,c("ColonyID")]
anot$AdmixGroup = gsub("V", "Admix_", anot$AdmixGroup)


mycolors <- c( turbo((5)), "white")
names(mycolors) <- unique(anot$AdmixGroup)
mycolors <- list(AdmixGroup = mycolors)

dev.off()
png("figures/PLOS_ITS2.heatmap.png", width = 6.5, height = 4, units = "in", res = 1000)
pheatmap(logfold_its2abun[,getpheatinfo$ColonyID],
         annotation_col = anot, cluster_cols = T,
         #gaps_col = cumsum(bks$tot),
         color = colorRampPalette(mako(100))(100), 
         annotation_colors = mycolors, fontsize = 4,
         legend = T,
         border_color = NA,
         angle_col = "45", lwd = .8)

dev.off()

allscore.its2 <- scores(its2.hell.cap, scaling="species", choices =c(1:5), tidy = T) %>% filter(score=="sites") %>% 
 left_join(bycolony.meta, by=c("label"="ColonyID"))
allscore.its2 %>% ggplot(aes(x=MDS1, y=MDS2, col=AdmixGroup)) + 
  geom_vline(xintercept = 0, lty='dashed', col ='gray') +
  geom_hline(yintercept = 0, lty='dashed', col ='gray') +
  geom_point(size=3, alpha = .7) + theme_classic() + 
  xlab("MDS1")+ylab("MDS2") + scale_color_viridis_d(option="turbo") + coord_equal() +
  guides(color=guide_legend(title="AdmixGroup")) +
  theme(text=element_text(size = 10), legend.position = "none")
#ggsave("figures/ITSwithsitelegend.png", width = .66*6.5, height = 3, units = "in", dpi = 500)
ggsave("figures/MEREV_ITS2_AdmixCol.png", height = 3, units = "in", dpi = 500)


#mds1, mds3
allscore.its2 %>% ggplot(aes(x=MDS1, y=MDS3, col=as.factor(AdmixGroup))) + 
  geom_vline(xintercept = 0, lty='dashed', col ='gray') +
  geom_hline(yintercept = 0, lty='dashed', col ='gray') +
  geom_point(size=3, alpha = .7) + theme_classic() + 
  xlab("MDS1")+ylab("MDS3") + scale_color_viridis_d(option="viridis") + coord_equal() +
  guides(color=guide_legend(title="Adgroup")) +
  theme(text=element_text(size = 10), legend.position = "none")

#keepsamps
keepITS2 <- allscore.its2 %>% drop_na()
keepITS2 <- keepITS2$label
its2.df.gen <- data.frame(its2.filtered.colony)
its2.df.gen <- its2.df.gen[keepITS2,]
its2.gen.meta <- allscore.its2 %>% filter(label %in% keepITS2)
its2.hell.gen <- data.frame(decostand(its2.df.gen, method = "hellinger"))
sampleDists.its2.gen <- vegdist(its2.hell.gen, method = "bray")
its2.shared.cap <- capscale(sampleDists.its2.gen~1)
plot(its2.shared.cap) 
sharedits.pca <- scores(its2.shared.cap, tidy = T, choices=c(1:3), scaling = "species") %>% filter(score=="sites") %>% left_join(bycolony.meta, by=c("label"="ColonyID"))

sharedits.pca %>% ggplot(aes(x=MDS1, y=MDS2, col=AdmixGroup)) + 
  geom_vline(xintercept = 0, lty='dashed', col ='gray') +
  geom_hline(yintercept = 0, lty='dashed', col ='gray') +
  geom_point(size=3, alpha = .7) + theme_classic() + 
  xlab("MDS1")+ylab("MDS2") + scale_color_viridis_d(option="turbo") + coord_equal() +
  guides(color=guide_legend(title="AdmixGroup")) +
  theme(text=element_text(size = 10), legend.position = "none")
#ggsave("figures/MEREV_ITS2_SHARED_AdmixCol.png", height = 3, units = "in", dpi = 500)

adonis2(sampleDists.its2.gen~as.factor(its2.gen.meta$AdmixGroup) + its2.gen.meta$Site+its2.gen.meta$Age, method = "margin")

library(gradientForest) #must be loaded after RF 
library(fastDummies)
source("scripts/R/RDA-forest_functions.R")

get_col <- rownames(as.matrix(ibs))
get_col <- meta %>% filter(ID %in% get_col) %>% dplyr::select(ColonyID)
ibs.ma <- as.matrix(ibs)
rownames(ibs.ma) <- get_col$ColonyID
colnames(ibs.ma) <- get_col$ColonyID
shared <- intersect(rownames(ibs.ma), rownames(as.matrix(its2.hell.dist)))
zoox.pro <- capscale(as.matrix(its2.hell.dist)[shared, shared]~1)
gen.pro <- capscale(ibs.ma[shared,shared]~1)

protest(zoox.pro, gen.pro) #not significant
plot(protest(zoox.pro, gen.pro))

#now subset your meta data to just the predictors we care about. 
zoox_preds <- its2.gen.meta %>% dplyr::select(Site, Age, AdmixGroup)

#dummify
zoox_covs<- dummy_cols(zoox_preds, select_columns = c('Site', 'Age', 'AdmixGroup'), remove_selected_columns = T)
zoox_covs <- zoox_covs[,colSums(zoox_covs)>0 & colSums(zoox_covs) < nrow(zoox_covs)]

plot(its2.shared.cap)
#Run gradientForest.
zoox.gf <- makeGF(its2.shared.cap, zoox_covs, ntrees = 1500)
importance(zoox.gf)
env=zoox_covs
plot(zoox.gf)

ii=data.frame(importance(zoox.gf))
names(ii)="importance"
# reordering by decreasing importance
ii$var=factor(row.names(ii),levels=rev(row.names(ii)))
#Pretty Plot :) 
ggplot(ii[1:min(10,nrow(ii)),],aes(var,importance))+
  geom_bar(stat="identity")+
  theme_bw()+
  coord_flip()+
  xlab("")

importance.cutoff=0
env=zoox_preds
is=c();ii0=ii[ii$importance>=importance.cutoff,]
for (f in colnames(env)) {
  i2=ii0[grep(paste("^",f,"_",sep=""),ii0$var),]
  if (nrow(i2)==0) {i2=ii0[grep(paste("^",f,sep=""),ii0$var),] }
  is=append(is,sum(i2$importance))
}
ii2=data.frame(cbind(importance=is,var=colnames(env)))
#ii2[ii2$var == "mclust",]$var <- "SubCluster"
ii2$importance=as.numeric(ii2$importance)
ii2$var=factor(ii2$var,levels=ii2$var[order(ii2$importance)])
#levels(ii2$var) <- c(levels(ii2$var)[1:2], "GeneticGroup")
zooxgf <- ggplot(ii2[1:min(100,nrow(ii2)),],aes(var,importance))+geom_bar(stat="identity")+theme_classic()+
  theme(text = element_text(size = 10)) +coord_flip()+xlab("") + ylab("Summed Variable Importance\n(Weighted)") 
zooxgf
sum(ii2$importance) #admix group = .175, about the same
zooximp <- ii2
zooximp$relimp <- zooximp$importance/sum(zooximp$importance)
#ggsave("figures/zooxgf10_euc.png", zooxgf, width = 3, height = 2, dpi = 500, units = "in")


########################################
########################################
#### MICROBE ASVS  ##############
load("data/Raegens_16S_ASVs.Rdata")
raw_microbe <- data.frame(t(as.matrix(data.frame(poritesASVs))))
raw_microbe$ID <- rownames(raw_microbe)
raw_microbe <- raw_microbe %>% pivot_longer(!ID, names_to = "ASV_group", values_to="readcount")
raw_microbe <- raw_microbe %>% left_join(meta) %>% group_by(ColonyID, ASV_group) %>% summarise(summed_count=sum(readcount)) %>%
  pivot_wider(names_from=ASV_group, values_from=summed_count) %>% drop_na()
raw_microbe <- data.frame(raw_microbe)
rownames(raw_microbe) <- raw_microbe$ColonyID
raw_microbe$ColonyID <- NULL
#filter by presense absense in colonies....
mic.pa <- raw_microbe
mic.pa[mic.pa>0]=1
hist(log10(rowSums(mic.pa)))
keeprows <- rownames(raw_microbe)[rowSums(mic.pa) >= 3] #must be present in at least 3 samps
keepcols <- colnames(raw_microbe)[colSums(raw_microbe) >=100] #trash samples with fewer than 100 total ASV counts
keepcols <- intersect(keepcols, colnames(raw_microbe)[colSums(mic.pa) > 1]) #gotta be present in more than one sample
keeprows <- intersect(keeprows, rownames(raw_microbe)[rowSums(raw_microbe) > 1000]) #gotta be present 

raw_mic.filtered <- raw_microbe[keeprows,keepcols] #got rid of 50 entries
range(rowSums(raw_mic.filtered))
mic.hell <- decostand(raw_mic.filtered, method = "hellinger")
mic.bray <- vegdist(mic.hell, method = "bray")
mic.cap <- capscale(mic.bray~1) #looks nice :) 
plot(mic.cap)
allscore.mic <- scores(mic.cap, scaling="species", choices =c(1:5), tidy = T) %>% filter(score=="sites") %>% 
  left_join(bycolony.meta, by=c("label"="ColonyID")) %>% 
  mutate(Site=ifelse(grepl("sOrpheus2", Site), "sOrpheus", Site))

allscore.mic %>% ggplot(aes(x=MDS1, y=MDS2, col=AdmixGroup)) + 
  geom_vline(xintercept = 0, lty='dashed', col ='gray') +
  geom_hline(yintercept = 0, lty='dashed', col ='gray') +
  geom_point(size=3, alpha = .7) + theme_classic() + 
  xlab("MDS1")+ylab("MDS2") + scale_color_viridis_d(option="turbo") + coord_equal() +
  guides(color=guide_legend(title="AdmixGroup")) +
  theme(text=element_text(size = 10), legend.position = "none")

get_col <- rownames(as.matrix(ibs))
get_col <- meta %>% filter(ID %in% get_col) %>% dplyr::select(ColonyID)
ibs.ma <- as.matrix(ibs)
rownames(ibs.ma) <- get_col$ColonyID
colnames(ibs.ma) <- get_col$ColonyID
shared <- intersect(rownames(ibs.ma), rownames(as.matrix(mic.hell)))
mic.pro <- capscale(as.matrix(mic.bray)[shared, shared]~1)
gen.pro <- capscale(ibs.ma[shared,shared]~1)
protest(mic.pro, gen.pro) #not significantly correlated in procrustes

#now subset your meta data to just the predictors we care about. 
getcols <- allscore.mic %>% drop_na() %>% dplyr::select(Site, Age, AdmixGroup, label)

mic_preds <- allscore.mic %>% drop_na() %>% dplyr::select(Site, Age, AdmixGroup)

#dummify
mic_covs<- dummy_cols(mic_preds, select_columns = c('Site', 'Age', 'AdmixGroup'), remove_selected_columns = T)
mic_covs <- mic_covs[,colSums(mic_covs)>0 & colSums(mic_covs) < nrow(mic_covs)]

mic.shared.cap <- capscale(as.matrix(mic.bray)[getcols$label, getcols$label]~1)
plot(mic.shared.cap)
sharedmic.pca <- scores(mic.shared.cap, tidy = T, choices=c(1:3), scaling = "species") %>% filter(score=="sites") %>% left_join(bycolony.meta, by=c("label"="ColonyID"))

sharedmic.pca %>% ggplot(aes(x=MDS1, y=MDS2, col=AdmixGroup)) + 
  geom_vline(xintercept = 0, lty='dashed', col ='gray') +
  geom_hline(yintercept = 0, lty='dashed', col ='gray') +
  geom_point(size=3, alpha = .7) + theme_classic() + 
  xlab("MDS1")+ylab("MDS2") + scale_color_viridis_d(option="turbo") + coord_equal() +
  guides(color=guide_legend(title="AdmixGroup")) +
  theme(text=element_text(size = 10), legend.position = "none")

adonis2(as.matrix(mic.bray)[getcols$label, getcols$label]~as.factor(mic_preds$AdmixGroup)+mic_preds$Site+mic_preds$Age, by = "margin")

#Run gradientForest.
mic.gf <- makeGF(mic.shared.cap, mic_covs, ntrees = 1500) #something is wrong with the split?
importance(mic.gf)
env=mic_covs

plot(mic.gf)

###yes! Do this at the start, just need to know colony ID, then calculate distances
#I've removed clones from the ibs matrices, so this is effectively the same thing.
#and we can boxplot within colony vs. within admix group centroid dists. 


ii=data.frame(importance(mic.gf))
names(ii)="importance"
# reordering by decreasing importance
ii$var=factor(row.names(ii),levels=rev(row.names(ii)))
#Pretty Plot :) 
ggplot(ii[1:min(10,nrow(ii)),],aes(var,importance))+
  geom_bar(stat="identity")+
  theme_bw()+
  coord_flip()+
  xlab("")

importance.cutoff=0
env=mic_preds
is=c();ii0=ii[ii$importance>=importance.cutoff,]
for (f in colnames(env)) {
  i2=ii0[grep(paste("^",f,"_",sep=""),ii0$var),]
  if (nrow(i2)==0) {i2=ii0[grep(paste("^",f,sep=""),ii0$var),] }
  is=append(is,sum(i2$importance))
}
ii2=data.frame(cbind(importance=is,var=colnames(env)))
ii2$importance=as.numeric(ii2$importance)
ii2$var=factor(ii2$var,levels=ii2$var[order(ii2$importance)])
micgf <- ggplot(ii2[1:min(100,nrow(ii2)),],aes(var,importance))+geom_bar(stat="identity")+theme_classic()+
  theme(text = element_text(size = 10)) +coord_flip()+xlab("") + ylab("Summed Variable Importance\n(Weighted)") 
micgf
micimp <- ii2 

micimp$relimp <- micimp$importance/sum(micimp$importance)

micimp$Community <- "Microbe"
zooximp$Community <- "Symbio."
allimp <- rbind(micimp, zooximp)
allimp$var <- as.character(allimp$var)
allimp[allimp$var == "Age",]$var <- "Size"
allimp[allimp$var == "AdmixGroup",]$var <- "AdmixGroup"
allimp$var <- factor(allimp$var, levels=c("AdmixGroup", "Site", "Size"))
com_imp <- allimp %>% drop_na() %>% ggplot(aes(x=var,y=relimp, fill=Community)) + geom_bar(stat="identity", width=.5, position = "dodge") +
  scale_fill_manual(values=c("forestgreen", "gold3")) + theme_classic() + xlab("") + ylab("Weighted Relative Importance") + 
  theme(text=element_text(size = 10),
        legend.position = c(.7, .85),
        #legend.position = "top",
        legend.key.height= unit(.1, 'in'),
        legend.key.width= unit(.1, 'in'),
        axis.text.x = element_text(angle = 35, vjust = .65, hjust=.5)) + labs(fill="")
com_imp


########################################################
###### Get microbe abundance plot and heatmap ##########
########################################################
load("data/asstax.Rdata")
ass_microbes <- data.frame(out)
ass_microbes$ASV <- rownames(ass_microbes)
mic_toplot <- data.frame(t(as.matrix(data.frame(poritesASVs))))
mic_toplot$ID <- rownames(mic_toplot)
mic_toplot <- mic_toplot %>% pivot_longer(!ID, names_to = "ASV", values_to = "count") %>%
  left_join(ass_microbes, by=c("ASV"="ASV"))
mic_toplot <- mic_toplot %>% filter(count>10) %>% group_by(ID, Class) %>% summarise(by_fam_count=sum(count))
mic_toplot <- mic_toplot %>% ungroup() %>% group_by(ID) %>% mutate(perctax=by_fam_count/sum(by_fam_count))
mic_toplot <- mic_toplot %>% filter(perctax > .05 & perctax < .9) %>% drop_na()
mic_toplot <- mic_toplot %>% group_by(ID) %>% mutate(perctax=by_fam_count/sum(by_fam_count)) 
mic_toplot <- mic_toplot %>% left_join(meta, by=c("ID"="ID")) 
order_by_site <- mic_toplot %>% ungroup() %>% arrange(Site)
getasv_ord <- mic_toplot %>% mutate(ID=factor(ID, levels=unique(order_by_site$ID))) %>%
  mutate(Class=as.factor(Class)) %>%
  mutate(Class=reorder(Class, perctax))

mic_toplot <- raw_mic.filtered
mic_toplot$ID <- rownames(mic_toplot)
mic_toplot <- mic_toplot %>% pivot_longer(!ID, names_to = "ASV", values_to = "count") %>%
  left_join(ass_microbes, by=c("ASV"="ASV"))
mic_toplot <- mic_toplot %>% filter(count>10) %>% group_by(ID, Class) %>% summarise(by_fam_count=sum(count))
mic_toplot <- mic_toplot %>% ungroup() %>% group_by(ID) %>% mutate(perctax=by_fam_count/sum(by_fam_count))
mic_toplot <- mic_toplot %>% filter(perctax > .05 & perctax < .9) %>% drop_na()
mic_toplot <- mic_toplot %>% group_by(ID) %>% mutate(perctax=by_fam_count/sum(by_fam_count)) 
mic_toplot <- mic_toplot %>% left_join(bycolony.meta, by=c("ID"="ColonyID")) 

order_by_site <- mic_toplot %>% ungroup() %>% arrange(Site)
getasv_ord <- mic_toplot %>% mutate(ID=factor(ID, levels=unique(order_by_site$ID))) %>%
  mutate(Class=as.factor(Class)) %>%
  mutate(Class=reorder(Class, perctax))

#this is a horrible plot, almost any other data vis would be better 
getasv_ord %>% mutate(Site=ifelse(Site=="sOrpheus2", "sOrpheus", Site))%>% 
  mutate(AdmixGroup=ifelse(is.na(AdmixGroup), "Unknown", AdmixGroup)) %>%
  mutate(AdmixGroup=gsub("V", "Admix_", AdmixGroup)) %>%
  mutate(AdmixGroup=factor(AdmixGroup, levels=c("Admix_1", "Admix_2", "Admix_3", "Admix_4", "Admix_5", "Unknown"))) %>%
  ggplot(aes(x=ID,y=perctax, fill=Class)) +
  geom_bar(stat="identity", width = 1, col = 'white', lwd = .09) + theme_classic() +
  theme(legend.position = "bottom") + scale_fill_viridis_d(option="turbo", begin = 0, direction = 1) +facet_wrap(~AdmixGroup, scales="free_x", nrow=1) +
  ylab("Proportion") +xlab("Sample") + theme(axis.text.x=element_blank(), text=element_text(size=8)) +
  guides(fill=guide_legend(nrow=3))

ggsave("figures/MicTaxa_StackedBar_Legend.png", width = 8, height = 3, units = "in", dpi = 500)


#### Pheatmap For Microbes (keep me)
p.long = data.frame(poritesASVs) %>% mutate(ASV = rownames(data.frame(poritesASVs))) %>% pivot_longer(!ASV, names_to = "ID", values_to = "count")
p.wide = p.long %>% left_join(meta) %>% 
  group_by(ColonyID, ASV) %>% summarise(count = mean(count)) %>%
  pivot_wider(names_from = ColonyID, values_from = count) %>% data.frame()
rownames(p.wide) = p.wide$ASV
p.wide$ASV = NULL
p.stand = decostand(data.frame(p.wide), "total", MARGIN =2)
keeps.p = head(sort(rowSums(p.stand), decreasing = T), 30)
toplot.microbes = data.frame(p.wide)[names(keeps.p),]
toplot.microbes$ASV = rownames(toplot.microbes)
toplot.microbes = toplot.microbes %>% left_join(ass_microbes)
toplot.microbes$r = 1:nrow(toplot.microbes)
table(toplot.microbes$Genus)
toplot.microbes = toplot.microbes %>% mutate(lab.part1 = ifelse(is.na(Species), paste("Unk_", r, sep = ""), 
                                  Species))
rownames(toplot.microbes) = paste(toplot.microbes$Genus, toplot.microbes$lab.part1)
toplot.microbes = toplot.microbes %>% dplyr::select(epa1:soj9)
#### do heatmap instead of stacked barplot
logfold_mic2abun = log10(toplot.microbes+1)

getpheatinfo = bycolony.meta
rownames(getpheatinfo) = getpheatinfo$ColonyID
getpheatinfo = getpheatinfo[colnames(logfold_mic2abun),]

getpheatinfo = getpheatinfo %>% arrange(AdmixGroup)
bks = getpheatinfo %>% group_by(AdmixGroup) %>% summarise(tot=n())

anot = data.frame(AdmixGroup=getpheatinfo[,c("AdmixGroup")])
rownames(anot) = getpheatinfo$ColonyID
anot$AdmixGroup = gsub("V", "Admix_", anot$AdmixGroup)

anot.cols = list(AdmixGroup = c(Admix_1 = turbo(5)[1],
                                Admix_2 = turbo(5)[2],
                                Admix_3 = turbo(5)[3],
                                Admix_4 = turbo(5)[4],
                                Admix_5 = turbo(5)[5],
                                "NA" = "white"))

dev.off()
png("figures/PLOS_MICOTU.heatmap.png", width = 6.5, height = 4, units = "in", res = 1000)
pheatmap(logfold_mic2abun[,getpheatinfo$ColonyID],
         border_color = NA,
         legend =T,
         annotation_col = anot, cluster_cols = T,
         gaps_col = cumsum(bks$tot),
         color = colorRampPalette(mako(100))(100), 
         annotation_colors = anot.cols, fontsize=4,
         angle_col = "45")

dev.off()


#############################################
### WITHIN VS. BETWEEN GROUPS
#check betgween within 
tmp.meta <- meta %>% left_join(max_ad, by=c("ColonyID")) %>% data.frame()
rownames(tmp.meta) <- tmp.meta$ID.x
get_samps <- tmp.meta %>% filter(Age.x=="adult") %>% rownames()

out.pa <- outdf
out.pa[out.pa>0]=1
hist(log10(rowSums(out.pa)))
keeprows <- rownames(outdf)[rowSums(out.pa) >= 3] #must be present in at least 3 samps
keepcols <- colnames(outdf)[colSums(outdf) >=100] #trash samples with fewer than 100 total ASV counts
keepcols <- intersect(keepcols, colnames(outdf)[colSums(out.pa) > 1])
keeprows <- intersect(keeprows, rownames(outdf)[rowSums(outdf) > 100])

its2.filtered <- outdf[keeprows,keepcols] #got rid of 50 entries
get_samps <- intersect(colnames(its2.filtered), get_samps)

byindivid.its2 <- decostand(data.frame(t(its2.filtered))[get_samps,], method = "hellinger")
byindivid.its2.dist <- vegdist(byindivid.its2, method = "bray", na.rm=T)
by_colony <- as.factor(tmp.meta[get_samps,]$ColonyID)
by_group <- as.factor(tmp.meta[get_samps,]$AdmixGroup)
by_site <- as.factor(tmp.meta[get_samps,]$Site.x)


bycol.disper <- betadisper(byindivid.its2.dist, by_colony)
bygroup.disper <- betadisper(byindivid.its2.dist, by_group)
bysite.disper <- betadisper(byindivid.its2.dist, by_site)
#bytr.disper <- betadisper(byindivid.its2.dist, by_tr[get_samps])


plot(bygroup.disper)
boxplot(bycol.disper$distances, bygroup.disper$distances, bysite.disper$distances, names=c("Within Colony", "Within Group", "Within Site"))
to.test <- rbind(data.frame(distance=bycol.disper$distances, group="Within Colony"),
                 data.frame(distance=bysite.disper$distances, group="Within Site"),
                 data.frame(distance=bygroup.disper$distances, group="Within Group"))
test.out <- pairwise.t.test(to.test$distance, to.test$group)

library(ggpubr)
library(rstatix)

pwc <- to.test %>% 
  pairwise_wilcox_test(
    distance ~ group, pool.sd = F,
    p.adjust.method = "holm", paired = F
  )
pwc
pwc <- pwc %>% add_xy_position(x = "group", step.increase = 1)
to.test$group <- factor(to.test$group, levels = c("Within Colony", "Within Group", "Within Site"))
ggboxplot(to.test, x = "group", y = "distance") +
  stat_pvalue_manual(pwc, hide.ns = TRUE, bracket.size = .2, bracket.nudge.y = -.5, tip.length = .02, step.increase = -.15) + ylab("Distance to Group Centroid") + xlab("")  +
  theme(text=element_text(size = 8))
ggsave("figures/ITS2WithinBetweenColony.png", width = 3, height = 3, units = "in", dpi = 500)



#############################################
### WITHIN VS. BETWEEN GROUPS MICROBES
#check betgween within 
tmp.meta <- meta %>% left_join(max_ad, by=c("ColonyID")) %>% data.frame()
rownames(tmp.meta) <- tmp.meta$ID.x
get_samps <- tmp.meta %>% filter(Age.x=="adult") %>% rownames()

tmpdf <- data.frame(as.matrix(data.frame(poritesASVs)))
out.pa <- tmpdf
out.pa[out.pa>0]=1
hist(log10(rowSums(out.pa)))
keeprows <- rownames(tmpdf)[rowSums(out.pa) >= 3] #must be present in at least 3 samps
keepcols <- colnames(tmpdf)[colSums(tmpdf) >=100] #trash samples with fewer than 100 total ASV counts
keepcols <- intersect(keepcols, colnames(tmpdf)[colSums(out.pa) > 1])
keeprows <- intersect(keeprows, rownames(tmpdf)[rowSums(tmpdf) > 100])

its2.filtered <- tmpdf[keeprows,keepcols] #got rid of 50 entries
get_samps <- intersect(colnames(its2.filtered), get_samps)

byindivid.its2 <- decostand(data.frame(t(its2.filtered))[get_samps,], method = "hellinger")
byindivid.its2.dist <- vegdist(byindivid.its2, method = "bray", na.rm=T)
by_colony <- as.factor(tmp.meta[get_samps,]$ColonyID)
by_group <- as.factor(tmp.meta[get_samps,]$AdmixGroup)
by_site <- as.factor(tmp.meta[get_samps,]$Site.x)


bycol.disper <- betadisper(byindivid.its2.dist, by_colony)
bygroup.disper <- betadisper(byindivid.its2.dist, by_group)
bysite.disper <- betadisper(byindivid.its2.dist, by_site)


plot(bygroup.disper)
plot(bysite.disper)
boxplot(bycol.disper$distances, bygroup.disper$distances, bysite.disper$distances, names=c("Within Colony", "Within Group", "Within Site"))
to.test <- rbind(data.frame(distance=bycol.disper$distances, group="Within Colony"),
                 data.frame(distance=bysite.disper$distances, group="Within Site"),
                 data.frame(distance=bygroup.disper$distances, group="Within Group"))
test.out <- pairwise.t.test(to.test$distance, to.test$group)


pwc <- to.test %>% 
  pairwise_wilcox_test(
    distance ~ group, pool.sd = F,
    p.adjust.method = "holm", paired = F, alternative="less"
  )
pwc
pwc <- pwc %>% add_xy_position(x = "group", step.increase = 1)
to.test$group <- factor(to.test$group, levels = c("Within Colony", "Within Group", "Within Site"))
ggboxplot(to.test, x = "group", y = "distance") +
  stat_pvalue_manual(pwc, hide.ns = TRUE, bracket.size = .2, bracket.nudge.y = 0, tip.length = .02, step.increase = -.05) + ylab("Distance to Group Centroid") + xlab("")  +
  theme(text=element_text(size = 8))

ggsave("figures/MICROBEWithinBetweenColony.png", width = 3, height = 3, units = "in", dpi = 500)

