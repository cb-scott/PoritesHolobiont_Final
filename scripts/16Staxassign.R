#I RAN THIS ON THE REMOTE CLUSTER - IT WILL TAKE TOO LONG ON YOUR LOCAL COMPUTER
library(dada2)
load("Raegens_16S_ASVs.Rdata")
out = assignTaxonomy(rownames(poritesASVs), "silvaref.fa.gz")
save(out, file = "asstax.Rdata")

### LOAD REMOTE RUN TO LOCAL 
load("data/asstax.Rdata")
#out is the taxa
ass_taxa = out

load("data/Raegens_16S_ASVs.Rdata")
load("data/MetadataMicrobial.Rdata") #gen ad is the metadata
raw_microbe <- as.matrix(data.frame(poritesASVs))
library(DESeq2)
library(vegan)
library(tidyverse)
library(phyloseq)
micmeta <- retain_samps_mic %>% select(ID, ColonyID, Site.x, Age.x, maxgroup) %>% unique()
micmeta <- data.frame(micmeta)
rownames(micmeta) = micmeta$ID
#genetically unique data 
raw2 = raw_microbe[colnames(raw_microbe) %in% rownames(micmeta),]
p <- phyloseq(otu_table(raw2, taxa_are_rows = T), tax_table(ass_taxa), sample_data=sample_data(micmeta))
diagdds = phyloseq_to_deseq2(p, ~ maxgroup + Site.x + Age.x)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric", sfType = "poscounts")

#get the basemean
base = phyloseq_to_deseq2(p, ~1)
base <- estimateSizeFactors(base, type = 'poscounts')
vsd.base <- varianceStabilizingTransformation(base, blind = T)
count.base = assay(vsd.base)
count.base[count.base <0] = 0 #vst counts <0 correspond to counts less than 1. Approximate as 0
count.base = data.frame(count.base)
count.base$ASV = rownames(count.base)
ass_taxa.df = data.frame(ass_taxa)
ass_taxa.df$ASV = rownames(ass_taxa.df)
count.taxa <- left_join(count.base, ass_taxa.df)
count.taxa <- count.taxa %>% pivot_longer(cols = E1:X96, names_to = "ID", values_to = "counts")
#genus level
count.family <- count.taxa %>% group_by(ID, Genus) %>% summarize(samp_gen_counts = sum(counts))
count.all <- count.family %>% filter(samp_gen_counts >= 5) %>% drop_na() %>% left_join(micmeta)
count.all$ID <- as.factor(count.all$ID)
count.all$order <- paste(count.all$maxgroup, count.all$Site.x)
count.all$ID <- reorder(count.all$ID, as.numeric(as.factor(count.all$order)))
a <- count.all %>% ggplot(aes(fill = Genus, y=samp_gen_counts, x=as.numeric(ID))) + geom_bar(position="fill", stat="identity")
selline <- count.all %>% select(ID, maxgroup) %>% unique() #will give you a vertical line at each genetic grouping
vline = cumsum(as.vector(table(selline$maxgroup)))
#vline = cumsum(as.vector(table(micmeta$Site.x)))
a2 <- a +  geom_vline(xintercept = vline[1:length(vline)-1]+.5, lwd = 1) + theme_classic() + 
  theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank()) + scale_fill_viridis_d(option = "turbo")  + xlab("Sample") + ylab("Relative Abundance") + 
  guides(fill=guide_legend(ncol=8))
a2
ggsave(filename  = "figures/mic_genus.png", a2, units = "in", dpi = 500, width = 6.5, height = 5)
#summarize by genus

######### SIGNIFICANT TAXA ###################

setseed(6)
getsig <- function(f, l1, l2){
  res = results(diagdds, contrast=c(f, l1, l2), cooksCutoff = FALSE)
  alpha = 0.01
  sigtab1 = res[which(res$padj < alpha), ]
  sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(p)[rownames(sigtab1), ], "matrix"))
  sigtab1$comparison = paste(l1, "vs", l2, sep = "_")
  return(sigtab1)
}
site1 <- getsig("Site.x", "nePelorus", "sePelorus")
site2 <- getsig("Site.x", "nePelorus", "sOrpheus")
site3 <- getsig("Site.x", "nePelorus", "pioneerBay")
site4 <- getsig("Site.x", "sOrpheus", "sePelorus")
site5 <- getsig("Site.x", "sOrpheus", "pioneerBay")
site6 <- getsig("Site.x", "pioneerBay", "sePelorus")
age <- getsig("Age.x", "adult", "juvenile")
#nosig dif site4? check on your end
allres <- rbind(site1, site2, site3, site5, site6, age)


# Genus order
x = tapply(allres$log2FoldChange, allres$Genus, function(x) max(x))
x = sort(x, TRUE)
allres$Genus = factor(as.character(allres$Genus), levels=names(x))
foldchange <- ggplot(allres, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3, alpha = .8) + theme_bw() + 
  ylim(c(-50, 60)) +
  geom_hline(yintercept = 0, col = 'grey', lty = 'dashed')+
  facet_grid(rows = vars(comparison)) +
  xlab("") + theme(strip.text = element_text(
    size = 6), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size = 8), legend.text = element_text(size = 8), legend.position = "right")
foldchange
ggsave(filename= "figures/allcomp_lfc.png", foldchange, width = 5, height = 9.5, units = "in", dpi = 500)

####################################
#### Get that nice heatmap :) ########
#########################
library(viridis)
set.seed(2)
rturbo <- sample(turbo(n = 10))
names(rturbo) <- c(1:10)

library("pheatmap")
ntd <-normTransform(diagdds)
select <- order(rowMeans(counts(diagdds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(diagdds)[,c( "Age.x","Site.x",  "maxgroup")])
df <- arrange(df, maxgroup, Site.x, Age.x)
colnames(df) <- c("Age", "Site", "GenGroup")
mydat <- data.frame(assay(ntd)[select,rownames(df)])
mydat$ASV <- rownames(mydat)
mydat = left_join(mydat, ass_taxa.df)
mydat$rowname = paste(mydat$Family, mydat$Species, sep = ",")
vird <- viridis(4)
bpal <- brewer.pal(4, "Accent")
names(bpal) <- c("nePelorus", "pioneerBay", "sePelorus", "sOrpheus")
ann_colors = list(
  GenGroup = rturbo,
  Site = bpal,
  Age = c(adult = "darkgrey", juvenile = "azure2")
 # CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
 # GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
)
png("figures/best20microbes.png", width = 6.5, height = 4.5, units = "in", res = 500)
#Family and Genus Names
#pheatmap(assay(ntd)[select,rownames(df)], cluster_rows=T, show_rownames=T, show_colnames = F,
 #        cluster_cols=F, annotation_col=df, labels_row = mydat$rowname, scale = "none",
  #       annotation_colors = ann_colors, gaps_col = cumsum(as.numeric(table(df$GenGroup))), fontsize = 8)
#Final = just family names
pheatmap(assay(ntd)[select,rownames(df)], cluster_rows=T, show_rownames=T, show_colnames = F,
        cluster_cols=F, annotation_col=df, labels_row = mydat$Family, scale = "none",
      annotation_colors = ann_colors, gaps_col = cumsum(as.numeric(table(df$GenGroup))), fontsize = 8)

dev.off()
