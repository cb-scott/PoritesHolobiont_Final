library(tidyverse)
library(RColorBrewer)
library(vegan)
library(randomForest)
library(caret)
library(gradientForest) #must be loaded after RF!!!! Masks functions!!!!
library(fastDummies)
library(mclust)
source("scripts/RDAforest/R/RDA-forest_functions.R")



#Read in your IBS matrix and data 
#names of your files
nam = read.table("data/clone_qual_filteredPorites.names")
colnames(nam) = "ID"

#Read in the metadata 
meta <- read.csv("data/Porites_Meta.csv")
#Clean it up a little bit
meta <- meta %>% mutate(Site=ifelse(Site=="sOrpheus2", "sOrpheus", Site))
meta <- meta %>% select(SampleID:Orientation)
meta <- meta[!apply(meta == "", 1, all),]
colID = left_join(nam, meta, by=c("ID" = "Edited.ID"))


#Load in your genetic distance Rdata object (or create it)
#all this does is load an object called gendist, which is a matrix of genetic
#distances with rownames and columnames corresponding to sample ID. 
#THIS IS THE SAME AS THE ANGSD IBS MATRIX - but saving myself the hassle of renaming 
#rows and columns every time I read it in.
load("data/gendist_convenient.Rdata")
#use bayesian clustering to find the optimal number of clusters
d_clust <- Mclust(as.matrix(gendist), G=1:15, 
                  modelNames = mclust.options("emModelNames"))

#Gen_ad actually contains the ADMIXTURE classifications for Ncluster = 4 (if you're intersted in this)
#It's not informative and does a bad job clustering because of the cryptic substructure
genad2 = cbind(gen_ad,d_clust$classification)

#I'm immediately going to trash this info - it's not informative and clutters the environment/analysis
gen_ad = genad2 %>% select(ID:`d_clust$classification`)
gen_ad$maxgroup <- gen_ad$`d_clust$classification` #rename this, bad convention
gen_ad$maxgroup <- as.factor(gen_ad$maxgroup) #genetic cluster should be a factor


#########################################################
##### RANDOM FOREST: HOST GENETICS ######################
#########################################################
set.seed(4)
train_group = floor(.7*nrow(gen_ad)) #grab 70% data as a train set
grab_rows <- sample(nrow(gen_ad), train_group)
train <- gen_ad[grab_rows,]
test <- gen_ad[!rownames(gen_ad)  %in% grab_rows,]

train <- train %>% drop_na() #drop NAs
gen.rf <- randomForest(train[,c("Site", "Age")], train$maxgroup, importance = T)
varImpPlot(gen.rf) #check your variable importances, accuracy is more intuitive than GINI to me, use that
gen.pred <- predict(gen.rf, test)

confmat = confusionMatrix(gen.pred, test$maxgroup) 
confmat #look at your accuracy, pval, etc. 
gen.imp <- data.frame(importance(gen.rf))
genrf.plot <- gen.imp %>% select(MeanDecreaseAccuracy, MeanDecreaseGini) %>% mutate(group = rownames(gen.imp)) %>% 
  ggplot(aes(x=group,y=MeanDecreaseAccuracy)) + 
  geom_bar(stat="identity") + theme_classic()+
  theme(text = element_text(size = 10)) +coord_flip()+xlab("") +
  geom_text(data = data.frame(), aes(label = paste("Accuracy:", round(confmat$overall[1], 2), "\np=", round(confmat$overall[6], 3)), x = 1, y = 12), size = 2)# + ylab("Summed Variable Importance\n(Weighted)") 
genrf.plot

#ggsave("figures/genrf.png", genrf.plot, dpi = 500, width = 2.5, height = 1.5, units = "in")


#########################################################
########### PCA/RDA: HOST GENETICS ######################
#########################################################

#Assign a color pallette for the genetic groups
library(viridis)
set.seed(2)
rturbo <- sample(turbo(n = 10))
names(rturbo) <-  as.character(c(1:10))

gen.cap <- capscale(gendist~1) #use vegan to make a PCA (~1 is an unconstrained ordination)
gen.score <- scores(gen.cap, choices = c(1:5), tidy=T) %>% filter(score == "sites")
gen.meta <- left_join(gen.score, gen_ad, by=c("label" = "ID"))

####PC1 vs PC2
gen.plot <- gen.meta %>% ggplot(aes(x=MDS1, y=MDS2, col = maxgroup)) + geom_vline(xintercept = 0, lty = 'dashed', col = 'gray') + 
  geom_hline(yintercept = 0, lty = 'dashed', col = 'gray') +
  geom_point(size = 3, alpha = .7) + coord_equal() + theme_classic() + scale_color_manual(values = rturbo) + 
  theme(text=element_text(size = 10), legend.position = "none")
gen.plot
#ggsave("figures/genplot10.png", gen.plot, units = "in", dpi = 500, height = 3, width = 3)

#####PC2 vs PC3
gen.plot <- gen.meta %>% ggplot(aes(x=MDS3, y=MDS2, col = maxgroup)) + geom_vline(xintercept = 0, lty = 'dashed', col = 'gray') + 
  geom_hline(yintercept = 0, lty = 'dashed', col = 'gray') +
  geom_point(size = 3, alpha = .7) + coord_equal() + theme_classic() + scale_color_manual(values = rturbo) + 
  theme(text=element_text(size = 10), legend.position = "none")
gen.plot
#ggsave("figures/genplotMDS32_10.png", gen.plot, units = "in", dpi = 500, height = 3, width = 4)

adonis2(gendist~Site+Age, data = gen_ad)


#######################################
##### SYMBOINTS: GF and RDA ###########
#######################################

load("data/its2_dist_euc.RData")
vst_dist = count.dist
meta = read.csv("~/Box Sync/Projects/Por_Bucket_AmpSeq/data/Updated_Porites_Meta.csv", header = T)
meta[meta$Site == "sOrpheus2",]$Site = "sOrpheus"
#subset the metadata for samples that are in vstdist 
rownames(gen_ad) <- gen_ad$ID
zoox_metasub <- meta[meta$ID %in% rownames(vst_dist),]
rownames(zoox_metasub) <- zoox_metasub$ID

#assign info by colony id (only one genotype point per colony but potentially 3 zoox replicates)
retain_samps_zoox <- left_join(zoox_metasub, gen_ad, by = c("ColonyID"="ColonyID")) %>% drop_na()
#check it. 
head(retain_samps_zoox)

#double check your data frame to make sure all columns exist 
zoox_dist <- vst_dist[rownames(vst_dist) %in% retain_samps_zoox$ID.x, colnames(vst_dist) %in% retain_samps_zoox$ID.x]
#now subset your meta data to just the predictors we care about. 
zoox_preds <- retain_samps_zoox %>% select(Site.x, SiteType.x, Age.x, maxgroup) %>% rename(Site=Site.x, SiteType=SiteType.x, Age=Age.x)
#dummify
zoox_covs<- dummy_cols(zoox_preds, select_columns = c('Site', 'SiteType', 'Age', 'maxgroup'), remove_selected_columns = T)
#check it
head(zoox_covs)

#Create a base ordination object for zoox community dist (a PCA)
zoox.cap <- capscale(zoox_dist ~1)
zoox.score <- scores(zoox.cap, choices = c(1:3), tidy = T) %>% filter(score == "sites")
zoox.pca <- cbind(zoox.score, zoox_preds) %>% ggplot(aes(x=MDS1, y=MDS2, col = maxgroup)) + geom_vline(xintercept = 0, lty = 'dashed', col = 'gray') + 
  geom_hline(yintercept = 0, lty = 'dashed', col = 'gray') +
  geom_point(size = 3, alpha = .7) + coord_fixed() + theme_classic() + scale_color_manual(values = rturbo) + 
  theme(text=element_text(size = 10), legend.position = "none")
zoox.pca

#procrustes analysis to see how much raw zoox PCA aligns with genetic PCA
shared <- intersect(rownames(gendist), rownames(zoox_dist))
zoox.pro <- zoox_dist[shared, shared]
gen.pro <- gendist[shared,shared]
plot(protest(zoox.pro, gen.pro))
protest(zoox.pro, gen.pro)

#Run Gradient Forest
zoox.gf <- makeGF(zoox.cap, zoox_covs, ntrees = 1500)
env=zoox_covs


ii=data.frame(importance(zoox.gf))
names(ii)="importance"
# reordering by decreasing importance
ii$var=factor(row.names(ii),levels=rev(row.names(ii)))
#Pretty Plot to check results :) 
ggplot(ii[1:min(10,nrow(ii)),],aes(var,importance))+
  geom_bar(stat="identity")+
  theme_bw()+
  coord_flip()+
  xlab("")

#Sum importances over factor levles 
importance.cutoff=0
env=zoox_preds[,c(1,3,4)]
is=c();ii0=ii[ii$importance>=importance.cutoff,]
for (f in colnames(env)) {
  i2=ii0[grep(paste("^",f,"_",sep=""),ii0$var),]
  if (nrow(i2)==0) {i2=ii0[grep(paste("^",f,sep=""),ii0$var),] }
  is=append(is,sum(i2$importance))
}

ii2=data.frame(cbind(importance=is,var=colnames(env)))
ii2$importance=as.numeric(ii2$importance)
ii2$var=factor(ii2$var,levels=ii2$var[order(ii2$importance)])
levels(ii2$var) <- c(levels(ii2$var)[1:2], "GeneticGroup")
#plot summed over levels
zooxgf <- ggplot(ii2[1:min(100,nrow(ii2)),],aes(var,importance))+geom_bar(stat="identity")+theme_classic()+
  theme(text = element_text(size = 10)) +coord_flip()+xlab("") + ylab("Summed Variable Importance\n(Weighted)") 
zooxgf
#get relative importances and save for later
zooximp <- ii2
zooximp$relimp <- zooximp$importance/sum(zooximp$importance)


#do the RDA.
zoox.rda <- capscale(zoox_dist~maxgroup+Site+Age, data = zoox_preds)
plot(zoox.rda)
#PERMANOVA:
adonis2(zoox_dist~maxgroup+Site+Age, data = zoox_preds) #get significance

######################################
##### MICROBES: GF and RDA ###########
######################################
#This will follow effectively the same procedure for the zoox,
#but I don't have an "easy object" saved, so we'll do the DESeq standardization here. 

load("data/Raegens_16S_ASVs.Rdata")
raw_microbe <- as.matrix(data.frame(poritesASVs))
library(DESeq2)

#create a deseq object with no formula
microbe.dds = DESeqDataSetFromMatrix(countData=raw_microbe, colData=data.frame(colnames(raw_microbe)),design= ~1)


startdim = dim(counts(microbe.dds))
keeprows <- rowSums(counts(microbe.dds)) >= 10 #remove any ASV's that don't have at least 10 reads mapping to them across all samples
keepcols <- colSums(counts(microbe.dds)) >=50 #trash samples with fewer than 50 total read count
filter = counts(microbe.dds)
filter[filter != 0] <- 1
keepcols <- colSums(filter) >= 5 #get rid of ASVS that don't appear in at least 5 samples

microbe.dds <- microbe.dds[keeprows,keepcols]
enddim = dim(counts(microbe.dds))
paste("Your starting number of ASVS was", startdim[1], ". Your ending number of ASVS is", enddim[1])
paste("Your starting number of samples was", startdim[2], ". Your ending number of samples is", enddim[2])

#Get stabilized counts as recommended by DESeq's documentation and waste not want not
microbe.dds <- estimateSizeFactors(microbe.dds, type = 'poscounts')
vsd <- varianceStabilizingTransformation(microbe.dds, blind = T)
count.tab = t(assay(vsd))
count.tab[count.tab <0] = 0 #vst counts <0 correspond to counts less than 1 in the stabilized object. Approximate as 0
microbe.dist <- vegdist(count.tab, method = "euclidean") # I think this is the best way, because now our zeroes are meaningful
microbe.dist <- as.matrix(microbe.dist)
#look at PCA
microbe.dds.cap <- capscale(microbe.dist~1)
plot(microbe.dds.cap)

###RUN A PROCRUSTES TEST ON IT
#how well does this align with our genetic distances?
shared <- intersect(rownames(gendist), rownames(microbe.dist))
mic.pro <- microbe.dist[shared, shared]
gen.pro <- gendist[shared,shared]
protest(mic.pro, gen.pro)

#Subset the metadata for samples that exist
microbe_metasub <- meta[meta$ID %in% rownames(microbe.dist),]
rownames(microbe_metasub) <- microbe_metasub$ID
microbe_metasub <- microbe_metasub[rownames(microbe.dist),]
retain_samps_mic <- left_join(microbe_metasub, gen_ad, by = c("ColonyID"="ColonyID")) %>% drop_na()
#output this to easily input to other scripts down the line
#save(retain_samps_mic, file = "data/MetadataMicrobial.Rdata")


#Subset your distance matrix for the samples you want to retain
mic_dist <- microbe.dist[rownames(microbe.dist) %in% retain_samps_mic$ID.x, colnames(microbe.dist) %in% retain_samps_mic$ID.x]
#now subset your meta data to just the predictors we care about. 
mic_preds <- retain_samps_mic %>% select(Site.x, SiteType.x, Age.x, maxgroup) %>% dplyr::rename(Site=Site.x, SiteType=SiteType.x, Age=Age.x)
#dummify
mic_covs<- dummy_cols(mic_preds, select_columns = c('Site', 'Age', 'maxgroup'), remove_selected_columns = T) %>% select(!SiteType)

mic.cap <- capscale(mic_dist ~1)

#Look at the colored PCA
mic.score <- scores(mic.cap, choices = c(1:3), tidy = T) %>% filter(score == "sites")
mic.pca <- cbind(mic.score, mic_preds) %>% ggplot(aes(x=MDS1, y=MDS2, col = maxgroup)) + geom_vline(xintercept = 0, lty = 'dashed', col = 'gray') + 
  geom_hline(yintercept = 0, lty = 'dashed', col = 'gray') +
  geom_point(size = 3, alpha = .7) + coord_fixed() + theme_classic() + scale_color_manual(values = rturbo) + 
  theme(text=element_text(size = 10), legend.position = "none") +ylim(c(-4,6)) + xlim(-5,4)
mic.pca
#ggsave("figures/micpca_10.png", mic.pca, dpi = 500, units = "in", height = 2.5, width = 2.25)

#Run Gradient Forest
env=mic_covs
mic.gf <- makeGF(mic.cap, mic_covs, ntrees = 1500)
plot(mic.gf)


#GET var importances - THIS WILL OVERWRITE THE ZOOX OBJECT 
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


#Combine by factor level 
importance.cutoff=0
env=mic_preds[,c(1,3,4)]
is=c();ii0=ii[ii$importance>=importance.cutoff,]
for (f in colnames(env)) {
  i2=ii0[grep(paste("^",f,"_",sep=""),ii0$var),]
  if (nrow(i2)==0) {i2=ii0[grep(paste("^",f,sep=""),ii0$var),] }
  is=append(is,sum(i2$importance))
}
ii2=data.frame(cbind(importance=is,var=colnames(env)))
ii2$importance=as.numeric(ii2$importance)
ii2$var=factor(ii2$var,levels=ii2$var[order(ii2$importance)])
ggplot(ii2[1:min(100,nrow(ii2)),],aes(var,importance))+geom_bar(stat="identity")+theme_classic()+
  theme(text = element_text(size = 20)) +coord_flip()+xlab("") + ylab("Summed Variable Importance-Weighted") 

#Combine with Zoox Variable Importances
levels(ii2$var) <- c(levels(ii2$var)[1:2], "GeneticGroup")
ii2$type = "Microbe"
zooximp$type = "Symbio."
ii2$relimp <- ii2$importance/sum(ii2$importance)
zooximp$relimp <- zooximp$importance/sum(zooximp$importance)
all <- data.frame(rbind(ii2, zooximp)) #doesn't work because I edited zooximp
mic_imp <- ii2
all$var <- ordered(all$var, c("GeneticGroup", "Site", "Age"))
all$type <- ordered(all$type, c("Symbio.", "Microbe"))
com_imp <- all %>% ggplot(aes(x=var,y=importance, fill=type)) + geom_bar(stat="identity", width=.5, position = "dodge") +
  scale_fill_manual(values=c("forestgreen", "gold3")) + theme_classic() + xlab("") + ylab("Weighted Importance") + 
  theme(text=element_text(size = 10),
        legend.position = c(.7, .85),
        #legend.position = "top",
        legend.key.height= unit(.1, 'in'),
        legend.key.width= unit(.1, 'in'),
        axis.text.x = element_text(angle = 35, vjust = .65, hjust=.5)) + labs(fill="")
com_imp
#ggsave("figures/compareGF.png", com_imp, width = 2, height = 2.75, units = "in", dpi = 500)

### DO RDA
mic.rda <- capscale(mic_dist~Age+Site+maxgroup, data = mic_preds)
plot(mic.rda)
#PERMANOVA
adonis2(mic_dist~maxgroup+Site+Age, data = mic_preds)#more overall variation in microbial community 


