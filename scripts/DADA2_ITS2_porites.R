library(dada2)
library(ShortRead)
library(Biostrings)

#check for primer seqs 
path <- "data/POR_ITS2/"
list.files(path)

fnFs <- sort(list.files(path, pattern = "1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "2.fastq.gz", full.names = TRUE))


#used cudadapt to remove primers already, so skip over

# Forward and reverse fastq filenames have the format:
cutFs <-fnFs
cutRs <- fnRs

#ITS2
FWD="GAATTGCAGAACTCCGTGAACC"
REV="CGGGTTCWCTTGTYTGACTTCATGC"

#Illumina
#FWD="CTACACGACGCTCTTCCGATCT"
#REV="CAGACGTGTGCTCTTCCGATCT"


allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "-")[[1]][2]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])


fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#check for primers
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
## YOU SHOULD HAVE ALL ZEROES HERE!
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#primers have been successfully trimmed....and illumina adapters....



path.cut = "data/POR_ITS2/filtN/"
cutFs <- sort(list.files(path.cut, pattern = "1.fastq.gz", full.names = T))
cutRs <- sort(list.files(path.cut, pattern = "2.fastq.gz", full.names = T))


filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
#takes forever!
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 5), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set mult

#CAN START HERE
path.filts = "data/POR_ITS2/filtN/filtered/"
filtFs <- sort(list.files(path.filts, pattern = "1.fastq.gz", full.names = T))
filtRs <- sort(list.files(path.filts, pattern = "2.fastq.gz", full.names = T))


#these all take forever - definitely needs to be on TACC?
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = T)

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
sample.names <- unname(sapply(filtFs, get.sample.name)) 

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#Considerations for your own data: Most of your reads should successfully merge.
#If that is not the case upstream parameters may need to be revisited: 
#Did you trim away the overlap between your reads? Potentially!!!
#Might have primers at start and end of read inflating sequences... may want to pass cutadapt again?
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, justConcatenate = T) #just concatenate doesn't consider paired end 
#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, justConcatenate = F) #just concatenate doesn't consider paired end 

#I don't think my reads overlap - after trimming primers read length only 126, target region is longer
#for more detail:
#https://github.com/benjjneb/dada2/issues/790
#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, minOverlap = 2)

seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
#write.table(colnames(seqtab.nochim), "zoox_asvs.table", quote = F, row.names = F, col.names = F)
seqout <- data.frame(seqtab.nochim)
seqout$SampleID <- rownames(seqtab.nochim)
saveRDS(seqout, file = "~/Box Sync/Projects/Porites Project/data/UpdatedSeqTabITS2.Rdata")
#write.table(seqout, "~/Box Sync/Projects/Porites Project/data/UpdatedSeqTabITS2.tab",  quote = F, row.names = F, col.names = T)

#grab a database from Danielle Claar
#https://github.com/baumlab/Claar_etal_2020_SciRep/blob/master/ITS2db_trimmed_derep_dada.fasta

its2.ref <- "data/SymPortal_unique_DIVs.fasta"  # CHANGE ME to location on your machine
taxa <- assignSpecies(seqtab.nochim, its2.ref, tryRC = TRUE)


taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

tp <- data.frame(taxa.print)
table(tp$Genus)
unique(tp$Species)

###########################################
## what if we don't merge? just classify F and R
# 
# seqtab <- makeSequenceTable(dadaRs)
# table(nchar(getSequences(seqtab)))
# 
# seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# dim(seqtab.nochim)
# sum(seqtab.nochim)/sum(seqtab)
# 
# #grab a database from Danielle Claar
# #https://github.com/baumlab/Claar_etal_2020_SciRep/blob/master/ITS2db_trimmed_derep_dada.fasta
# 
# its2.ref <- "data/ITS2db_trimmed_derep_dada.fasta"  # CHANGE ME to location on your machine
# taxa <- assignTaxonomy(seqtab.nochim, its2.ref, multithread = TRUE, tryRC = TRUE)
# 
# 
# taxa.print <- taxa  # Removing sequence rownames for display only
# rownames(taxa.print) <- NULL
# head(taxa.print)
# 
# tp <- data.frame(taxa.print)
# tp$Genus
# table(tp$Species)
# table(tp$Genus)
# 
# #we get better resolution?!
# 


##########################################
#### GO WITH THE MERGED PAIRS FOR NOW#####
##########################################
#see: https://www.biostars.org/p/9498358/

#BiocManager::install("phyloseq")
library(phyloseq)
library(tidyverse)

samples.out <- rownames(seqtab.nochim)
#did I do replication in my experiment:
length(samples.out)
length(unique(samples.out))
#yes, not all files names are unique

#get your metadata - this processing will be custom for your data format!
meta <- read_csv("~/Box Sync/Projects/Porites Project/data/Porites_Meta.csv", col_names = T)
meta <- meta %>% select(`Edited ID`:Orientation)
meta <- meta %>% rename(`Edited ID`="ID")

#rename my duplicated samples
sam_names = data.frame(samples.out)
colnames(sam_names) = c("ID")
#duplicate metadata rows for duplicated samples
meta <- left_join(sam_names, meta)
new_ids <- c()
for(i in 1:length(meta$ID)){
  new_ids <- c(new_ids, ifelse(meta$ID[i] %in% new_ids, paste(meta$ID[i], "_1", sep = ""), meta$ID[i]))
}
meta$ID <- new_ids
#do it twice to get rid of triplicate values ^^

rownames(seqtab.nochim) <- new_ids
rownames(meta) <- new_ids
write.csv(meta, "data/Updated_Porites_Meta.csv", quote = F)
saveRDS(seqtab.nochim, file = "data/ITS2_ASVS.rds")
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa))

saveRDS(ps, file = "data/ITS2_ps.rds")


#Convert DNA strings to more readable names 
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
write.csv(otu_table(ps), "data/namedASVtable_ITS2.csv")
plot_richness(ps, x="Age", measures=c("Shannon", "Simpson"), color="Age")

ps.prop <-transform_sample_counts(ps, function(otu) otu/sum(otu))

#get the relative abundance of otus by sample using ps.prop
abun_mat <- as.matrix(na.omit(otu_table(ps.prop)))
library(vegan)
dist_mat <- vegdist(abun_mat, method = "bray")
pheatmap(dist_mat)

abun.cap <- capscale(abun_mat~1, dist = "bray")
plot(abun.cap)

score_dist <- as.data.frame(abun.cap$CA$u)
score_dist$ID = rownames(score_dist)
#str(score_dist)
#str(meta)
score_meta <- left_join(score_dist, meta)

ggplot(score_meta, aes(x=MDS1, y =MDS2, col = as.factor(Age))) + geom_point(size = 3, alpha = .8) + theme_classic() +
  coord_equal() + geom_hline(yintercept = 0, col = 'grey', lty = 'dashed') + geom_vline(col = 'grey', xintercept = 0, lty = 'dashed')

#remove Na's from meta:
meta_nona = meta[rownames(meta) %in% rownames(abun_mat),]

#age significant?
abun.cap <- capscale(abun_mat~Age, meta_nona, dist = "bray")
plot(abun.cap)

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]

#ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

#DESeq with this?
write.csv(abun_mat, file = "~/Box Sync/Projects/Porites Project/data/ITS2_ASV_Matrix.csv", quote = F)
