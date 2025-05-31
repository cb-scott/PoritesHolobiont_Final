######################################################
######### Create new porites genome #################
######################################################
#Porites lobata
#https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02960-7#availability-of-data-and-materials
cd /scratch/06909/cbscott/new_plob_genome

export PGEN=/scratch/06909/cbscott/new_plob_genome/Porites_lobata_v3.fa
cat $PGEN /work/06909/cbscott/sym_references/symABCD_longref.fasta > P_lobata_zoox.fasta

conda activate map_software

echo "bowtie2-build P_lobata_zoox.fasta P_lobata_zoox.fasta" > btb 
python3 $HOME/bin/ls6_launcher_creator.py -j btb -n btb -a IBN21018 -e cbscott@utexas.edu -t 02:00:00 -N 1 -w 1 -q development
sbatch btb.slurm
############ TRIM AND MAP READS #############
conda activate map_software 
#Trim fastqs
#These files exist in: /work/06909/cbscott/PoritesProject/reads/PoritesHostReads/trimmed_fastqs
>trimse
for file in *.fq; do
echo "cutadapt -q 15,15 -m 25 -o ${file/.fq/}.trim $file > ${file}_trimlog.txt" >> trimse;
done
python3 $HOME/bin/ls6_launcher_creator.py -j trimse -n trimse -a IBN21018 -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 24 -q development

#map them to the gen we've already built
#mapped to P. Lutea Genome
#export GEN=/scratch/06909/cbscott/new_plob_genome/P_lobata_zoox.fasta
export GEN=/work/06909/cbscott/PoritesProject/genome/porties_sym_genome.fasta

>map_porites
for F in /work/06909/cbscott/PoritesProject/reads/PoritesHostReads/trimmed_fastqs/*.trim; do
  echo "bowtie2 --no-unal --local -x $GEN -U $F -S ${F#/work/06909/cbscott/PoritesProject/reads/PoritesHostReads/trimmed_fastqs/}.master.local.sam && \
  samtools sort -O bam -o ${F#/work/06909/cbscott/PoritesProject/reads/PoritesHostReads/trimmed_fastqs/}.master.local.sorted.bam ${F#/work/06909/cbscott/PoritesProject/reads/PoritesHostReads/trimmed_fastqs/}.master.local.sam && samtools index -c ${F#/work/06909/cbscott/PoritesProject/reads/PoritesHostReads/trimmed_fastqs/}.master.local.sorted.bam" >> map_porites; done
python3 $HOME/bin/ls6_launcher_creator.py -j map_porites -n map_porites -a IBN21018 -e cbscott@utexas.edu -t 02:00:00 -N 1 -w 24 -q development
sbatch map_porites.slurm

#coral scaffolds exists in work dir
cp /work/06909/cbscott/PoritesProject/genome/coral_scaffolds .
#coral scaffolds needs carrot removed
sed -i 's/>//' coral_scaffolds

#Get rid of zooxanthellae reads
>coral_sep
for file in *.sorted.bam; do
echo "cat coral_scaffolds | tr '\n' ' ' | xargs samtools view -bh $file > ${file/.sorted.bam}.coral.sorted.bam && samtools index -c ${file/.sorted.bam}.coral.sorted.bam" >> coral_sep; done
python3 $HOME/bin/ls6_launcher_creator.py -j coral_sep -n coral_sep -a IBN21018 -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 24 -q development
sbatch coral_sep.slurm

###################################################
##### REMOVE LOW QUALITY BAMS/LOW COV SITES########
###################################################

conda activate angsd
# quality assessment, removing bams with log(coverage)<3SD
# also imposng minimum number of individuals(MI) a locus must be seen in (genotyping rate cutoff - 50%)
export MinIndPerc=0.5
FILTERSQ='-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minInd $MI'
TODOQ="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
echo 'export NIND=`cat bams | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc1
echo "ls *coral*.bam > bams && source calc1 && angsd -b bams -GL 1 $FILTERSQ $TODOQ -P 12 -out dd">a0
python3 $HOME/bin/ls6_launcher_creator.py -j a0 -n a0 -a IBN21018 -e cbscott@utexas.edu -t 2:00:00 -w 1 -q development 
  sbatch --dependency=afterok:2366516 a0.slurm 

conda activate R_TACC
Rscript ~/bin/plotQC.R prefix=dd

#get rid of poor quality samples, in addition to species conflicts, rerun
wc -l quality.txt #83
#retain samples with at least 10% sites at 5 coverage (why not)
tail -79 quality.txt | awk '{print $1}' > goodqual.bams #bads is a file I made based on the first hclust tree


###################################################
##### CLUSTER AND REMOVE CLONES ###################
###################################################
conda activate angsd

export MinIndPerc=0.5 #cut off looks ok
# initial IBS production, detecting and removing clones (see hctree.pdf and resulting bams.nr)
FILTERS0='-minInd $MI -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -doHWE 1 -hetbias_pval 1e-3 -skipTriallelic 1'
TODO0='-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2'
echo 'export NIND=`cat goodqual.bams | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc1
echo "source calc1 && angsd -b goodqual.bams -GL 1 $FILTERS0 $TODO0 -P 12 -out result_rem_lq" >a1
python3 $HOME/bin/ls6_launcher_creator.py -j a1 -n a1 -a IBN21018 -e cbscott@utexas.edu -t 2:00:00 -w 1 -q development

#actually, let's pull these guys locally and filter on my end
cp result_rem_lq.ibsMat all_highq_Porites.ibsMat
awk -F. '{print $1}' goodqual.bams > all_highq_Porites.names

conda activate R_TACC
#set cut off at .10 - the level of true duplicates

Rscript ~/bin/detect_clones.R goodqual.bams result_rem_lq.ibsMat 0.10


###################################################
##### FINAL IBS MAT GENERATION  ###################
###################################################

#we've removed clonal individuals based on our tecnical and colony replicates.
#additionally, we've removed all low quality and low coverage BAMS
#total: retain 67 bams
conda activate angsd

export MinIndPerc=0.5 #cut off looks ok
# initial IBS production, detecting and removing clones (see hctree.pdf and resulting bams.nr)
FILTERS0='-minInd $MI -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -doHWE 1 -hetbias_pval 1e-3 -skipTriallelic 1'
TODO0='-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2 -doBcf 1'
echo 'export NIND=`cat bams.nr | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc1
echo "source calc1 && angsd -b bams.nr -GL 1 $FILTERS0 $TODO0 -P 12 -out clone_qual_filteredPorites_" >a2
python3 $HOME/bin/ls6_launcher_creator.py -j a2 -n a2 -a IBN21018 -e cbscott@utexas.edu -t 1:00:00 -w 1 -q development
sbatch a2.slurm


###################################################
##### RUN ADMIXTURE ON RESULTING FILE #############
###################################################
conda activate pcangsd
>check_admix
for i in {2..12}; do 
echo "pcangsd -b clone_qual_filteredPorites_.beagle.gz --admix_K $i --admix -o clone_qual_filteredPorites_newgen_${i} --tree" >> check_admix; done 
python3 $HOME/bin/ls6_launcher_creator.py -j check_admix -n check_admix -a IBN21018 -e cbscott@utexas.edu -t 2:00:00 -w 12 -q normal
sbatch check_admix.slurm

awk -F. '{print $1}' bams.nr > clone_qual_filteredPorites.names

pcangsd -b clone_qual_filteredPorites_.beagle.gz --admix --admix_iter 1000 --admix_batch 15 --admix_tol 1e-10 -o clone_qual_filteredPorites_newgen_strict --tree

#pcangsd -b clone_qual_filteredPorites.beagle.gz --admix --admix_K 4 -o clone_qual_filter4

####now do angsd again for every admixture group.

clone_qual_filteredPorites.bamlist.V1.admix.beagle.gz
conda activate angsd
>all.angsd
for file in bamlist*.admix; do #created these locally in R, need them for FST
export MinIndPerc=0.5 #cut off looks ok
# initial IBS production, detecting and removing clones (see hctree.pdf and resulting bams.nr)
FILTERS0='-minInd $MI -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -snp_pval 1e-5 -minMaf 0.05 -dosnpstat 1 -doHWE 1 -hetbias_pval 1e-3 -skipTriallelic 1' &&\
TODO0='-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2' &&\
echo 'export NIND=`cat $file | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' > $file.calc1 &&\
echo "source $file.calc1 && angsd -b $file -GL 1 $FILTERS0 $TODO0 -P 12 -out clone_qual_filteredPorites.${file}" >> all.angsd; done
python3 $HOME/bin/ls6_launcher_creator.py -j all.angsd -n all.angsd -a IBN21018 -e cbscott@utexas.edu -t 2:00:00 -w 4 -q normal
sbatch all.angsd.slurm

conda activate pcangsd
for file in clone_qual_filteredPorites.bamlist.*.admix.beagle.gz; do
  pcangsd -b $file --admix -o ${file/.beage.gz}.ADMIXAGAIN; done


###################################################################
###### GET PAIRWISE FST BETWEEN ADMIXTURE GROUPS ##################
###################################################################
conda activate angsd 


export POR_GEN=/work/06909/cbscott/PoritesProject/genome/porties_sym_genome.fasta

#filter out poor quality reads and calculate the .idx files
>dofst_ad
for file in bamlist*admix; do
echo "angsd -b $file -anc $POR_GEN -out ${file/.bams} -dosaf 1 -gl 1 -minMapQ 30 -minQ 20" >> dofst_ad; done
python3 $HOME/bin/ls6_launcher_creator.py -j dofst_ad -n dofst_ad -a IBN21018 -e cbscott@utexas.edu -t 1:00:00 -w 10 -q development


#I think I need these pairwise
>do_sfs
for i in {1..5}; do
  for j in {1..5}; do
  echo "realSFS bamlist.V$i.admix.saf.idx bamlist.V$j.admix.saf.idx > ad2ad.${i}_${j}.ml.small" >> do_sfs; done; done

python3 $HOME/bin/ls6_launcher_creator.py -j do_sfs -n do_sfs -a IBN21018 -e cbscott@utexas.edu -t 12:00:00 -w 12 -N 1 -q vm-small

>do_fst
for i in {1..5}; do
  for j in {1..5}; do
  echo "realSFS fst index bamlist.V$i.admix.saf.idx bamlist.V$j.admix.saf.idx -sfs ad2ad.${i}_${j}.ml.small -fstout ad2ad.${i}_${j}.fst" >> do_fst; done; done
#grep -v "$(cat remove)" do_fst > redFST
python3 $HOME/bin/ls6_launcher_creator.py -j do_fst -n do_fst -a IBN21018 -e cbscott@utexas.edu -t 04:00:00 -w 12 -N 1 -q vm-small

>calcfst
for file in *.fst.fst.idx; do
  echo "realSFS fst stats $file > ${file/.fst.fst.idx}.fst_stat" >>calcfst ; done
python3 $HOME/bin/ls6_launcher_creator.py -j calcfst -n calcfst -a IBN21018 -e cbscott@utexas.edu -t 1:00:00 -w 12 -N 1 -q vm-small

ls *.fst_stat >allfst

>Merge.fst.admix5
for file in `cat allfst`; do
  cat $file >> Merge.fst.admix5; done

paste allfst Merge.fst.admix5 > named.fst.admix5 #then pull local into R

