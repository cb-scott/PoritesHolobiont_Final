######## REMOVE CORAL CONTAMINATION IN READS #########
#This only needs to be done for 16S reads, not ITS2

#Use MetaTAXA to blast and classify
conda create -n metaxa # you need to create a new conda environment, this software has coflicts with your previous environment
conda activate metaxa  #activate new enviornment
conda install -c bioconda blast #install blast
#conda install -c bioconda blast-legacy #need the legacy version of blast
conda install hmmer
conda install -c bioconda mafft
conda install -c bioconda metaxa #install metataxa
conda install seqtk


#note: I tested this on Jumana's reads becasue there are less files, so be sure to
#check your particular file paths!

cd BUCKET #change directories into where your reads are. For Jumana, this should just be the 16S reads directory
ls
#make sure you see a bunch of files named *.1.fastq.gz and *.2.fastq.gz

#need to unzip your reads. This will take a while
>unzipper
for file in *.fastq.gz; do
echo "gunzip $file" >> unzipper; done

ls6_launcher_creator.py -j unzipper -n unzipper -a IBN21018 -e YOUREMAIL@utexas.edu -t 00:30:00 -N 1 -w 24 -q vm-small
sbatch unzipper.slurm

#Wait for this to run.

#let's figure out what output we need
#Yes. Use the “--plus T” option to specify use of BLAST+ instead of the “old-school” BLAST.
#megablast will speed it up but reduce accuracy. I don't think we care for this preliminary search, so go ahead and leave it
>class
for file in *1.fastq; do
echo "metaxa2 -1 $file -2 ${file/.1.fastq}.2.fastq -o ${file/.1.fastq} --plus T --split-pairs T --cpu 16 --graphical F --table F" >>class; done
#change this based on the number of jobs you need! Do
wc -l class
#Take the result of this and divide it by 24
#if is <1 set N = 1
#if it is >1 set N = 2 below
ls6_launcher_creator.py -j class -n class -a IBN21018 -e YOUREMAIL@utexas.edu -t 03:00:00 -N 2 -w 24 -q normal
sbatch class.slurm

#This is going to create a ton of files in your directory:
#-rw-r--r-- 1 cbscott G-801666 1861559 Mar 27 12:47 Buck-A20-90-16S.extraction.results
#-rw-r--r-- 1 cbscott G-801666 5204074 Mar 27 12:47 Buck-A20-90-16S.extraction.fasta
#-rw-r--r-- 1 cbscott G-801666       0 Mar 27 12:52 Buck-A20-90-16S.eukaryota.fasta
#-rw-r--r-- 1 cbscott G-801666    3692 Mar 27 12:52 Buck-A20-90-16S.uncertain.fasta
#-rw-r--r-- 1 cbscott G-801666 1248408 Mar 27 12:52 Buck-A20-90-16S.taxonomy.txt
#-rw-r--r-- 1 cbscott G-801666  453449 Mar 27 12:52 Buck-A20-90-16S.summary.txt
#-rw-r--r-- 1 cbscott G-801666     943 Mar 27 12:52 Buck-A20-90-16S.mitochondria.fasta
#-rw-r--r-- 1 cbscott G-801666   39908 Mar 27 12:52 Buck-A20-90-16S.chloroplast.fasta
#-rw-r--r-- 1 cbscott G-801666 4593579 Mar 27 12:52 Buck-A20-90-16S.bacteria.fasta
#-rw-r--r-- 1 cbscott G-801666     932 Mar 27 12:52 Buck-A20-90-16S.archaea.fasta

#For each of your samples, grab the bacterial and archaeal reads:
#(this might take a second, but should be ok to do on login node - doesn't need to be a job)
for file in *.summary.txt; do
>${file/.summary.txt}.roi &&\
sed '/Sequences of archaeal origin/,/---/!d;//d' $file >> ${file/.summary.txt}.roi &&\
sed '/Sequences of bacterial origin/,/---/!d;//d' $file >> ${file/.summary.txt}.roi; done

>subset
for file in *.roi; do
echo "seqtk subseq ${file/.roi}.1.fastq $file > ${file/.roi}.filtered.1.fastq &&\
seqtk subseq ${file/.roi}.2.fastq $file > ${file/.roi}.filtered.2.fastq " >> subset; done
ls6_launcher_creator.py -j subset -n subset -a IBN21018 -e YOUREMAIL@utexas.edu -t 01:00:00 -N 1 -w 24 -q vm-small

#determine percent read loss after filtering
for file in *.filtered.1.fastq; do
  echo ${file/.filtered.1.fastq} >> names &&\
  orig=`cat ${file/.filtered.1.fastq}.1.fastq | wc -l` &&\
  origreads=`expr $orig / 4` &&\
  new=`cat $file | wc -l` &&\
  newreads=`expr $new / 4` &&\
  echo $origreads >> startreads &&\
  echo $newreads >> endreads; done

paste names startreads endreads > filtering_read_retention

#does this look reasonable? The values in the third column should be close to
#the values in the second column.
head filtering_read_retention


#Rezip all of your files to save space.
>zipper
for file in *.fastq; do
echo "gzip $file">> zipper; done
ls6_launcher_creator.py -j zipper -n zipper -a IBN21018 -e YOUREMAIL@utexas.edu -t 00:30:00 -N 1 -w 24 -q vm-small
sbatch zipper.slurm

##### Now, let's move these files back to your local computer.
#To do this, first remove the files you downloaded last time. Let me know if you have questions on this.
#Go to your working directory (on your local computer), then to "Data/YourProjectReads" and delete all
#of the .gz files as well as any folders that were created by DADA2 (e.g., filtered)
#The only things you should keep are any metadata files or scripts that might be here.

#In your TACC terminal, do
pwd

#copy this path.

#then in ANOTHER terminal window (get one using command+N)
scp yourusername@ls6.tacc.utexas.edu:PATHYOUCOPIED/*.filtered*gz ~/The/Path/To/Your/Project/Folder/Data

#Remember, you can get the path to your project folder from finder. It might be helpful to copy it here
#in the script for future reference.


#Then, run DADA2 exactly as we discussed last week!
