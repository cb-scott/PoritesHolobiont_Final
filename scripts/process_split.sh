conda activate map_software
bs list projects
bs download project -i 384004281 -o amp_seq

git clone https://github.com/Y-Lammers/Split_on_Primer.git

#needs a primer file
#let's make one:


cutadapt \
  -g ^AAGGCC \
  -G ^TTGGAA \
  --untrimmed-output=untrimmed.1.fastq.gz \
  --untrimmed-paired-output=untrimmed.2.fastq.gz
  -o trimmed.2.fastq.gz \
  -p trimmed.2.fastq.gz \
  input.1.fastq.gz \
  input.2.fastq.gz

cutadapt -a ^FWDPRIMER...RCREVPRIMER -A ^REVPRIMER...RCFWDPRIMER --discard-untrimmed -o out.1.fastq.gz -p out.2.fastq.gz in.1.fastq.gz in.2.fastq.gz


>sep_files
for file in *AND*R1*; do
name=${file##*AND-}
#ITS2
echo "cutadapt -a ^GAATTGCAGAACTCCGTGAACC...GCATGAAGTCAYACAAGWGAACCCG -A ^CGGGTTCWCTTGTYTGACTTCATGC...GGTTCACGGAGTTCTGCAATTC --discard-untrimmed -o ${name/_*_L001_R1_001.fastq.gz}.1.fastq.gz -p ${name/_*_L001_R1_001.fastq.gz}.2.fastq.gz $file ${file/R1_001.fastq.gz}R2_001.fastq.gz" >> sep_files &&
#16S
echo "cutadapt -a ^CCTACGGGNGGCCTACGGGNGGCWGCAG...GGATTAGATACCCVHGTAGTC -A ^GACTACHVGGGTATCTAATCC...CTGCWGCCNCCCGTAGGCCNCCCGTAGG --discard-untrimmed -o ${file%-AND*}.1.fastq.gz -p ${file%-AND*}.2.fastq.gz $file ${file/R1_001.fastq.gz}R2_001.fastq.gz" >> sep_files; done
python3 $HOME/bin/ls6_launcher_creator.py -j sep_files -n sep_files -a IBN21018 -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 14 -q development

mkdir sourcefiles
mv *AND* sourcefiles
#all files should have the same naming convention
>trim16s
for file in *16S_S*R1*; do
  echo "cutadapt -a ^CCTACGGGNGGCCTACGGGNGGCWGCAG...GGATTAGATACCCVHGTAGTC -A ^GACTACHVGGGTATCTAATCC...CTGCWGCCNCCCGTAGGCCNCCCGTAGG --discard-untrimmed -o ${file/_S*_L001_R1_001.fastq.gz}.1.fastq.gz -p ${file/_S*_L001_R1_001.fastq.gz}.2.fastq.gz $file ${file/R1_001.fastq.gz}R2_001.fastq.gz" >> trim16s; done
python3 $HOME/bin/ls6_launcher_creator.py -j trim16s -n trim16s -a IBN21018 -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 14 -q development

>trimITS2
for file in *ITS2_S*R1*; do
  echo "cutadapt -a ^CCTACGGGNGGCCTACGGGNGGCWGCAG...GGATTAGATACCCVHGTAGTC -A ^GACTACHVGGGTATCTAATCC...CTGCWGCCNCCCGTAGGCCNCCCGTAGG --discard-untrimmed -o ${file/_S*_L001_R1_001.fastq.gz}.1.fastq.gz -p ${file/_S*_L001_R1_001.fastq.gz}.2.fastq.gz $file ${file/R1_001.fastq.gz}R2_001.fastq.gz" >> trimITS2; done
python3 $HOME/bin/ls6_launcher_creator.py -j trimITS2 -n trimITS2 -a IBN21018 -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 14 -q development

#RUN CUTADAPT AGAIN ON SAMPLES TO REMOVE REMAINING PRIMER SEQS
