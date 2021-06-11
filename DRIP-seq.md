# DRIP-seq

## Check FASTQ quality
``` Bash
module load fastqc/0.11.8

for i in *.fastq*;
do  
fastqc -t 8 $i  
done  
```

## Trim adapters if necessary
``` bash
module load trimmomatic/0.36  

IN=/path/to/fastq  
ADAPTERS=/path/to/trimmomatic/0.36/adapters  
R1=$IN/*R1.fastq.gz  
R2=$IN/*R2.fastq.gz  

trimmomatic PE -phred33 $R1 $R2 ${R1}_pairedout ${R1}_unpairedout ${R2}_pairedout ${R2}_unpairedout ILLUMINACLIP:$ADAPTERS/TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 AVGQUAL:20 MINLEN:36\
```

## Alignment
``` Bash
module load bowtie2 samtools igvtools bedtools

directories='*'
for folder in $directories
do

cd $folder
FQFILES='Data/*R1.fastq*'
BOWTIE2INDEX='/path/to/bowtie2/hg19'
mkdir Alignments
mkdir FastQC

for FQFILE in $FQFILES;
do
   echo "Processing file: $FQFILE"
   FQFILE2=${FQFILE/R1/R2}
   SAMFILE=${FQFILE/Data/Alignments}
   SAMFILE=${SAMFILE/fastq/sam} 
   echo $"1:START bowtie alignment with output file $SAMFILE"
   bowtie2 -x $BOWTIE2INDEX -1 $FQFILE -2 $FQFILE2 > $SAMFILE
   echo $'1:FINISHED bowtie alignment\n'

   BAMFILE=${SAMFILE/sam/bam}
   echo '2:START SAM to BAM conversion'
   samtools view -b -h -S -F 4 $SAMFILE > $BAMFILE
   echo $'2:FINISHED SAM to BAM conversion\n'

   echo '3:START fastqc'
   fastqc -t 8 $BAMFILE > fastqc
   echo $'3:FINISHED fastqc\n'

   BAMSORTEDPREFIX=${BAMFILE/bam/sorted}
   echo '4:START BAM sorting'
   samtools sort -m 4G $BAMFILE $BAMSORTEDPREFIX
   echo $'4:FINISHED BAM sorting\n'

   BAMSORTEDFILE=${BAMFILE/bam/sorted.bam} 

   echo '5:START BAM indexing'
   samtools index $BAMSORTEDFILE
   echo $'5:FINISHED BAM indexing\n'

   RMDUPFILE=${BAMSORTEDFILE/bam/rmdup.bam}
   echo '6:START removal of PCR duplicates'
   samtools rmdup -s $BAMSORTEDFILE $RMDUPFILE
   echo $'6:FINISHED removal of PCR duplicates\n'

   echo '7:START remove intermediate files'
   echo "   rm $SAMFILE $BAMFILE $BAMSORTEDFILE"
   echo $'7:FINISHED remove intermediate files\n'

   FLAGSTATFILE=${RMDUPFILE/bam/flagstat}
   echo '8:START flagstat report'
   samtools flagstat $RMDUPFILE > $FLAGSTATFILE
   echo $'8:FINISHED flagstat report\n'

   TDFFILE=${RMDUPFILE/bam/tdf}
   echo '9:START TDF conversion'
   igvtools count $RMDUPFILE $TDFFILE hg19
   echo $'9:FINISHED TDF conversion\n'
done

echo 'ALL PROCESSING COMPLETED'

done
```

## Remove ENCODE Blacklist and MAPQ Score < 30
``` Bash
module load samtools bedtools

OUT=/path/to/out
BLACKLIST=/path/to/hg19-blacklist.bed

for i in *.bam;
do
samtools view -q 30 -b $i | bedtools intersect -a stdin -b $BLACKLIST -v > $OUT/${i}_BLrm_MAPQ30.bam
done
```

## Create Tag Directory
``` Bash
module load homer

OUT=/path/to/out

for j in *.bam;
do
makeTagDirectory $OUT/${j}_tag
done
```

## Call Peaks
``` Bash
module load homer

OUT=/path/to/out
control=Input.Sample.bam_tag

for j in *_tag;
do
findPeaks $j -style histone -o $OUT/${j}.fdr0.05.peaks.txt -i $control -fdr 0.05
done
```

## Diffbind (R)
``` r
library("DiffBind")

#Make metadata file
samples <- read.csv("metadata.csv", header=TRUE, sep=",")
DRIP <- dba(sampleSheet=samples)
DRIP <- dba.count(DRIP)
DRIP <- dba.contrast(DRIP, group1 = DRIP$masks$GM, group2=DRIP$masks$CP)
DRIP_de <- dba.analyze(DRIP,method=c(DBA_DESEQ2), bFullLibrarySize=F) #Use bFullLibrarySize=T for RNase H1-expressing neurons and controlds
DRIP_de$config$th <- 0.05
dba.plotVolcano(DRIP_de, contrast=1, method=DBA_DESEQ2)

DRIP.DB <- dba.report(DRIP_de)
out <- as.data.frame(DRIP.DB)
write.table(out, file="DRIP.txt", sep="\t", quote=F, row.names=F)
```

## Annotate
``` r
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peak <-readPeakFile("DRIP.txt")
peakAnno <- annotatePeak(peak, TxDb=txdb, annoDb="org.Hs.eg.db", 
                         tssRegion=c(-3000, 3000), verbose=FALSE)

plotAnnoBar(peakAnno)
annotation <- as.data.frame(peakAnno)
write.table(annotation, file="DRIP.annotated.txt", sep="\t", quote=F, row.names=F, col.names=F)
```
