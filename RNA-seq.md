# RNA-seq

## Check FASTQ quality
```bash  
module load fastqc/0.11.8

for i in *.fastq*;
do  
fastqc -t 8 $i  
done  
```

## Trim adapters if necessary
```bash  
module load trimmomatic/0.36  

IN=/path/to/fastq  
ADAPTERS=/path/to/trimmomatic/0.36/adapters  
R1=$IN/*R1.fastq.gz  
R2=$IN/*R2.fastq.gz  

trimmomatic PE -phred33 $R1 $R2 ${R1}_pairedout ${R1}_unpairedout ${R2}_pairedout ${R2}_unpairedout ILLUMINACLIP:$ADAPTERS/TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 AVGQUAL:20 MINLEN:36  
```

## Generate STAR genome directory
```bash
module load star/2.5.3a  

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir STAR --genomeFastaFiles GRCh37.p13.genome.fa --sjdbGTFfile gencode.v19.annotation.gtf --sjdbOverhang 149
```

## Alignment
```bash
genomeDir=/path/to/STAR/  
REF=/path/to/GRCh37.p13.genome.fa  
GTF=/path/to/gencode.v19.annotation.gtf  

FASTQ1='fastq/*L005_001.R1.fastq.gz'  
FASTQ2='fastq/*L005_001.R2.fastq.gz'  

STAR --runThreadN 8 --genomeDir $genomeDir --sjdbGTFfile $GTF --sjdbOverhang 149 --bamRemoveDuplicatesType UniqueIdentical --readFilesIn $FASTQ1 $FASTQ2 --twopassMode Basic --outSAMtype BAM SortedByCoordinate Unsorted --quantMode TranscriptomeSAM GeneCounts --readFilesCommand zcat
```

## Merge BAMS if necessary
```bash
BAM1='Run3/Aligned.sortedByCoord.out.bam'  
BAM2='Run4/Aligned.sortedByCoord.out.bam'  
BAM3='Run5/Aligned.sortedByCoord.out.bam'  
samtools merge RH1_merged.bam $BAM1 $BAM2 $BAM3
```

## Remove ribosomal and mitochondrial RNA
```bash
module load bedtools/2.29.2  

bedtools intersect -a Aligned.sortedByCoord.out.bam -b hg19_rRNA.bed -v > Sample_X_alignedSorted_rRNArm.bam  
bedtools intersect -a Sample_L1_alignedSorted_rRNArm.bam -b hg19_MT.bed -v > Sample_X_alignedSorted_rRNArm_MTrm.bam  
```

## Get Count Matrix
```bash
module load subread/1.5.2  

GTF=/path/to/gencode.v19.annotation.gtf  
IN=/path/to/BAMS  
OUT=/path/to/counts  

featureCounts -T 8 -p -t gene -s 1 -O -F GTF -a $GTF -o $OUT/Counts_gene_stranded.txt $IN/*.bam
```

# Differential Gene Expression
```R
library(DESeq2)  
library(dplyr)  
library(ggplot2)    

#Read in counts and define meta table, refine sample list if needed**  
counttab <- read.table("/path/to/counts/Counts_gene_stranded.txt", header=T, sep="\t", stringsAsFactors=F, check.names=F, row.names=1)  
counttab <- na.omit(counttab)  

#metatab <- read.table("~/path/to/counts/metatable.txt", header=T, sep="\t", stringsAsFactors=F, check.names=F)  
metatab <- data.frame(matrix(ncol = 3, nrow = 24))   
row.names(metatab) <- names(counttab)  
names(metatab) <- c("celltype", "group", "batch")  
metatab$celltype <- as.factor(c(rep("Neuron", 12), rep("NPC", 12)))  
metatab$group<- as.factor(rep(c("BFP", "D145N", "RH1"), 1, each = 4))   
metatab$batch <- as.factor(c(rep(1:2, 3, each = 2)))  

refinedSamples <- metatab$celltype == "Neuron" & metatab$group != "D145N"  
counttab <- counttab[,refinedSamples]  
metatab <- metatab[refinedSamples,]  

all(rownames(metatab) == colnames(counttab))

#Create DESeq2 object

dds <- DESeqDataSetFromMatrix(countData = counttab,
                                 colData = metatab,
                                 design = ~ group + batch)  
dds <- dds[rowSums(counts(dds) >= 10) >= 2,]  
```

## Stabilize variance for multidimensional scaling
```r
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)    
library(genefilter)  
 
dds <- estimateSizeFactors(dds)  

options(repr.plot.height=3, repr.plot.width=5)  
head(assay(dds), 3)  
par(mfrow=c(1,2))  
meanSdPlot(assay(dds), ranks=FALSE)  

options(repr.plot.height=3, repr.plot.width=5)  
rld <- rlog(dds, blind = FALSE)  
head(assay(rld), 3)  
meanSdPlot(assay(rld), ranks=FALSE)  

vsd <- vst(dds, blind = FALSE)  
head(assay(vsd), 3)  
meanSdPlot(assay(vsd), ranks=FALSE)  

df <- bind_rows(  
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%  
         mutate(transformation = "log2(x + 1)"),  
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),  
  as_data_frame(assay(rlog)[, 1:2]) %>% mutate(transformation = "rlog"))  

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")  
df$transformation <- factor(df$transformation, levels=lvls)  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +  
  coord_fixed() + facet_grid( . ~ transformation)  

par(mfrow=c(1,3))  
boxplot(log2(assay(dds)+1), las=2, main="log2(x+1)")  
boxplot(assay(rld), las=2, main="rld")  
boxplot(assay(vsd), las=2, main="vsd")  

sampleDists <- dist(t(assay(rld)))  
options(repr.plot.height=3, repr.plot.width=5)  
sampleDistMatrix <- as.matrix( sampleDists )  
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep = " - " )  
colnames(sampleDistMatrix) <- NULL  
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)  
pheatmap(sampleDistMatrix,  
         clustering_distance_rows = sampleDists,  
         clustering_distance_cols = sampleDists,  
         col = colors)  

poisd <- PoissonDistance(t(counts(dds)))  
options(repr.plot.height=3, repr.plot.width=5)  
samplePoisDistMatrix <- as.matrix( poisd$dd )  
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep=" - " )  
colnames(samplePoisDistMatrix) <- NULL  
pheatmap(samplePoisDistMatrix,  
         clustering_distance_rows = poisd$dd,  
         clustering_distance_cols = poisd$dd,  
         col = colors)  

plotPCA(rld, intgroup = c("group", "celltype"))  

mds <- cbind(as.data.frame(colData(rld)), cmdscale(sampleDistMatrix))  
ggplot(mds, aes(x = `1`, y = `2`, color = group, shape = celltype)) +  
  geom_point(size = 3) + coord_fixed()  
```

## DESEQ2
```r
dds <- DESeq(dds)  
res <- results(dds, contrast=c("group","BFP","RH1"))    
sig_res <- results(dds, alpha = 0.05)  
table(sig_res$padj < 0.05)  
summary(sig_res)  
```

## Annotation
```r  
library(AnnotationDbi)  
library(org.Hs.eg.db)  

ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

write.csv(res, "/path/to/deseq.csv", quote=FALSE)
```

# Visualization of DEG results
```r
library(ggbeeswarm)  

#Top gene  
topGene <- rownames(res)[which.min(res$padj)]
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("group", "celltype"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = group, y = count, color = group)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

#Specific gene  
rnaseh1 <- plotCounts(dds, gene='ENSG00000171865.5', intgroup=c("group"), returnData=TRUE)  
ggplot(rnaseh1, aes(x=group, y=count, col=group)) +  
  geom_point(position=position_jitter(width=.1,height=0)) +  
  scale_y_log10()  
```

### Heatmap
```r
library(circlize)  
library(ComplexHeatmap)  

colors <- colorRamp2(seq(min(mat), max(mat), length = 3), c("blue", "#EEEEEE", "red"))

resSort <- res[order(res$padj),]  
topgenes <- head(rownames(resSort),100)  
mat <- assay(rld)[topgenes,]  
mat <- mat - rowMeans(mat)  
Heatmap(mat, name = "Normalized Counts", col = colors, column_title = "BFP vs. RNase H1")  
```
### Volcano Plot
```r  
library(EnhancedVolcano)  

EnhancedVolcano(res,
                lab = res$symbol,  
                x = 'log2FoldChange',  
                y = 'pvalue',  
                title = 'RNase H1 vs. BFP',  
                pCutoff = 10e-6,  
                FCcutoff = 1,  
                cutoffLineType = 'twodash',  
                cutoffLineWidth = 0.5,   
                pointSize = 2.0,  
                labSize = 5.0,  
                colAlpha = 1,  
                legendLabels=c('NS','Log2FC','p-value', 'p-value & Log2FC'),  
                col=c('cornsilk4', 'darkblue', 'darkorchid', 'deeppink2'),  
                legendPosition = 'right',  
                legendLabSize = 16,  
                legendIconSize = 5.0,  
                drawConnectors = TRUE,  
                widthConnectors = 0.1)  
```

# Gene Set Enrichment Analysis

```r
library(clusterProfiler)  
library(enrichplot)  
library(organism, character.only = TRUE)  
organism = org.Hs.eg.db  
degenes = read.csv("/path/to/deseq.csv ", header=TRUE)  

original_gene_list <- degenes$log2FoldChange  
names(original_gene_list) <- degenes$symbol  
gene_list<-na.omit(original_gene_list)  
gene_list = sort(gene_list, decreasing = TRUE)  
gse <- gseGO(geneList=gene_list,  
             ont ="ALL",  
             keyType = "SYMBOL",  
             nPerm = 10000,  
             minGSSize = 3,  
             maxGSSize = 800,  
             pvalueCutoff = 0.05,   
             verbose = TRUE,   
             OrgDb = organism,   
             pAdjustMethod = "none")  

#Dotplot  
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)  

#Enrichment Map   
emapplot(gse, showCategory = 10)   

#Category Map: categorySize can be either 'pvalue' or 'geneNum'  
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)  
```
# Integrate RNA-seq & ChIP-seq data  

```bash
module load python/2.7.16 py_packages anaconda2

DEEPTOOLS=/path/to/programs/deepTools/bin  
IN=/path/to/bigwigs  
OUT=/path/to/out  
BLACKLIST=/path/to/hg19-blacklist.bed  
GTF=/path/to/hg19_gtf_gene.bed  

$DEEPTOOLS/multiBigwigSummary BED-file \  
 --BED $GTF \  
 -b $IN/*.bw \  
 --blackListFileName $BLACKLIST \  
 -o $OUT/ChIPvRNA_GTF_multibigwigsummary_blrm_rmdup_min30.npz  

 $DEEPTOOLS/plotPCA \  
 -in $OUT/ChIPvRNA_GTF_multibigwigsummary_blrm_rmdup_min30.npz \  
 -o $OUT/ChIPvRNA_GTF_multibigwigsummary_blrm_rmdup_min30_PCA.pdf \  
 -T "correlation"  

 $DEEPTOOLS/plotCorrelation \  
 -in $OUT/ChIPvRNA_GTF_multibigwigsummary_blrm_rmdup_min30.npz \  
  --skipZeros \  
 --removeOutliers \  
 --corMethod spearman \  
 --whatToPlot heatmap \  
 --colorMap Spectral_r \  
 --skipZeros \  
 --plotTitle "Spearman Correlation" \  
 -o $OUT/ChIPvRNA_GTF_multibigwigsummary_blrm_rmdup_min30_SPEARMAN_HEATMAP.pdf \  
 --outFileCorMatrix $OUT/ChIPvRNA_GTF_multibigwigsummary_blrm_rmdup_min30_SPEARMAN_HEATMAP_CORRELATION.tab  

 $DEEPTOOLS/plotCorrelation \  
 -in $OUT/ChIPvRNA_GTF_multibigwigsummary_blrm_rmdup_min30.npz \  
  --skipZeros \  
 --removeOutliers \  
 --corMethod pearson \  
 --whatToPlot heatmap \  
 --colorMap Spectral_r \  
 --skipZeros \  
 --plotTitle "Pearson Correlation" \  
 -o $OUT/ChIPvRNA_GTF_multibigwigsummary_blrm_rmdup_min30_PEARSON_HEATMAP.pdf \  
 --outFileCorMatrix $OUT/ChIPvRNA_GTF_multibigwigsummary_blrm_rmdup_min30_PEARSON_HEATMAP_CORRELATION.tab  

 $DEEPTOOLS/plotCorrelation \  
 -in $OUT/ChIPvRNA_GTF_multibigwigsummary_blrm_rmdup_min30.npz \  
  --skipZeros \  
 --removeOutliers \  
 --corMethod spearman \  
 --whatToPlot scatterplot \  
 --colorMap Spectral_r \  
 --skipZeros \  
 --plotTitle "Spearman Correlation" \  
 -o $OUT/ChIPvRNA_GTF_multibigwigsummary_blrm_rmdup_min30_SPEARMAN_SCATTER.pdf \  
 --outFileCorMatrix $OUT/ChIPvRNA_GTF_multibigwigsummary_blrm_rmdup_min30_SPEARMAN_SCATTER_CORRELATION.tab  

 $DEEPTOOLS/plotCorrelation \  
 -in $OUT/ChIPvRNA_GTF_multibigwigsummary_blrm_rmdup_min30.npz \  
  --skipZeros \  
 --removeOutliers \  
 --corMethod pearson \  
 --whatToPlot scatterplot \  
 --colorMap Spectral_r \  
 --skipZeros \  
 --plotTitle "Pearson Correlation" \  
 -o $OUT/ChIPvRNA_GTF_multibigwigsummary_blrm_rmdup_min30_PEARSON_SCATTER.pdf \  
 --outFileCorMatrix $OUT/ChIPvRNA_GTF_multibigwigsummary_blrm_rmdup_min30_PEARSON_SCATTER_CORRELATION.tab  
```
