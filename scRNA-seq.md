# scRNA-seq

## Load in required packages in R
```r
library(Seurat)
library(dplyr)
library(patchwork)
library(cellranger)
library(ggplot2)
library(Matrix)
library(reshape2)
library(WebGestaltR)
library(stringr)
library(clustree)
```

## Read in 10X data
```r
rloopumi <- Read10X(data.dir = "./rloop/")
colnames(rloopumi) <- do.call("c",lapply(strsplit(colnames(rloopumi),"-",fixed=TRUE),function(x) x[1]))
mtx <- readMM("./HTO/matrix.mtx.gz")
tag <- read.table("./HTO/features.tsv.gz", header=F)
barcode <- read.table("./HTO/barcodes.tsv.gz", header=F)
rownames(mtx) <- lab[match(tag[,1],lab[,1]),7]
colnames(mtx) <- barcode[,1]
rloophto <- mtx
joint.bcs <- intersect(colnames(rloopumi), colnames(rloophto))

#Subset RNA and HTO counts by joint cell barcodes
rloopumi <- rloopumi[, joint.bcs]
rloophto <- as.matrix(rloophto[, joint.bcs])

#Confirm that the HTO have the correct names
rownames(rloophto)
```

## Create Seurat object from loaded 10X Data
```r
rloop <- CreateSeuratObject(counts = rloopumi)
```

## Extract the singlets, filter data, & single cell transform
```r
rloop[["percent.mt"]] <- PercentageFeatureSet(rloop, pattern = "^MT-")
rloop.singlet <- subset(rloop, idents = "Singlet")
rloop <- subset(rloop.singlet, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000  & percent.mt < 25)
list <- SplitObject(rloop, split.by = "hash.ID")
for (i in 1:length(list)) {
  list[[i]] <- SCTransform(list[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
}
```

## Integration
```r
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
options(future.globals.maxSize = 73400032000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features,
                                    verbose = FALSE)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                           anchor.features = features, verbose = FALSE)
rloop <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                                     verbose = FALSE)

DefaultAssay(rloop) <- "integrated"
rloop <- RunPCA(rloop, verbose = FALSE)
rloop <- RunUMAP(rloop, dims = 1:30)
rloop <- FindNeighbors(rloop, dims = 1:30)
```

## Test best resolution with Clustree
```r
res.range<-seq(from=0, to=1, by=0.1)
rloop <- FindClusters(rloop, resolution = res.range, verbose=F)
clustree(rloop, prefix = "integrated_snn_res.")
rloop <- FindClusters(rloop, resolution = 0.35)
DefaultAssay(rloop) <- "RNA"
rloop <- NormalizeData(rloop, verbose = FALSE)
```

## Find cluster marker genes and DEGs
```r
DimPlot(rloop, label =T)
rloop.markers<-FindAllMarkers(rloop, only.pos = T)
write.csv(rloop.markers, "Cluster_markers_logfc0.csv")
```

## Name cluster markers using CSEA
```r
Idents(rloop)<-rloop$integrated_snn_res.0.35
rloop <- RenameIdents(rloop, `0` = "Glia", `1` = "Glia",
                      `2` = "Glia", `3` = "Neurons", `4` = "Unident",
                      `5` = "Glia", `6` = "Glia",
                      `7` = "Unident", `8` = "NPC",
                      `9` = "NPC", `10` = "Neurons",
                      `11` = "Unident", `12` = "Unident",
                      `13` = "Neurons")
rloop$celltype <- Idents(rloop)
```

## Break into groups
```r
rloop$treatment<-ifelse(rloop$hash.ID=="553-RH-rep1", "RH", "BFP")
rloop$treatment2<-ifelse(rloop$hash.ID=="2607-RH-rep2", "RH", rloop$treatment)
rloop$treatment3<-ifelse(rloop$hash.ID=="553-RH-rep2", "RH", rloop$treatment2)
rloop$treatment4<-ifelse(rloop$hash.ID=="2607-RH-rep1", "RH", rloop$treatment3)
Idents(rloop)<-rloop$treatment4
```

## Assign RNASEH1 level
```r
RH <- as.data.frame(rloop@assays$SCT@data)
RH2 <- RH["RNASEH1",]
RH2 <- as.data.frame(t(RH2))
rloop$RH_level <- RH2
```

## Subgroup based on RNASEH1 expression level
```r
RH_hi <- WhichCells(rloop, expression = RNASEH1 > 1.25)
rloop$RH_hi <- ifelse(rownames(rloop@meta.data) %in% RH_hi, "hi", "lo")
```

## DEGs based on RNASEH1 and split into appropriate groups
```r
rloop$level <- paste(rloop$RH_hi, rloop$treatment4, sep = "_")
Idents(rloop)<-rloop$level

#Overall
RH_hiVlo <- FindMarkers(rloop, ident.1 = "RH_hi", ident.2 = "BFP_lo", assay = "RNA", slot = "data", logfc.threshold = 0)
write.csv(RH_hiVlo, "DEG_RHhivBFPlo_logfc0.csv")

#Glia
RHvBFP_gli<-FindMarkers(rloop, ident.1 = "Glia_RH_hi", ident.2 = "Glia_BFP_lo", assay = "RNA", slot = "data", logfc.threshold = 0)
write.csv(RHvBFP_gli, "DEG_Glia_RHhivBFPlo_logfc.csv")

#Neurons
RHvBFP_neu<-FindMarkers(rloop, ident.1 = "Neurons_RH_hi", ident.2 = "Neurons_BFP_lo", assay = "RNA", slot = "data", logfc.threshold = 0)
write.csv(RHvBFP_neu, "DEG_Neuron_RHhivBFPlo_logfc.csv")

#NPCs
RHvBFP_npc<-FindMarkers(rloop, ident.1 = "NPC_RH_hi", ident.2 = "NPC_BFP_lo", assay = "RNA", slot = "data", logfc.threshold = 0,)
write.csv(RHvBFP_npc, "DEG_NPC_RHhivBFPlo_logfc.csv")
```

## Cell Phase
```r
s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
rloop <- CellCycleScoring(rloop, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ggplot(rloop[[]]) +
  aes(x=celltype, fill=factor(Phase)) +
  geom_bar(position = "fill")
Idents(rloop)<-rloop$Phase
DimPlot(rloop)

cycle<-as.data.frame(table(rloop@meta.data$group, rloop@meta.data$Phase))
cycle<-as.data.frame(table(rloop@meta.data$type, rloop@meta.data$Phase))
write.csv(cycle, "group_phase.csv")
```

# Data Visualization

## Volcano Plot
```r
library(EnhancedVolcano)
degenes <- read.csv("DEG_RHhivBFPlo_logfc0.csv", header=TRUE)  

keyvals <- ifelse(
  RHvBFP$p_val_adj > 0.5, '#2E008B',
  ifelse(RHvBFP_npc$avg_log2FC > 0, '#CC339C',
         '#0DDE31'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == '#CC339C'] <- 'high'
names(keyvals)[keyvals == '#0DDE31'] <- 'mid'
names(keyvals)[keyvals == '#2E008B'] <- 'low'

p1 <- EnhancedVolcano(RHvBFP,
                      lab = rownames(RHvBFP),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      selectLab = c('RNASEH1'),
                      #selectLab = rownames(RHvBFP)[which(names(keyvals) %in% c('high', 'low'))],
                      colAlpha = 1,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      cutoffLineType = 'blank',
                      cutoffLineCol = 'black',
                      #colCustom = keyvals
                      )

p1 +
  ggplot2::coord_cartesian(xlim=c(-4, 4))
```

## Gene Set Enrichment Analysis and Visualization
```r
library(clusterProfiler)  
library(enrichplot)  
library(organism, character.only = TRUE)  
organism = org.Hs.eg.db
degenes = read.csv("DEG_RHhivBFPlo_logfc0.csv", header=TRUE)  

original_gene_list <- degenes$avg_log2FC  
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
```

### Dotplot  
```r
dotplot(gse, showCategory=5, split=".sign") + facet_grid(.~.sign)  
```

### Enrichment Map   
```r
emapplot(gse, showCategory = 10)   
```

### Category Map
```r
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)  
```
