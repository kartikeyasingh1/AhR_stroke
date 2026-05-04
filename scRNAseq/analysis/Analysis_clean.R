---
title: "AHR_scRNAseq_analysis"
authors: "Alba Simats & Corinne Benakis"
date: "25/04/2026"
---

library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
library(openxlsx)
library(RColorBrewer)
library(patchwork)
library(readxl)
library(RColorBrewer)
library(ggplot2)

setwd("...")
getwd()

matrix_dir = "/Volumes/.../filtered_feature_bc_matrix/"
list.files(matrix_dir)

barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- Matrix::readMM(file = matrix.path)

barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

GEX <- mat[0:32285,] #check the feature.name file and adjust the rows here
GEX_feature.names <- feature.names[0:32285,]

HTO <- mat[32286:32289,]
HTO_feature.names <- feature.names[32286:32289,]
HTO_feature.names

colnames(GEX) = barcode.names$V1
colnames(HTO) = barcode.names$V1
rownames(GEX) = GEX_feature.names$V2
rownames(HTO) = HTO_feature.names$V1

rownames(GEX) <- make.unique(rownames(GEX))
rownames(HTO) <- paste0(rownames(HTO), "-HTO")


# Create seurat object 
s1 <- CreateSeuratObject (counts = GEX, project = "Ahr")
s1[["HTO"]] <- CreateAssayObject(counts = HTO)
s1 <- NormalizeData(s1, assay = "HTO", normalization.method = "CLR")
s1 <- HTODemux(s1, assay = "HTO", positive.quantile = 0.99)
table(s1$HTO_classification.global)

Idents(s1) <- "HTO_maxID"
RidgePlot(s1, assay = "HTO", features = rownames(s1[["HTO"]])[1:4], ncol = 2)

Idents(s1) <- "HTO_classification.global"
VlnPlot(s1, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
head(s1)

# Extract singlets
Ahr2 <- subset(s1, idents = c("Doublet", "Negative"), invert = TRUE)
table(Ahr2$orig.ident)

# Get batches based on cell names
sample <- sapply(colnames(GetAssayData(object = Ahr2, slot = "counts")),
                 FUN=function(x){substr(x,18,18)})
sample <- as.numeric(as.character(sample))
names(sample) <- colnames(GetAssayData(object = Ahr2, slot = "counts"))
Ahr2 <- AddMetaData(Ahr2, sample, "sample")

new.grouping <- c("organ")
Ahr2[[new.grouping]] <- new.grouping
colnames(Ahr2@meta.data)
Ahr2$organ[Ahr2$HTO_classification == "B0301-HTO" ] <- "Sp" #spleen
Ahr2$organ[Ahr2$HTO_classification == "B0302-HTO" ] <- "Br" #brain
Ahr2$organ[Ahr2$HTO_classification == "B0303-HTO" ] <- "Lp" #ileal lamina propria
Ahr2$organ[Ahr2$HTO_classification == "B0304-HTO" ] <- "Bl" #blood
table(Ahr2$organ)
saveRDS(Ahr2, "Ahr.rds")

Ahr2[["percent.mt"]] <- PercentageFeatureSet(Ahr2, pattern = "^mt-")
VlnPlot(object = Ahr2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(Ahr2, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Ahr2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Ahr2 <- subset(Ahr2, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(object = Ahr2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(Ahr2, "Ahr.rds")

library(sctransform)
Ahr2 <- SCTransform(Ahr2, vars.to.regress = "percent.mt", verbose = FALSE)
Ahr2 <- RunPCA(object = Ahr2, features = VariableFeatures(object = Ahr2), verbose = FALSE)
Ahr2 <- ProjectDim(object = Ahr2)
ElbowPlot(object = Ahr2)
DimHeatmap(object = Ahr2, dims = 8:16, cells = 500, balanced = TRUE)

Ahr2 <- FindNeighbors(object = Ahr2, dims = 1:10) 
Ahr2 <- FindClusters(object = Ahr2, resolution = 0.1) 
Ahr2 <- RunUMAP(object = Ahr2, dims = 1:10)
DimPlot(Ahr2, reduction = 'umap', label = TRUE, group.by = "seurat_clusters")
saveRDS(Ahr2, file = "Ahr.rds")

Ahr <-Ahr2
DimPlot(Ahr)
new.grouping <- c("experiment")
Ahr[[new.grouping]] <- new.grouping
colnames(Ahr@meta.data)
Ahr$experiment[Ahr$sample == "1" ] <- "WT"
Ahr$experiment[Ahr$sample == "2" ] <- "WT"
Ahr$experiment[Ahr$sample == "3" ] <- "WT"
Ahr$experiment[Ahr$sample == "4" ] <- "KO"
Ahr$experiment[Ahr$sample == "5" ] <- "KO"
Ahr$experiment[Ahr$sample == "6" ] <- "KO"
DimPlot(object = Ahr, group.by="experiment")


Ahr
Ahr@active.assay = "SCT"
Ahr <- FindVariableFeatures(object = Ahr, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
Ahr[["SCT"]] <- split(Ahr[["SCT"]], f = Ahr$experiment)
Ahr <- IntegrateLayers(object = Ahr, method = HarmonyIntegration, assay = "SCT", orig.reduction = "pca", 
                       new.reduction = 'harmony', verbose = FALSE)
Ahr[["RNA"]] <- JoinLayers(Ahr[["RNA"]])
Ahr <- FindNeighbors(object = Ahr, reduction = "harmony", dims = 1:14) #15
Ahr <- FindClusters(object = Ahr, resolution = 0.1) #0.2
Ahr <- RunUMAP(object = Ahr, dims = 1:14, reduction = "harmony")
saveRDS(Ahr, "Ahr_integrated.rds")

Ahr@active.assay = "RNA"
Ahr <- NormalizeData(object = Ahr, normalization.method = "LogNormalize", scale.factor = 1e4)
Ahr <- FindVariableFeatures(object = Ahr, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(object = Ahr))
Ahr <- ScaleData(object = Ahr, features = rownames(x = Ahr), vars.to.regress = c("nCount_RNA", "percent.mito"))
saveRDS(Ahr, "Ahr_integrated.rds")

## further analyses were performed on organs of interest separately ##
