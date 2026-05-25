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

# GSEA using the Gene Ontology (GO) database. 
Ahr_lp -> Ahr 

Layers(Ahr[["RNA"]])
options(spe = c("mouse"))
Ahr <- GeneSetAnalysisGO(Ahr, parent = "GO:0002376")
matr <- Ahr@misc$AUCell$GO$"GO:0002376"
matr <- RenameGO(matr)
head(matr, 4:3)

GeneSetAnalysisGO()
SeuratExtend::Heatmap(CalcStats(matr, f = Ahr_lp_noNKnoRibo$seurat_clusters, order = "p", n = 3), lab_fill = "zscore")

stats_cluster <- CalcStats(
  matr,
  f = Ahr$seurat_clusters,
  order = "p",
  n = 3
)

features_keep <- unique(rownames(stats_cluster))
cluster_cond <- interaction(
        Ahr$experiment,
        Ahr$seurat_clusters,
  sep = "_"
)

stats_final <- CalcStats(
  matr[features_keep, ],
  f = cluster_cond
)
SeuratExtend::Heatmap(
  stats_final,
  lab_fill = "zscore"
)

WaterfallPlot(matr, f = Ahr$experiment, ident.1 = "KO", ident.2 = "WT", top.n = 5)
p <- WaterfallPlot(matr, f = Ahr$experiment,ident.1 = "KO",ident.2 = "WT",style = "segment", color_theme = "D", top.n = 5, len.threshold = 2)

p + 
  ggtitle(" Enriched Pathways (immune_system_process) in Ahr lp ") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
    axis.text.y = element_text(size = 8)  # Smaller y-axis font
  )


# module score on a specific gene set #
AddModuleScore(
        Ahr,
  features=gene_DCTolerance_merge_GOandLiterature_indata,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "Cluster",
  seed = 1,
  search = FALSE,
  layer = "data"
)

Ahr <- AddModuleScore(
        Ahr,
  features = list(gene_DCTolerance_merge_GOandLiterature_indata),  
  name = "gene_DCTolerance_merge_GOandLiterature_indata_Module"
)

VlnPlot(Ahr, features = "gene_DCTolerance_merge_GOandLiterature_indata_Module1", 
        group.by = "seurat_clusters", cols=color_palette_lp_noNK, split.by = "experiment")


# DELTA AHR module (KO-WT)
df <- FetchData(
        Ahr,
  vars = c(
    "gene_DCTolerance_merge_GOandLiterature_indata_Module1", 
    "seurat_clusters",
    "experiment", 
    "sample"     
  )
)

df_mouse <- df %>%
  group_by(sample, experiment, seurat_clusters) %>%
  summarise(
    module_mean = mean(gene_DCTolerance_merge_GOandLiterature_indata_Module1),
    .groups = "drop"
  )

sum_stats <- df_mouse %>%
  group_by(seurat_clusters, experiment) %>%
  summarise(
    mean = mean(module_mean, na.rm = TRUE),
    sd   = sd(module_mean, na.rm = TRUE),
    n    = n(),
    se   = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  select(seurat_clusters, experiment, mean, se) %>%
  pivot_wider(names_from = experiment, values_from = c(mean, se))

df_effect2 <- sum_stats %>%
  mutate(
    delta_KO_WT = mean_KO - mean_WT,
    se_delta    = sqrt(se_KO^2 + se_WT^2)
  ) %>%
  select(seurat_clusters, delta_KO_WT, se_delta)

df_stats <- df_mouse %>%
  group_by(seurat_clusters) %>%
  wilcox_test(module_mean ~ experiment, exact = FALSE) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  select(seurat_clusters, p, p_adj)

df_plot <- df_effect2 %>%
  left_join(df_stats, by = "seurat_clusters") %>%
  mutate(
    seurat_clusters = factor(seurat_clusters),
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )
df_plot <- df_plot %>%
  arrange(delta_KO_WT) %>%
  mutate(seurat_clusters = factor(seurat_clusters, levels = seurat_clusters))

ggplot(df_plot, aes(x = seurat_clusters, y = delta_KO_WT)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = delta_KO_WT - se_delta,
                    ymax = delta_KO_WT + se_delta),
                width = 0.25) +
  geom_point(size = 3) +
  geom_text(aes(label = sig),
            vjust = -1.1, size = 4) +
  labs(
    x = "Cluster",
    y = expression(Delta*" AHR module (KO - WT)")
  ) +
  theme_classic()+ labs(
    subtitle = "Mean ± SE of KO−WT difference per cluster; n = 3 mice per genotype\nWilcoxon tests per cluster; BH-FDR across clusters"
  )

# Sample-level stats for module score per cluster (KO vs WT)
obj <- Ahr
score_col   <- "gene_DCTolerance_merge_GOandLiterature_indata_Module1" 
cluster_col <- "seurat_clusters"
cond_col    <- "experiment"   
sample_col  <- "sample"       

md <- obj@meta.data

df_sample <- md %>%
  dplyr::select(all_of(c(sample_col, cond_col, cluster_col, score_col))) %>%
  dplyr::filter(!is.na(.data[[score_col]])) %>%
  dplyr::group_by(
    sample    = .data[[sample_col]],
    condition = .data[[cond_col]],
    cluster   = .data[[cluster_col]]
  ) %>%
  dplyr::summarise(
    score   = mean(.data[[score_col]], na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

sample_counts <- df_sample %>%
  group_by(cluster, condition) %>%
  summarise(n_samples = n_distinct(sample), .groups = "drop") %>%
  pivot_wider(
    names_from = condition,
    values_from = n_samples,
    values_fill = 0
  ) %>%
  arrange(cluster)
print(sample_counts)

eligible_clusters <- sample_counts %>%
  filter(KO >= 3, WT >= 3) %>%
  pull(cluster)

message("Clusters retained for analysis: ",
        paste(eligible_clusters, collapse = ", "))

df_sample_filt <- df_sample %>%
  filter(cluster %in% eligible_clusters)

stats_tbl <- map_dfr(sort(unique(df_sample_filt$cluster)), function(cl) {
  df_cl <- df_sample_filt %>% filter(cluster == cl)
  tst <- wilcox_test(df_cl, score ~ condition)
  tibble(
    cluster = cl,
    p = tst$p
  )
}) %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  add_significance("p.adj")

print(stats_tbl)

p <- ggplot(
  df_sample_filt,
  aes(x = factor(cluster), y = score, fill = condition)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.35,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    aes(size = n_cells, color = condition),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    alpha = 0.85
  ) +
  scale_fill_manual(values = c("KO" = "darkgrey", "WT" = "firebrick")) +
  scale_color_manual(values = c("KO" = "darkgrey", "WT" = "firebrick")) +
  theme_classic() +
  labs(
    x = "Seurat cluster",
    y = "Mean module score per sample",
    size = "Cells in cluster"
  )

print(p)
stats_tbl

