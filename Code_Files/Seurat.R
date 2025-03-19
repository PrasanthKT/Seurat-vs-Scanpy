####################################################
# scRNA-Seq Analysis with Seurat
####################################################

## Load the Libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)

## Path to the Directory
setwd("~/Desktop/Multi-omics_Project/Multiome")

## Load the scRNA_Seq data
counts <- Read10X_h5("filtered_feature_bc_matrix.h5")
names(counts)

## Extract the gene expression matrix from the counts object
gene_counts <- counts$`Gene Expression`

## Create a Seurat object using the gene expression counts
seurat_obj <- CreateSeuratObject(
  counts = gene_counts, 
  project = "Multiome_RNA"
)

## Summary of the Seurat object
print(seurat_obj)

# Calculate the percentage of mitochondrial genes per cell
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

## Visualize quality control (QC) metrics:
# - nFeature_RNA: number of genes detected per cell
# - nCount_RNA: total UMI counts per cell
# - percent.mt: percentage of mitochondrial gene counts
VlnPlot(
  object = seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

# Summary statistics of key QC metrics
summary(seurat_obj$nFeature_RNA)
summary(seurat_obj$nCount_RNA)
summary(seurat_obj$percent.mt)

## Define quality thresholds:
# - min_genes: Minimum number of genes per cell (to exclude low-quality or empty droplets)
# - max_genes: Maximum number of genes (to remove doublets or multiplets)
# - mt_cutoff: Maximum percentage of mitochondrial genes.
min_genes <- 200
max_genes <- 6000
mt_cutoff <- 15

## Filter the Seurat object based on these QC criteria
seurat_obj <- subset(
  x = seurat_obj,
  subset = nFeature_RNA > min_genes &
    nFeature_RNA < max_genes &
    percent.mt < mt_cutoff
)

## Summaries after Quality filtering
dim(seurat_obj)
summary(seurat_obj$nFeature_RNA)
summary(seurat_obj$nCount_RNA)
summary(seurat_obj$percent.mt)

## QC metrics after filtering
VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

## Normalization
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

## Variable Feature Selection

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 4000,
  verbose = FALSE
)

## Top 20 variable features
top20 <- head(VariableFeatures(seurat_obj), 20)
print(top20)

## Scaling the Data and PCA
seurat_obj <- ScaleData(
  seurat_obj,
  features = rownames(seurat_obj),
  verbose = FALSE
)
## Run PCA on the variable features
seurat_obj <- RunPCA(
  seurat_obj, 
  features = VariableFeatures(seurat_obj),
  npcs = 30,
  verbose = FALSE
)
# ElbowPlot
ElbowPlot(seurat_obj, ndims = 30)

# Set the number of principal components.
dims_use <- 1:12

# Clustering and UMAP
seurat_obj <- FindNeighbors(
  seurat_obj,
  dims = dims_use
)
# 2. Cluster cells 
seurat_obj <- FindClusters(
  seurat_obj,
  resolution = 1.0,
  algorithm = 4,
  random.seed = 42
)
# 3. Run UMAP 
seurat_obj <- RunUMAP(
  seurat_obj,
  dims = dims_use
)
# UMAP-Plot
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)

## Marker Identification
markers <- FindAllMarkers(
  object = seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# View the top markers
head(markers)

# Identify markers excluding ribosomal genes:
markers_no_ribo <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  assay = "RNA",
  features = setdiff(
    rownames(seurat_obj),
    grep("^RP[SL]", rownames(seurat_obj), value = TRUE)
  )
)

# Visualize expression of specific marker genes using FeaturePlot
FeaturePlot(seurat_obj, features = c("KRT1", "NOTCH3", "ANO1", "LY6D"))

# Rank genes within each cluster based on their average log2 fold-change
markers <- markers %>%
  group_by(cluster) %>%
  mutate(rank_gene = rank(-avg_log2FC, ties.method = "first")) %>%
  ungroup()

# Limit the markers to the top 100 genes per cluster
markers_limited <- markers %>%
  group_by(cluster) %>%
  slice_min(order_by = rank_gene, n = 100) %>%
  ungroup()

# Select the top 5 genes per cluster for labeling
markers_label <- markers_limited %>%
  group_by(cluster) %>%
  slice_min(order_by = rank_gene, n = 5) %>%
  ungroup()

# Plot the top markers using ggplot2
ggplot(markers_limited, aes(x = rank_gene, y = avg_log2FC)) +
  geom_point(aes(color = as.factor(cluster)), size = 1, alpha = 0.5) +
  geom_text_repel(
    data = markers_label,
    aes(label = gene),
    size = 3,
    box.padding = 0.3,
    max.overlaps = Inf
  ) +
  facet_wrap(~ cluster, scales = "free", ncol = 4) +
  labs(
    x = "Ranking (by avg_log2FC)",
    y = "avg_log2FC",
    title = "Cluster Markers (Seurat): Score vs. Rank"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10)
  )

# Extract the top 5 markers per cluster based on avg_log2FC
cluster_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  ungroup()

# Create a dot plot of the top 5 markers per cluster
DotPlot(seurat_obj, features = unique(cluster_markers$gene)) +
  RotatedAxis() +
  ggtitle("Top 5 Markers per Cluster")


# Save the processed RNA-seq Seurat object
saveRDS(seurat_obj, file = "processed_rna_seurat.rds")
