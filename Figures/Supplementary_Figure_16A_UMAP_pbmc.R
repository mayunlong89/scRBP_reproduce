# --- 2025-08-07
# PBMC clustering and UMAP visualization (10,393 cells after QC)

# ============================================================
# Load required packages
# ============================================================
library(Seurat)
library(ggplot2)

# ============================================================
# Step 1: Load data and filter out "Other" cell type
# ============================================================
df <- readRDS("pbmc_merged.rds")
DefaultAssay(df) <- "RNA"

# Remove "Other" cell type, retaining 10,393 cells
keep <- setdiff(unique(df$merged_celltype), "Other")
Idents(df) <- df$merged_celltype
df_subset <- subset(df, idents = keep)
table(Idents(df_subset))

# Quick UMAP overview with original annotations
p <- DimPlot(
  df_subset,
  reduction = "umap",
  group.by  = "merged_celltype",
  label     = TRUE,
  repel     = TRUE
)

pdf("pbmc_predict_ct.pdf")
p
dev.off()

# ============================================================
# Step 2: Refine monocyte subtype annotations
# ============================================================
source("annotate_monocyte_subtypes.R")
df_subset <- annotate_monocyte_subtypes(
  df_subset,
  celltype_col   = "merged_celltype",
  monocyte_label = "Monocytes",
  assay          = "RNA",
  delta          = 0.25,
  min_expr       = 0.1,
  out_col        = "merged_celltype_refined"
)

# Reclassify ambiguous monocytes and standardize labels
cell_chr <- as.character(df_subset$merged_celltype_refined)
cell_chr[cell_chr == "Monocytes_ambiguous"] <- "CD14+ Monocytes"
cell_chr[cell_chr == "FCGR3A+ Monocytes"]  <- "CD16+ Monocytes"
df_subset$merged_celltype_refined <- factor(cell_chr)

# ============================================================
# Step 3: Re-cluster with standard Seurat workflow
# ============================================================
pbmc <- df_subset
DefaultAssay(pbmc) <- "RNA"

# Normalize, find variable features, scale
pbmc <- NormalizeData(pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc, nfeatures = 3000, verbose = FALSE)
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc), verbose = FALSE)

# PCA (50 components)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 50, verbose = FALSE)
ElbowPlot(pbmc)

# Neighbor graph and clustering (using top 30 PCs)
pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.4)

# ============================================================
# Step 4: UMAP with cosine distance
# ============================================================
pbmc <- RunUMAP(
  pbmc,
  reduction      = "pca",
  dims           = 1:30,
  n.neighbors    = 40,
  min.dist       = 0.25,
  metric         = "cosine",
  seed.use       = 12345,
  reduction.name = "umap_cosine",
  reduction.key  = "UMAPC_",
  verbose        = FALSE
)

# Ensure monocyte labels are consistent after clustering
cell_chr <- as.character(pbmc$merged_celltype_refined)
cell_chr[cell_chr == "Monocytes_ambiguous"] <- "CD14+ Monocytes"
cell_chr[cell_chr == "FCGR3A+ Monocytes"]  <- "CD16+ Monocytes"
pbmc$merged_celltype_refined <- factor(cell_chr)

# ============================================================
# Step 5: Plot UMAP colored by refined cell type
# ============================================================
ct_colors <- c(
  "CD14+ Monocytes" = "#4A774D",
  "CD16+ Monocytes" = "#8DC890",
  "mDCs"            = "#C3D5DE",
  "pDCs"            = "#A02D28",
  "NK cells"        = "#225EA8",
  "CD4+ T cells"    = "#FB6A4A",
  "CD8+ T cells"    = "#6A51A3",
  "B cells"         = "#E69F00"
)

p <- DimPlot(
  pbmc,
  reduction = "umap_cosine",
  group.by  = "merged_celltype_refined",
  cols      = ct_colors,
  label     = TRUE,
  repel     = TRUE
) + theme(legend.title = element_blank())

pdf("pbmc_predict_ct2.pdf")
p
dev.off()

# ============================================================
# Step 6: Save outputs
# ============================================================
saveRDS(pbmc, file = "PBMC_subset10393cells_reclustered.rds")
ggsave("PBMC_subset_umap_clusters.png",  width = 7, height = 6, dpi = 300)
ggsave("PBMC_subset_umap_celltypes.png", width = 7, height = 6, dpi = 300)
