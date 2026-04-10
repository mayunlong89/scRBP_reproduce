# Convert Seurat object to h5ad format for Python/Scanpy compatibility

# Load required packages
library(Seurat)
library(SeuratDisk)

# Read the Seurat object
seurat_obj <- readRDS("PBMC_healthy_subset_morethan500cell.rds")

# Copy annotation_broad3 to a standardized column name "celltype"
# Verify the column exists: colnames(seurat_obj@meta.data)
# Convert factor/integer to character to preserve actual cell type names
seurat_obj$annotation_broad3 <- as.character(seurat_obj$annotation_broad3)
seurat_obj$celltype <- seurat_obj$annotation_broad3

# Ensure the default assay is RNA
DefaultAssay(seurat_obj) <- "RNA"

# Strip scale.data to reduce file size (keep counts and normalized data only)
seurat_obj <- DietSeurat(
  seurat_obj,
  counts     = TRUE,
  data       = TRUE,
  scale.data = FALSE,
  dimreducs  = c("pca", "umap")
)

# Save as .h5Seurat format (intermediate format required by SeuratDisk)
SaveH5Seurat(
  seurat_obj,
  filename  = "PBMC_healthy_subset_morethan500cell.h5Seurat",
  assay     = "RNA",
  layer     = "counts",
  overwrite = TRUE
)

# Convert .h5Seurat to .h5ad format (readable by Scanpy)
Convert(
  "PBMC_healthy_subset_morethan500cell.h5Seurat",
  dest      = "h5ad",
  assay     = "RNA",
  overwrite = TRUE,
  verbose   = TRUE,
  options   = list(compression = "gzip", compression_opts = 9)
)
