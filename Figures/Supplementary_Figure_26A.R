# --- 2026-02-10
# TRS UMAP FeaturePlot and cell-type boxplot for significant RBPs
# across 8 brain/psychiatric traits and 4 mRNA regions (590K fetal brain cells)

# ============================================================
# Load required packages
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(grid)
  library(dplyr)
})

# ============================================================
# Paths (edit these for your environment)
# ============================================================
SEURAT_RDS <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds"
TRS_BASE   <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes"
CELLTYPE_COL <- "Type.v2"
UMAP_REDUCTION <- "umap_continuum"

# ============================================================
# Define all disease x region x RBP combinations to plot
# ============================================================
plot_config <- list(
  # --- SCZ (Schizophrenia) ---
  list(disease = "SCZ", dir_id = "07_SCZ", region = "3UTR",
       rbps = c("ELAVL4", "CELF2", "CPEB4", "EIF3B", "RBFOX2", "HNRNPR", "HNRNPU", "SRRM4", "SRSF6", "YBX3")),
  list(disease = "SCZ", dir_id = "07_SCZ", region = "5UTR",
       rbps = c("CELF4", "EIF4A1", "EIF4G2", "RBFOX2", "RBM17", "HNRNPA1", "HNRNPH1")),
  list(disease = "SCZ", dir_id = "07_SCZ", region = "CDS",
       rbps = c("HNRNPH1", "NOVA1", "PABPC1", "RBFOX1", "RPS3", "SRRM4", "SRSF9", "YWHAE", "YWHAH")),
  list(disease = "SCZ", dir_id = "07_SCZ", region = "Introns",
       rbps = c("CELF2", "DDX5", "KHSRP", "RBFOX2", "RBFOX3", "TRA2A")),

  # --- BIP (Bipolar disorder) ---
  list(disease = "BIP", dir_id = "04_BIP", region = "3UTR",
       rbps = c("CELF2", "SRRM4", "KHSRP")),

  # --- ASD (Autism spectrum disorder) ---
  list(disease = "ASD", dir_id = "01_ASD", region = "3UTR",
       rbps = c("CELF2", "SRRM4")),
  list(disease = "ASD", dir_id = "01_ASD", region = "5UTR",
       rbps = c("CELF4", "DDX24", "FTO", "ILF3", "RBMX")),
  list(disease = "ASD", dir_id = "01_ASD", region = "CDS",
       rbps = c("DHX9", "ELAVL3", "HNRNPM")),

  # --- ADHD (Attention deficit hyperactivity disorder) ---
  list(disease = "ADHD", dir_id = "02_ADHD", region = "3UTR",
       rbps = c("CELF2", "SRRM4", "RBFOX2")),
  list(disease = "ADHD", dir_id = "02_ADHD", region = "5UTR",
       rbps = c("CELF4", "HNRNPH1", "RPS19", "EIF4G2", "RBFOX2")),
  list(disease = "ADHD", dir_id = "02_ADHD", region = "Introns",
       rbps = c("CELF2", "RPS27A")),

  # --- TS (Tourette syndrome) ---
  list(disease = "TS", dir_id = "08_TS", region = "3UTR",
       rbps = c("ELAVL3", "ELAVL4", "CPEB4")),

  # --- AN (Anorexia nervosa) ---
  list(disease = "AN", dir_id = "03_AN", region = "5UTR",
       rbps = c("HNRNPA1", "HNRNPH1")),

  # --- MDD (Major depressive disorder) ---
  list(disease = "MDD", dir_id = "05_MDD", region = "5UTR",
       rbps = c("ERI2", "CELF4")),

  # --- OCD (Obsessive-compulsive disorder) ---
  list(disease = "OCD", dir_id = "06_OCD", region = "Introns",
       rbps = c("LSM11", "RBFOX2", "RBFOX3"))
)

# ============================================================
# Reusable function: UMAP FeaturePlot for one RBP TRS
# ============================================================
plot_trs_umap <- function(seurat_obj, trs_matrix, rbp,
                          disease, region,
                          reduction = "umap_continuum",
                          out_dir   = ".",
                          pt.size   = 0.3,
                          pal = colorRampPalette(c("grey85", "white", "#F87C56"))(256)) {

  # Subset to shared cells
  common <- intersect(Cells(seurat_obj), rownames(trs_matrix))
  obj <- subset(seurat_obj, cells = common)

  # Extract TRS values with quantile winsorization (1st-99th percentile)
  v <- trs_matrix[common, rbp, drop = TRUE]
  q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
  v <- pmin(pmax(v, q[1]), q[2])

  # Add TRS to metadata
  feature_name <- paste0(rbp, "_TRS")
  obj[[feature_name]] <- v

  # Generate UMAP FeaturePlot
  p <- FeaturePlot(
    obj,
    features  = feature_name,
    reduction = reduction,
    pt.size   = pt.size,
    order     = TRUE,
    raster    = FALSE
  ) +
    scale_color_gradientn(
      colours = pal,
      name    = "TRS",
      oob     = scales::squish,
      guide   = guide_colorbar(
        nbin      = 256,
        raster    = TRUE,
        barheight = unit(35, "mm"),
        barwidth  = unit(4,  "mm")
      )
    ) +
    theme_void() +
    ggtitle(paste(rbp, "TRS on UMAP"))

  # Save as TIFF (publication-quality)
  outfile <- file.path(out_dir, paste0(rbp, "_", region, "_", disease, ".umap_trs.tiff"))
  ggsave(outfile, p, width = 7, height = 5, units = "in", dpi = 600, compression = "lzw")
  message("[OK] Saved: ", outfile)

  invisible(p)
}

# ============================================================
# Reusable function: TRS boxplot by cell type
# ============================================================
plot_trs_box_by_celltype <- function(seurat_obj, trs_matrix, rbp,
                                     celltype_col,
                                     disease, region,
                                     out_dir       = ".",
                                     pal_cols      = c("#F26B4F", "#F7C4B2", "#FCEFEA"),
                                     sort_by       = "median",
                                     show_outliers = FALSE,
                                     box_width     = 0.65) {

  # Extract TRS values aligned to Seurat cell order
  v <- trs_matrix[Cells(seurat_obj), rbp, drop = TRUE]
  names(v) <- Cells(seurat_obj)

  # Quantile winsorization (1st-99th percentile)
  q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
  v <- pmin(pmax(v, q[1]), q[2])
  names(v) <- Cells(seurat_obj)

  # Build data frame with cell type annotations
  df <- data.frame(
    cell      = names(v),
    TRS       = as.numeric(v),
    cell_type = seurat_obj@meta.data[names(v), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  # Sort cell types by median/mean TRS (descending)
  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # Gradient color palette mapped to sorted cell types
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  p <- ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width         = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size          = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp, " TRS across cell types"),
      x     = NULL,
      y     = "scRBP TRS"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )

  # Save as TIFF
  outfile <- file.path(out_dir, paste0(rbp, "_", region, "_", disease, ".TRS_box_by_celltype.tiff"))
  ggsave(outfile, p, width = 10, height = 4.2, units = "in", dpi = 600, compression = "lzw")
  message("[OK] Saved: ", outfile)

  invisible(p)
}

# ============================================================
# Batch execution: loop over all configurations
# ============================================================
# Load Seurat object once
seurat_obj <- readRDS(SEURAT_RDS)

for (cfg in plot_config) {
  disease <- cfg$disease
  dir_id  <- cfg$dir_id
  region  <- cfg$region
  rbps    <- cfg$rbps

  # Construct TRS matrix path
  trs_file <- file.path(
    TRS_BASE, dir_id, paste0("z_", disease, "_TRS"),
    paste0(disease, "_processed_magma_results_590Kcells_", region, "_REAL_only.sc.TRS_matrix.csv")
  )

  if (!file.exists(trs_file)) {
    warning("[SKIP] TRS file not found: ", trs_file)
    next
  }

  # Output directory (same as TRS directory)
  out_dir <- dirname(trs_file)

  # Read TRS matrix once per disease x region
  trs_matrix <- read.csv(trs_file, row.names = 1, check.names = FALSE)

  # Generate UMAP and boxplot for each RBP
  for (rbp in rbps) {
    if (!rbp %in% colnames(trs_matrix)) {
      warning("[SKIP] RBP not found in TRS matrix: ", rbp, " (", disease, " ", region, ")")
      next
    }

    # UMAP FeaturePlot
    plot_trs_umap(
      seurat_obj = seurat_obj,
      trs_matrix = trs_matrix,
      rbp        = rbp,
      disease    = disease,
      region     = region,
      reduction  = UMAP_REDUCTION,
      out_dir    = out_dir
    )

    # Cell-type boxplot
    plot_trs_box_by_celltype(
      seurat_obj   = seurat_obj,
      trs_matrix   = trs_matrix,
      rbp          = rbp,
      celltype_col = CELLTYPE_COL,
      disease      = disease,
      region       = region,
      out_dir      = out_dir
    )
  }
}
