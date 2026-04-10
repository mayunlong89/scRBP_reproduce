# --- 2025-12-20 / 2025-12-23
# TRS UMAP FeaturePlot and cell-type boxplot for significant RBPs
# across 8 autoimmune diseases and 4 mRNA regions (170K PBMC cells)

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
SEURAT_RDS <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_covid19_PBMC/PBMC_healthy_subset_morethan500cell.rds"
TRS_BASE   <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_covid19_PBMC/z_GRNBoost2_15times_results/z_scRBP_results_for_8autoimmune_diseases"
CELLTYPE_COL <- "annotation_broad3"

# ============================================================
# Define all disease x region x RBP combinations to plot
# ============================================================
# Each row: disease, disease_dir_id, region, vector of RBPs
plot_config <- list(
  # --- RA ---
  list(disease = "RA",  dir_id = "03_RA",  region = "3UTR",    rbps = c("ZFP36L1", "RPS11", "PCBP1", "DDX3X", "ZFP36")),
  list(disease = "RA",  dir_id = "03_RA",  region = "5UTR",    rbps = c("MAP1LC3B", "RPS11", "PRPF8", "HNRNPF", "EIF4G2")),
  list(disease = "RA",  dir_id = "03_RA",  region = "CDS",     rbps = c("SSB", "WTAP", "GRSF1", "MKRN1", "MAP1LC3B", "DDX3X")),
  list(disease = "RA",  dir_id = "03_RA",  region = "Introns",  rbps = c("YWHAE", "DDX3X", "PCBP1", "HNRNPUL1", "PABPC1", "G3BP1", "NIPBL")),
  # --- MS ---
  list(disease = "MS",  dir_id = "08_MS",  region = "3UTR",    rbps = c("RPS11", "ZFP36L1", "HNRNPUL1", "ZFP36", "PCBP1", "DDX3X")),
  list(disease = "MS",  dir_id = "08_MS",  region = "5UTR",    rbps = c("RPS15A", "RPS11")),
  list(disease = "MS",  dir_id = "08_MS",  region = "CDS",     rbps = c("WTAP", "MAP1LC3B")),
  list(disease = "MS",  dir_id = "08_MS",  region = "Introns",  rbps = c("YWHAE", "DDX3X", "PCBP1", "NIPBL", "HNRNPUL1")),
  # --- IBD ---
  list(disease = "IBD", dir_id = "01_IBD", region = "3UTR",    rbps = c("RBM17", "PCBP1", "ZFP36")),
  list(disease = "IBD", dir_id = "01_IBD", region = "5UTR",    rbps = c("APOBEC3C", "PUM1", "EIF4G2", "RPS15A")),
  list(disease = "IBD", dir_id = "01_IBD", region = "CDS",     rbps = c("PTBP1")),
  list(disease = "IBD", dir_id = "01_IBD", region = "Introns",  rbps = c("PCBP1")),
  # --- UC ---
  list(disease = "UC",  dir_id = "05_UC",  region = "3UTR",    rbps = c("HNRNPUL1", "SBDS", "PCBP1", "ZFP36")),
  list(disease = "UC",  dir_id = "05_UC",  region = "5UTR",    rbps = c("EIF4G2", "HNRNPM", "APOBEC3C", "PUM1", "RPS15A")),
  list(disease = "UC",  dir_id = "05_UC",  region = "CDS",     rbps = c("WTAP")),
  list(disease = "UC",  dir_id = "05_UC",  region = "Introns",  rbps = c("RPS11", "HNRNPC", "DDX3X")),
  # --- CD ---
  list(disease = "CD",  dir_id = "07_CD",  region = "3UTR",    rbps = c("ZFP36L1", "HNRNPUL1", "RBM17", "PCBP1", "ZFP36")),
  list(disease = "CD",  dir_id = "07_CD",  region = "5UTR",    rbps = c("MAP1LC3B", "HNRNPF", "RPS15A", "PUM1", "RPS11")),
  list(disease = "CD",  dir_id = "07_CD",  region = "Introns",  rbps = c("DDX3X", "PCBP1", "HNRNPUL1", "RPS11")),
  # --- T1D ---
  list(disease = "T1D", dir_id = "06_T1D", region = "3UTR",    rbps = c("DDX3X", "PCBP1", "NOP56", "ZFP36L1", "RPS11")),
  list(disease = "T1D", dir_id = "06_T1D", region = "5UTR",    rbps = c("HNRNPF", "SF3B1", "ZNF800", "TOP1", "RPS11")),
  list(disease = "T1D", dir_id = "06_T1D", region = "CDS",     rbps = c("MAP1LC3B", "MKRN1", "VIM", "DDX3X", "WTAP")),
  list(disease = "T1D", dir_id = "06_T1D", region = "Introns",  rbps = c("HNRNPUL1")),
  # --- PBC ---
  list(disease = "PBC", dir_id = "02_PBC", region = "3UTR",    rbps = c("ZFP36L1")),
  list(disease = "PBC", dir_id = "02_PBC", region = "5UTR",    rbps = c("SF3B1", "LIN28B")),
  # --- SLE ---
  list(disease = "SLE", dir_id = "04_SLE", region = "5UTR",    rbps = c("ZNF800")),
  list(disease = "SLE", dir_id = "04_SLE", region = "CDS",     rbps = c("WTAP", "MAP1LC3B"))
)

# ============================================================
# Reusable function: UMAP FeaturePlot for one RBP TRS
# ============================================================
plot_trs_umap <- function(seurat_obj, trs_matrix, rbp,
                          disease, region,
                          out_dir  = ".",
                          pt.size  = 0.3,
                          pal = colorRampPalette(c("grey85", "white", "#F87C56"))(256)) {

  # Subset to shared cells
  common <- intersect(Cells(seurat_obj), rownames(trs_matrix))
  obj <- subset(seurat_obj, cells = common)

  # Extract TRS values with quantile winsorization (1st–99th percentile)
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
    reduction = "umap",
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

  # Quantile winsorization (1st–99th percentile)
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
    paste0(disease, "_processed_magma_results_170Kcells_", region, "_REAL_only.sc.TRS_matrix.csv")
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
