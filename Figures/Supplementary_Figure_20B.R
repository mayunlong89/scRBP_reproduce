# --- 2026-11-13
# Regulon activity trend heatmap across 9 developmental stages
# Rows = regulons (split by trend cluster), columns = 100 grid points (pseudotime)

# ============================================================
# Load required packages
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ============================================================
# Input / Output
# ============================================================
trend_file   <- "RegulonTrends_global_v3.trend01_matrix.csv"
cluster_file <- "RegulonTrends_global_v3.clusters.csv"
out_pdf      <- "RegulonTrends_global_v3.stageBinned.pdf"

# ============================================================
# 1) Read trend matrix: regulon x grid points (values scaled 0-1)
# ============================================================
mat_df <- fread(trend_file, data.table = FALSE, check.names = FALSE)

# Handle first column as row names (compatible with or without "regulon" header)
if (!("regulon" %in% colnames(mat_df))) {
  rn <- mat_df[[1]]
  mat_df <- mat_df[, -1, drop = FALSE]
  rownames(mat_df) <- rn
} else {
  rownames(mat_df) <- mat_df$regulon
  mat_df$regulon <- NULL
}

mat <- as.matrix(mat_df)
storage.mode(mat) <- "numeric"
stopifnot(all(grepl("^g\\d+$", colnames(mat))))

# ============================================================
# 2) Read cluster assignments (for row splitting)
# ============================================================
cl <- fread(cluster_file, data.table = FALSE)
cl$cluster <- as.factor(cl$cluster)
cl <- cl[match(rownames(mat), cl$regulon), ]
stopifnot(all(cl$regulon == rownames(mat)))
row_split <- cl$cluster

# ============================================================
# 3) Build continuous stage annotation across grid points
# ============================================================
n_grid    <- ncol(mat)
g_idx     <- seq_len(n_grid)
stage_cont <- 1 + (g_idx - 1) * (9 - 1) / (n_grid - 1)  # Continuous mapping: 1..9

# Discrete stage bins (for labeling)
stage_bin <- cut(stage_cont,
                 breaks = seq(1, 9, length.out = 10),
                 include.lowest = TRUE,
                 labels = paste0("Stage", 1:9))

# Continuous color gradient for stage annotation bar
stage_col_fun <- colorRamp2(
  c(1, 5, 9),
  c("#3B4CC0", "#a2d5f2", "#FDE725")  # Blue -> light blue -> yellow
)

# Top annotation: continuous color strip representing developmental stage
ha_top <- HeatmapAnnotation(
  Stage = stage_cont,
  col   = list(Stage = stage_col_fun),
  annotation_name_side = "left",
  annotation_legend_param = list(
    Stage = list(title = "Stage", at = 1:9, labels = paste0("Stage", 1:9))
  ),
  height = unit(4, "mm")
)

# Stage boundary positions and labels for annotation
stage_break_pos <- round(seq(1, n_grid, length.out = 9))
stage_labels    <- paste0("Stage", 1:9)

# ============================================================
# 4) Heatmap color scale (scaled activity: 0-1)
# ============================================================
hm_col_fun <- colorRamp2(c(0, 0.5, 1), c("white", "#9BD39B", "#1F7A1F"))

# ============================================================
# 5) Draw heatmap
# ============================================================
pdf(out_pdf, width = 16, height = 7)

ht <- Heatmap(
  mat,
  name               = "trend01",
  col                = hm_col_fun,
  cluster_rows       = TRUE,
  cluster_columns    = FALSE,       # Preserve temporal order
  show_column_names  = FALSE,       # g1..g100 labels too dense to display
  row_split          = row_split,   # Split rows by trend cluster
  top_annotation     = ha_top,
  row_names_gp       = gpar(fontsize = 6),
  heatmap_legend_param = list(title = "Scaled activity", at = c(0, 0.5, 1))
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")

# Overlay stage labels above the annotation bar
decorate_annotation("Stage", {
  for (i in seq_along(stage_break_pos)) {
    x <- stage_break_pos[i] / n_grid
    grid.text(stage_labels[i], x = unit(x, "npc"), y = unit(1.6, "npc"),
              just = "center", gp = gpar(fontsize = 9))
  }
})

dev.off()
cat("Wrote:", out_pdf, "\n")
