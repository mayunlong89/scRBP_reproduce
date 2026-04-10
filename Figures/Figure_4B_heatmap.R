# --- 2025-12-13
# TRS heatmap for 8 autoimmune diseases (Figure 4B)
# Rows grouped by cell type with within-block hierarchical clustering

# ============================================================
# Load required packages
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})

# ============================================================
# Step 1: Read data and split into TRS / P-value matrices
# ============================================================
heatmap_8autoimmunediseases <- read.csv(
  "01_heatmap_8autoimmunediseases.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Set row names from the ID column
rownames(heatmap_8autoimmunediseases) <- heatmap_8autoimmunediseases$ID
data_h1 <- heatmap_8autoimmunediseases[, -1]

# TRS score matrix (odd-numbered columns: 1,3,5,...,15)
data_h1_score <- data_h1[, c(1, 3, 5, 7, 9, 11, 13, 15)]
colnames(data_h1_score) <- c("CD", "IBD", "MS", "PBC", "RA", "SLE", "T1D", "UC")

# P-value matrix (even-numbered columns: 2,4,6,...,16)
data_h1_p <- data_h1[, c(2, 4, 6, 8, 10, 12, 14, 16)]
colnames(data_h1_p) <- paste0(colnames(data_h1_score), "_P")

# ============================================================
# Step 2: Filter rows with at least one significant disease (P < 0.05)
# ============================================================
sig_rows <- rowSums(data_h1_p < 0.05) > 0

heatmap_data <- data_h1_score[sig_rows, , drop = FALSE]
p_vals       <- data_h1_p[sig_rows, , drop = FALSE]

# ============================================================
# Step 3: Extract cell type from row IDs and define block order
# ============================================================
id_vec        <- rownames(heatmap_data)
cell_type_raw <- sub(".*_", "", id_vec)

# Map raw cell type labels to display names
cell_type_map <- c(
  "Monocyte CD14" = "CD14+ monocytes",
  "Monocyte CD16" = "CD16+ monocytes",
  "mDC"           = "mDC",
  "pDC"           = "pDC",
  "NK"            = "NK cells",
  "T CD4+"        = "CD4+ T cells",
  "T CD8+"        = "CD8+ T cells",
  "T g/d"         = "T gamma-delta cells",
  "T reg"         = "T reg cells",
  "Cycling"       = "Cycling",
  "MAIT"          = "MAIT",
  "Platelets"     = "Platelets & B cells",
  "B"             = "Platelets & B cells"
)

cell_type <- ifelse(cell_type_raw %in% names(cell_type_map),
                    cell_type_map[cell_type_raw],
                    cell_type_raw)

# Define display order for cell type blocks
cell_type_levels <- c(
  "CD14+ monocytes", "CD16+ monocytes", "mDC", "pDC",
  "NK cells", "CD4+ T cells", "CD8+ T cells",
  "T gamma-delta cells", "T reg cells",
  "Cycling", "MAIT", "Platelets & B cells"
)

cell_type_fac <- factor(cell_type, levels = cell_type_levels)

# ============================================================
# Step 4: Hierarchical clustering within each cell type block
# ============================================================
split_idx <- split(seq_len(nrow(heatmap_data)), cell_type_fac)

row_order <- unlist(lapply(split_idx, function(idx) {
  idx <- idx[!is.na(idx)]
  if (length(idx) <= 1) return(idx)
  submat <- heatmap_data[idx, , drop = FALSE]
  hc <- hclust(dist(submat))
  idx[hc$order]
}))

# Reorder all matrices by the computed row order
heatmap_data_ord <- heatmap_data[row_order, , drop = FALSE]
p_vals_ord       <- p_vals[row_order, , drop = FALSE]
cell_type_ord    <- cell_type_fac[row_order]

# ============================================================
# Step 5: Build row annotation (cell type color bar) and block gaps
# ============================================================
annotation_row <- data.frame(CellType = cell_type_ord)
rownames(annotation_row) <- rownames(heatmap_data_ord)

# Assign a distinct color to each cell type
ct_cols <- setNames(
  brewer.pal(length(cell_type_levels), "Set3"),
  cell_type_levels
)
anno_colors <- list(CellType = ct_cols)

# Compute gap positions between cell type blocks
ct_tab   <- table(cell_type_ord)
gaps_row <- cumsum(ct_tab)[-length(ct_tab)]

# ============================================================
# Step 6: Draw heatmap with significance annotations
# ============================================================
# Red-blue diverging color palette
reds <- colorRampPalette(rev(c(
  "firebrick3", "#E77A77", "#ECBA84", "white", "lightblue", "#336699"
)))(102)

# Significance symbols: ** P < 0.001, * P < 0.01, # P < 0.05
sig_symbols <- matrix(
  ifelse(p_vals_ord < 0.001, "**",
         ifelse(p_vals_ord < 0.01, "*",
                ifelse(p_vals_ord <= 0.05, "#", ""))),
  nrow = nrow(p_vals_ord),
  byrow = FALSE
)

pearson_heatmap <- pheatmap(
  heatmap_data_ord,
  cellwidth         = 5,
  cellheight        = 5,
  fontsize          = 6,
  fontsize_row      = 6,
  border_color      = "black",
  scale             = "none",
  color             = reds,
  cluster_rows      = FALSE,
  cluster_cols      = TRUE,
  annotation_row    = annotation_row,
  annotation_colors = anno_colors,
  gaps_row          = gaps_row,
  display_numbers   = sig_symbols
)

pearson_heatmap
