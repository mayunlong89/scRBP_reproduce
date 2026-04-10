# --- 2025-12-14
# TRS heatmap for 8 autoimmune diseases — validation cohort (10K PBMC, no row filtering)

# ============================================================
# Load required packages
# ============================================================
library(pheatmap)
library(RColorBrewer)

# ============================================================
# Step 1: Read data and split into TRS / P-value matrices
# ============================================================
heatmap_8autoimmunediseases <- read.csv(
  "01_heatmap_8autoimmunediseases_10K_validation.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Set row names from the ID column
rownames(heatmap_8autoimmunediseases) <- heatmap_8autoimmunediseases$ID
data_h1 <- heatmap_8autoimmunediseases[, -1, drop = FALSE]

# Auto-detect TRS and P-value columns by suffix pattern (robust to column reordering)
trs_cols <- grep("_TRS$", colnames(data_h1), value = TRUE)
p_cols   <- grep("_P$",   colnames(data_h1), value = TRUE)

data_h1_score <- data_h1[, trs_cols, drop = FALSE]
data_h1_p     <- data_h1[, p_cols,   drop = FALSE]

# Standardize column names: TRS -> disease name; P -> disease_P
colnames(data_h1_score) <- sub("_TRS$", "", colnames(data_h1_score))
colnames(data_h1_p)     <- paste0(sub("_P$", "", colnames(data_h1_p)), "_P")

# Ensure disease order is consistent between TRS and P matrices
common_traits <- intersect(colnames(data_h1_score), sub("_P$", "", colnames(data_h1_p)))
data_h1_score <- data_h1_score[, common_traits, drop = FALSE]
data_h1_p     <- data_h1_p[, paste0(common_traits, "_P"), drop = FALSE]

# ============================================================
# Step 2: Validation set — use all rows without significance filtering
# ============================================================
heatmap_data <- data_h1_score
p_vals       <- data_h1_p

# ============================================================
# Step 3: Extract cell type from row IDs and define block order
# ============================================================
id_vec        <- rownames(heatmap_data)
cell_type_raw <- sub(".*_", "", id_vec)

# Standardize cell type labels if needed
cell_type <- cell_type_raw
cell_type[cell_type == "mDCs"] <- "mDC"
cell_type[cell_type == "pDCs"] <- "pDC"

cell_type_levels <- c(
  "CD14+ Monocytes", "CD16+ Monocytes",
  "mDC", "pDC", "NK cells",
  "CD4+ T cells", "CD8+ T cells", "B cells"
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

# Reorder all matrices by computed row order
heatmap_data_ord <- heatmap_data[row_order, , drop = FALSE]
p_vals_ord       <- p_vals[row_order, , drop = FALSE]
cell_type_ord    <- cell_type_fac[row_order]

# Build row annotation and gap positions
annotation_row <- data.frame(CellType = cell_type_ord)
rownames(annotation_row) <- rownames(heatmap_data_ord)

ct_cols     <- setNames(brewer.pal(length(cell_type_levels), "Set3"), cell_type_levels)
anno_colors <- list(CellType = ct_cols)

ct_tab   <- table(cell_type_ord)
gaps_row <- cumsum(ct_tab)[-length(ct_tab)]

# ============================================================
# Step 5: Draw heatmap (rows = regulons, columns = diseases)
# ============================================================
reds <- colorRampPalette(rev(c(
  "firebrick3", "#E77A77", "#ECBA84", "white", "lightblue", "#336699"
)))(102)

# Significance symbols: ** P < 0.001, * P < 0.01, # P < 0.05
sig_symbols <- matrix(
  ifelse(p_vals_ord <= 0.001, "**",
         ifelse(p_vals_ord <= 0.01, "*",
                ifelse(p_vals_ord <= 0.05, "#", ""))),
  nrow = nrow(p_vals_ord), byrow = FALSE
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

# ============================================================
# Step 6: Transposed heatmap (rows = diseases, columns = regulons)
# ============================================================
heatmap_data_t <- t(heatmap_data_ord)
p_vals_t       <- t(p_vals_ord)

sig_symbols_t <- matrix(
  ifelse(p_vals_t <= 0.001, "**",
         ifelse(p_vals_t <= 0.01, "*",
                ifelse(p_vals_t <= 0.05, "#", ""))),
  nrow = nrow(p_vals_t), byrow = FALSE
)

# Row annotation becomes column annotation after transpose
pheatmap(
  heatmap_data_t,
  cellwidth         = 5,
  cellheight        = 5,
  fontsize          = 6,
  fontsize_row      = 6,
  border_color      = "black",
  scale             = "none",
  color             = reds,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  annotation_col    = annotation_row,
  annotation_colors = anno_colors,
  gaps_col          = gaps_row,
  display_numbers   = sig_symbols_t
)
