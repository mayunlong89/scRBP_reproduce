
#2026-02-22
# ====== Packages ======
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

# =========================
# 0) Input
# =========================

#version 3 ---> 217 regulon-celltype pairs
#infile <- "01_heatmap_8braindiseases.csv"  # <- 改成你的文件路径

infile <- "01_heatmap_9metabolic_traits_Heatmap_TRS.csv"  # <- 改成你的文件路径

#version 1 ---> 128 regulon-celltype pairs
#infile <- "01_heatmap_8braindiseases_v1_high_priority.csv"  # <- 改成你的文件路径

heatmap_df <- read.csv(infile, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

stopifnot("ID" %in% colnames(heatmap_df))
rownames(heatmap_df) <- heatmap_df$ID
data_h1 <- heatmap_df %>% dplyr::select(-ID)

# =========================
# 1) Auto-detect TRS / P columns
#    Expect pairs like: SCZ_TRS, SCZ_P, ...
# =========================
trs_cols <- grep("_TRS$", colnames(data_h1), value = TRUE)
p_cols   <- grep("_P$",   colnames(data_h1), value = TRUE)

if (length(trs_cols) == 0 || length(p_cols) == 0) {
  stop("No *_TRS or *_P columns detected. Please check column names.")
}

# traits inferred from TRS columns
traits <- sub("_TRS$", "", trs_cols)

# keep only traits with both TRS and P
keep_traits <- traits[paste0(traits, "_P") %in% p_cols]
if (length(keep_traits) == 0) stop("No matched TRS/P pairs found.")

# Order traits (use your preferred order when present)
preferred_trait_order <- c("AKP","ALT","HDL","LDL","SHBG","TBIL","TC","TG","TST")
trait_order <- c(preferred_trait_order[preferred_trait_order %in% keep_traits],
                 setdiff(keep_traits, preferred_trait_order))

data_h1_score <- data_h1 %>% dplyr::select(all_of(paste0(trait_order, "_TRS")))
data_h1_p     <- data_h1 %>% dplyr::select(all_of(paste0(trait_order, "_P")))

colnames(data_h1_score) <- trait_order
colnames(data_h1_p)     <- paste0(trait_order, "_P")

# =========================
# 2) Keep rows with any significant disease (P < 0.05)
# =========================
sig_rows <- apply(data_h1_p, 1, function(x) any(x < 0.05, na.rm = TRUE))
heatmap_data <- data_h1_score[sig_rows, , drop = FALSE]
p_vals       <- data_h1_p[sig_rows, , drop = FALSE]

if (nrow(heatmap_data) == 0) stop("No rows pass significance filter (any P < 0.05).")

# =========================
# 3) Parse cell type from ID
#    Default: take substring after the last "_"
# =========================
id_vec <- rownames(heatmap_data)

cell_type_raw <- sub(".*_", "", id_vec)
cell_type_raw <- trimws(cell_type_raw)

# (可选) 如果你们新数据有需要统一命名，在这里加映射规则（示例保留）

cell_type <- cell_type_raw


# =========================
# 4) Define cell type block order
#    Use your biological order; append unseen types at the end
# =========================
# =========================
# 4) Define cell type order
#    Order by regulon counts per cell type (descending),
#    BUT keep hepatocyte always at the top.
# =========================

# counts are based on the currently retained (significant) rows
celltype_counts <- table(cell_type)

# sort: primary = -count, secondary = cell type name (stable)
celltype_sorted <- names(sort(celltype_counts, decreasing = TRUE))

# force hepatocyte first (if present)
if ("hepatocyte" %in% celltype_sorted) {
  celltype_sorted <- c("hepatocyte", setdiff(celltype_sorted, "hepatocyte"))
} else {
  warning("Cell type 'hepatocyte' not found in current heatmap rows.")
}

# set factor levels using this order
cell_type_levels <- celltype_sorted
cell_type_fac <- factor(cell_type, levels = cell_type_levels)




preferred_celltype_order <-  celltype_sorted

cell_type_levels <- c(preferred_celltype_order[preferred_celltype_order %in% unique(cell_type)],
                      setdiff(unique(cell_type), preferred_celltype_order))

cell_type_fac <- factor(cell_type, levels = cell_type_levels)



# =========================
# 5) Row order: cluster within each cell type block
# =========================
split_idx <- split(seq_len(nrow(heatmap_data)), cell_type_fac)

row_order <- unlist(lapply(split_idx, function(idx) {
  idx <- idx[!is.na(idx)]
  if (length(idx) <= 1) return(idx)
  
  submat <- as.matrix(heatmap_data[idx, , drop = FALSE])
  
  # correlation distance within block (more stable than Euclidean for TRS)
  # fall back to Euclidean if cor fails
  ord <- tryCatch({
    cm <- suppressWarnings(cor(t(submat), method = "pearson", use = "pairwise.complete.obs"))
    d  <- as.dist(1 - cm)
    hc <- hclust(d, method = "average")
    idx[hc$order]
  }, error = function(e) {
    d  <- dist(submat)
    hc <- hclust(d, method = "average")
    idx[hc$order]
  })
  ord
}))

heatmap_data_ord <- heatmap_data[row_order, , drop = FALSE]
p_vals_ord       <- p_vals[row_order, , drop = FALSE]
cell_type_ord    <- cell_type_fac[row_order]

# =========================
# 6) Annotation + gaps
# =========================
annotation_row <- data.frame(CellType = cell_type_ord)
rownames(annotation_row) <- rownames(heatmap_data_ord)

# More colors: use hcl.colors (works for many categories)
ct_cols <- setNames(grDevices::hcl.colors(length(cell_type_levels), palette = "Dynamic"),
                    cell_type_levels)
anno_colors <- list(CellType = ct_cols)

ct_tab   <- table(cell_type_ord)
gaps_row <- cumsum(ct_tab)[-length(ct_tab)]  # block separators

# =========================
# 7) Heatmap colors + significance symbols
# =========================
reds <- colorRampPalette(rev(c("firebrick3","firebrick3","#E77A77","#E77A77","#E77A77","#ECBA84","#ECBA84","white","lightblue","#336699")))(102)

sig_mat <- matrix(
  ifelse(p_vals_ord <= 0.001, "**",
         ifelse(p_vals_ord <= 0.01, "*",
                ifelse(p_vals_ord <= 0.05, "#", ""))),
  nrow = nrow(p_vals_ord),
  ncol = ncol(p_vals_ord),
  byrow = FALSE
)

# =========================
# 8) Plot
# =========================
pheatmap(
  mat = as.matrix(heatmap_data_ord),
  color = reds,
  scale = "none",
  border_color = "black",
  fontsize = 6,
  fontsize_row = 6,
  cellwidth = 5,
  cellheight = 5,
  cluster_rows = FALSE,   # 我们已手动排好 row order
  cluster_cols = TRUE,    # 疾病列聚类（如果要固定顺序改 FALSE）
  annotation_row = annotation_row,
  annotation_colors = anno_colors,
  gaps_row = gaps_row,
  display_numbers = sig_mat,
  number_color = "black"
)
