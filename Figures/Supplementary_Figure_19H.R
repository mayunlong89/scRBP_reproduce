
#2026-01-09

#Step 2.3: cell-subtype-by-regulon activity score--heatmap plotting

suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
  library(dplyr)
})

# =========================
# User knobs
# =========================

# regulon  cell subtype matrix
infile  <- "10X_fetal_brain_590Kcells_ras_ct_all4regions.symbol_cell_subtype_rss.csv"
#infile  <- "fetal_brain_scRBP_ras_RSS_mean_20reps.csv"
outfile <- "fetal_brain_590K_regulon_specificity_RSS_heatmap_cell_subtype_rss.pdf"


#RBP mean expression
infile  <- "celltype_by_RBP_mean_expression.csv"
outfile <- "fetal_brain_590K_RBP_mean_expression_heatmap.pdf"

# mode:
#   "raw" : plot the input RSS/JSD values directly (recommended if already RSS)
#   "mss" : normalize each regulon across cell types to [0,1] by dividing by max (Figure-1-like)
mode <- "mss"

# optional: keep top-variance regulons (set NA to keep all)
top_k_regulons_by_var <- NA

# clustering
cluster_rows <- TRUE
cluster_cols <- TRUE
dist_method  <- "euclidean"
clust_method <- "complete"

# display settings
show_colnames_if_leq <- 80

# color control
# for raw RSS, you may cap at a high percentile to avoid a few extreme values dominating
cap_quantile <- 0.995   # set NULL to disable, e.g., NULL
# =========================

# -------------------------
# Read
# -------------------------
df <- fread(infile, data.table = FALSE, check.names = FALSE)
stopifnot(ncol(df) >= 2)

rownames(df) <- df[[1]]
df <- df[, -1, drop = FALSE]

mat0 <- as.matrix(df)
storage.mode(mat0) <- "numeric"

# Basic sanity: RSS/JSD should be non-negative in most implementations
# (If your file includes z-scores, negatives can appear; then use mode="colZ" instead.)
# message("Range: ", paste(range(mat0, na.rm = TRUE), collapse = " ~ "))

# Optional: keep top-variance regulons (columns)
mat <- mat0
if (!is.na(top_k_regulons_by_var) && top_k_regulons_by_var < ncol(mat)) {
  v <- apply(mat, 2, stats::var, na.rm = TRUE)
  keep <- names(sort(v, decreasing = TRUE))[seq_len(top_k_regulons_by_var)]
  mat <- mat[, keep, drop = FALSE]
}

# -------------------------
# Transform for visualization (NO z-score scaling for RSS/JSD)
# -------------------------
if (mode == "mss") {
  # Normalize each regulon to [0,1] across cell types (Figure-1-like: 0 -> subtype max)
  # If a regulon is all zeros, keep it as zeros.
  m <- apply(mat, 2, max, na.rm = TRUE)
  m[m == 0 | is.na(m)] <- 1
  mat <- sweep(mat, 2, m, "/")
}

# Optional cap for raw/mss to stabilize colors (winsorize)
if (!is.null(cap_quantile) && is.finite(cap_quantile)) {
  cap <- stats::quantile(mat, probs = cap_quantile, na.rm = TRUE, names = FALSE)
  mat[mat > cap] <- cap
}

# Replace NA (if any) with 0 for plotting
mat[is.na(mat)] <- 0

# -------------------------
# Row annotation: lineage from cell type names (edit patterns as needed)
# -------------------------
celltype <- rownames(mat)
lineage <- dplyr::case_when(
  grepl("^RG$|RG\\.div|Radial", celltype, ignore.case = TRUE) ~ "RG",
  grepl("Div|Prog|IPC", celltype, ignore.case = TRUE)        ~ "Progenitor/IPC",
  grepl("OPC", celltype, ignore.case = TRUE)                 ~ "OPC",
  grepl("Oligo", celltype, ignore.case = TRUE)               ~ "Oligodendrocyte",
  grepl("Astro", celltype, ignore.case = TRUE)               ~ "Astrocyte",
  grepl("Micro", celltype, ignore.case = TRUE)               ~ "Microglia",
  grepl("Endo|Peri|Vasc", celltype, ignore.case = TRUE)      ~ "Vascular",
  grepl("Excit|Inhib|Neuron|Newborn", celltype, ignore.case = TRUE) ~ "Neuron",
  TRUE                                                       ~ "Other"
)

anno_row <- data.frame(Lineage = factor(lineage))
rownames(anno_row) <- rownames(mat)

#2171b5
# -------------------------
# Palette
# -------------------------
if (mode == "mss") {
  # 0 -> max (white -> blue) similar to many "module specificity score" plots
  pal <- colorRampPalette(c("white","#466812"))(256)
  bk  <- seq(0, max(mat, na.rm = TRUE), length.out = 257)
} else {
  # raw RSS/JSD: white -> blue (higher = more specific)
  pal <- colorRampPalette(c("white", "#466812"))(256)
  bk  <- seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = 257)
}

# -------------------------
# Plot
# -------------------------
if (grepl("\\.pdf$", outfile, ignore.case = TRUE)) {
  pdf(outfile, width = 12, height = 7)
  on.exit(dev.off(), add = TRUE)
} else if (grepl("\\.png$", outfile, ignore.case = TRUE)) {
  png(outfile, width = 2400, height = 1400, res = 200)
  on.exit(dev.off(), add = TRUE)
}


pheatmap(
  mat,
  color = pal,
  breaks = bk,
  border_color = NA,
  scale = "none",   # IMPORTANT: do not rescale RSS/JSD again
  
  cluster_rows = cluster_rows,
  cluster_cols = cluster_cols,
  clustering_distance_rows = dist_method,
  clustering_distance_cols = dist_method,
  clustering_method = clust_method,
  
  annotation_row = anno_row,
  
  show_rownames = TRUE,
  show_colnames = (ncol(mat) <= show_colnames_if_leq),
  fontsize_row = 9,
  fontsize_col = 6
)



lineage <- dplyr::case_when(
  grepl("^RG$|RG\\.div|Radial", celltype, ignore.case = TRUE) ~ "RG",
  grepl("Div|Prog|IPC", celltype, ignore.case = TRUE)        ~ "Progenitor/IPC",
  grepl("OPC", celltype, ignore.case = TRUE)                 ~ "OPC",
  grepl("Oligo", celltype, ignore.case = TRUE)               ~ "Oligodendrocyte",
  grepl("Astro", celltype, ignore.case = TRUE)               ~ "Astrocyte",
  grepl("Micro", celltype, ignore.case = TRUE)               ~ "Microglia",
  grepl("Endo|Peri|Vasc", celltype, ignore.case = TRUE)      ~ "Vascular",
  grepl("Excit|Inhib|Neuron|Newborn", celltype, ignore.case = TRUE) ~ "Neuronal",
  TRUE                                                       ~ "Other"
)


stage_col_fun <- colorRamp2(
  c(1, 5, 9),
  c("#3B4CC0", "#FDE725", "#B40426")  # 蓝->黄->红（示例）
)

# =========================
# 1) For each cell type, pick top-K "most specific" regulons
#    specificity = x_ct - max(other_ct)
# =========================
top_k <- 30

mat_use <- mat   # use the matrix you plotted (after mode="mss" normalization)
stopifnot(nrow(mat_use) >= 2)

cell_types <- rownames(mat_use)
regulons   <- colnames(mat_use)

# Compute delta matrix: rows=celltypes, cols=regulons
# delta[i,j] = mat[i,j] - max(mat[-i, j])
delta <- matrix(NA_real_, nrow = nrow(mat_use), ncol = ncol(mat_use),
                dimnames = dimnames(mat_use))

for (i in seq_len(nrow(mat_use))) {
  others_max <- apply(mat_use[-i, , drop = FALSE], 2, max, na.rm = TRUE)
  delta[i, ] <- mat_use[i, ] - others_max
}

# For each cell type, select top-K regulons by delta (largest advantage)
top_list <- lapply(cell_types, function(ct) {
  d <- delta[ct, ]
  # optional: require that ct is also the absolute max cell type for that regulon
  # is_winner <- mat_use[ct, ] == apply(mat_use, 2, max, na.rm = TRUE)
  # d[!is_winner] <- -Inf
  
  ord <- order(d, decreasing = TRUE)
  keep <- ord[seq_len(min(top_k, length(ord)))]
  data.frame(
    cell_type = ct,
    regulon   = regulons[keep],
    score_ct  = mat_use[ct, keep],
    delta     = d[keep],
    stringsAsFactors = FALSE
  )
})

top_tab <- dplyr::bind_rows(top_list) %>%
  dplyr::arrange(cell_type, dplyr::desc(delta))

# Save the table for inspection
data.table::fwrite(top_tab, file = "celltype_top10_regulons_by_delta.csv")

# =========================
# 2) Build the union regulon set and make a subset heatmap matrix
# =========================
regs_union <- unique(top_tab$regulon)

mat_sub <- mat_use[, regs_union, drop = FALSE]

# Option A (recommended): keep your "mss" semantics (0-1), but reorder columns by winner cell type
winner_ct <- apply(mat_sub, 2, function(x) rownames(mat_sub)[which.max(x)])
winner_ct <- factor(winner_ct, levels = cell_types)

ord_cols <- order(winner_ct, apply(mat_sub, 2, max, na.rm = TRUE), decreasing = TRUE)
mat_sub  <- mat_sub[, ord_cols, drop = FALSE]
winner_ct <- winner_ct[ord_cols]

anno_col <- data.frame(WinnerCellType = winner_ct)
rownames(anno_col) <- colnames(mat_sub)

# =========================
# 3) Plot subset heatmap
# =========================
# white -> blue, for mss in [0,1]
#pal <- colorRampPalette(c("white", "#2171b5"))(256)
#pal <- colorRampPalette(c("#466812","white", "#FF837D"))(256)
pal <- colorRampPalette(c("white", "#466812"))(256)


bk  <- seq(0, 1, length.out = 257)

pheatmap::pheatmap(
  mat_sub,
  color = pal,
  breaks = bk,
  border_color = NA,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,         # keep our ordering by winner cell type (more interpretable)
  annotation_row = anno_row,
  annotation_col = anno_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 9,
  fontsize_col = 3
)

 

# Option B (if you prefer clustering within the subset):
# set cluster_cols = TRUE and remove the manual ordering above.



write.csv(mat_sub,file="celltype_by_RBP_mean_expression_top.csv", quote = F)
write.csv(mat_sub,file="10X_fetal_brain_590Kcells_ras_ct_all4regions.symbol_main_celltype_rss_TOP.csv", quote = F)



