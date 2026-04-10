# --- 2026-01-09
# Cell-type-by-regulon activity score (RSS/MSS) heatmap
# Identifies top cell-type-specific regulons and generates subset heatmaps

# ============================================================
# Load required packages
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
  library(dplyr)
})

# ============================================================
# User configuration
# ============================================================
# Input: cell type x regulon specificity matrix (RSS or mean expression)
# Uncomment ONE input/output pair:

# Option A: Main cell type RSS matrix
infile  <- "10X_fetal_brain_590Kcells_ras_ct_all4regions.symbol_main_celltype_rss.csv"
outfile <- "fetal_brain_590K_regulon_specificity_RSS_heatmap.pdf"

# Option B: Cell subtype RSS matrix
# infile  <- "10X_fetal_brain_590Kcells_ras_ct_all4regions.symbol_cell_subtype_rss.csv"
# outfile <- "fetal_brain_590K_regulon_specificity_RSS_heatmap_cell_subtype_rss.pdf"

# Option C: RBP mean expression matrix
# infile  <- "celltype_by_RBP_mean_expression.csv"
# outfile <- "fetal_brain_590K_RBP_mean_expression_heatmap.pdf"

# Normalization mode:
#   "raw": plot input RSS/JSD values directly (recommended if already RSS)
#   "mss": normalize each regulon to [0,1] across cell types (module specificity score)
mode <- "mss"

# Optional: keep only top-variance regulons (set NA to keep all)
top_k_regulons_by_var <- NA

# Clustering parameters
cluster_rows <- TRUE
cluster_cols <- TRUE
dist_method  <- "euclidean"
clust_method <- "complete"

# Display settings
show_colnames_if_leq <- 80

# Color winsorization: cap at this quantile to avoid extreme values dominating
cap_quantile <- 0.995  # Set NULL to disable

# ============================================================
# Part 1: Full heatmap of all regulons
# ============================================================

# Read input matrix (rows = cell types, columns = regulons)
df <- fread(infile, data.table = FALSE, check.names = FALSE)
stopifnot(ncol(df) >= 2)

rownames(df) <- df[[1]]
df <- df[, -1, drop = FALSE]

mat0 <- as.matrix(df)
storage.mode(mat0) <- "numeric"

# Optional: filter to top-variance regulons
mat <- mat0
if (!is.na(top_k_regulons_by_var) && top_k_regulons_by_var < ncol(mat)) {
  v <- apply(mat, 2, stats::var, na.rm = TRUE)
  keep <- names(sort(v, decreasing = TRUE))[seq_len(top_k_regulons_by_var)]
  mat <- mat[, keep, drop = FALSE]
}

# Normalize for visualization (no z-score scaling for RSS/JSD)
if (mode == "mss") {
  # Scale each regulon to [0,1] across cell types (0 = min, 1 = most specific cell type)
  m <- apply(mat, 2, max, na.rm = TRUE)
  m[m == 0 | is.na(m)] <- 1
  mat <- sweep(mat, 2, m, "/")
}

# Optional winsorization to stabilize color range
if (!is.null(cap_quantile) && is.finite(cap_quantile)) {
  cap <- stats::quantile(mat, probs = cap_quantile, na.rm = TRUE, names = FALSE)
  mat[mat > cap] <- cap
}

mat[is.na(mat)] <- 0

# Row annotation: assign lineage from cell type names
celltype <- rownames(mat)
lineage <- dplyr::case_when(
  grepl("^RG$|RG\\.div|Radial", celltype, ignore.case = TRUE)         ~ "RG",
  grepl("Div|Prog|IPC", celltype, ignore.case = TRUE)                 ~ "Progenitor/IPC",
  grepl("OPC", celltype, ignore.case = TRUE)                           ~ "OPC",
  grepl("Oligo", celltype, ignore.case = TRUE)                         ~ "Oligodendrocyte",
  grepl("Astro", celltype, ignore.case = TRUE)                         ~ "Astrocyte",
  grepl("Micro", celltype, ignore.case = TRUE)                         ~ "Microglia",
  grepl("Endo|Peri|Vasc", celltype, ignore.case = TRUE)                ~ "Vascular",
  grepl("Excit|Inhib|Neuron|Newborn", celltype, ignore.case = TRUE)    ~ "Neuronal",
  TRUE                                                                  ~ "Other"
)

anno_row <- data.frame(Lineage = factor(lineage))
rownames(anno_row) <- rownames(mat)

# Color palette (white -> green, for MSS in [0,1])
pal <- colorRampPalette(c("white", "#466812"))(256)
if (mode == "mss") {
  bk <- seq(0, max(mat, na.rm = TRUE), length.out = 257)
} else {
  bk <- seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = 257)
}

# Draw full heatmap
if (grepl("\\.pdf$", outfile, ignore.case = TRUE)) {
  pdf(outfile, width = 12, height = 7)
  on.exit(dev.off(), add = TRUE)
} else if (grepl("\\.png$", outfile, ignore.case = TRUE)) {
  png(outfile, width = 2400, height = 1400, res = 200)
  on.exit(dev.off(), add = TRUE)
}

pheatmap(
  mat,
  color  = pal,
  breaks = bk,
  border_color = NA,
  scale  = "none",
  cluster_rows = cluster_rows,
  cluster_cols = cluster_cols,
  clustering_distance_rows = dist_method,
  clustering_distance_cols = dist_method,
  clustering_method = clust_method,
  annotation_row  = anno_row,
  show_rownames   = TRUE,
  show_colnames   = (ncol(mat) <= show_colnames_if_leq),
  fontsize_row    = 9,
  fontsize_col    = 6
)

# ============================================================
# Part 2: Identify top-K cell-type-specific regulons
# ============================================================
top_k <- 30
mat_use    <- mat
cell_types <- rownames(mat_use)
regulons   <- colnames(mat_use)

# Compute delta matrix: delta[i,j] = score[i,j] - max(score[-i,j])
# Measures how much regulon j is more active in cell type i vs the next-best
delta <- matrix(NA_real_, nrow = nrow(mat_use), ncol = ncol(mat_use),
                dimnames = dimnames(mat_use))

for (i in seq_len(nrow(mat_use))) {
  others_max <- apply(mat_use[-i, , drop = FALSE], 2, max, na.rm = TRUE)
  delta[i, ] <- mat_use[i, ] - others_max
}

# Select top-K regulons per cell type by specificity delta
top_list <- lapply(cell_types, function(ct) {
  d   <- delta[ct, ]
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

fwrite(top_tab, file = "celltype_top30_regulons_by_delta.csv")

# ============================================================
# Part 3: Subset heatmap of top cell-type-specific regulons
# ============================================================
regs_union <- unique(top_tab$regulon)
mat_sub    <- mat_use[, regs_union, drop = FALSE]

# Order columns by winner cell type (most specific), then by max score
winner_ct <- apply(mat_sub, 2, function(x) rownames(mat_sub)[which.max(x)])
winner_ct <- factor(winner_ct, levels = cell_types)

ord_cols  <- order(winner_ct, apply(mat_sub, 2, max, na.rm = TRUE), decreasing = TRUE)
mat_sub   <- mat_sub[, ord_cols, drop = FALSE]
winner_ct <- winner_ct[ord_cols]

# Column annotation: which cell type "wins" each regulon
anno_col <- data.frame(WinnerCellType = winner_ct)
rownames(anno_col) <- colnames(mat_sub)

# Draw subset heatmap (columns ordered by winner, not clustered)
pal_sub <- colorRampPalette(c("white", "#466812"))(256)
bk_sub  <- seq(0, 1, length.out = 257)

pheatmap::pheatmap(
  mat_sub,
  color  = pal_sub,
  breaks = bk_sub,
  border_color = NA,
  scale  = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,   # Preserve winner-cell-type ordering
  annotation_row = anno_row,
  annotation_col = anno_col,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  fontsize_row   = 9,
  fontsize_col   = 3
)

# Save the subset matrix
write.csv(mat_sub, file = "celltype_regulon_specificity_top_subset.csv", quote = FALSE)
