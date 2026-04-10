


#Step 2.2: cell-type-by-RBP mean expression--heatmap-plotting
#Step 2.2: cell-type-by-RBP mean expression--heatmap-plotting

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# =========================
# User knobs
# =========================
expr_file <- "celltype_by_RBP_mean_expression.csv"
rss_file  <- "10X_fetal_brain_590Kcells_ras_ct_all4regions.symbol_main_celltype_rss.csv"

out_main_pdf <- "Fig_main_RBPlevel_exprZ_plus_RSSdot.pdf"
out_supp_pdf <- "Fig_supp_RegulonLevel_exprZ_plus_RSSdot.pdf"

# Expression z-score:
#   "rbp" = for each (RBP gene), z-score across cell types  (recommended)
zscore_by <- "rbp"

# Main fig: collapse multiple regulons per RBP
#   "max" (recommended) / "mean" / "median"
collapse_fun <- "max"

# Select top features to plot (avoid overcrowding)
top_n_rbps_main     <- 60   # main figure RBPs
top_n_regulons_supp <- 60  # supplementary regulons

select_mode_main <- "maxRSS"     # "maxRSS" or "varRSS" based on collapsed RBP-level RSS
select_mode_supp <- "maxRSS"     # "maxRSS" or "varRSS" based on regulon-level RSS

# clustering orders
cluster_celltypes <- TRUE
cluster_features  <- TRUE
dist_method  <- "euclidean"
clust_method <- "complete"

# aesthetics
dot_size_range <- c(0, 5)
tile_limits <- c(-2, 2)   # z-score color range (clip). set NULL to disable clipping

# Regexp to map regulon -> base RBP gene
# Works for: AGO2_Intron, AGO2_Introns, AGO2_CDS, AGO2_3'UTR, AGO2_5'UTR, etc.
region_pattern <- "(?:_)?(Intron|Introns|CDS|3'?UTR|5'?UTR)$"
# =========================

# -------------------------
# Helpers
# -------------------------
zscore_vec <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

collapse_vec <- function(x, fun = "max") {
  fun <- match.arg(fun, c("max", "mean", "median"))
  if (fun == "max")    return(max(x, na.rm = TRUE))
  if (fun == "mean")   return(mean(x, na.rm = TRUE))
  if (fun == "median") return(stats::median(x, na.rm = TRUE))
}

# regulon -> base gene
regulon_to_gene <- function(reg_names) {
  # strip trailing region label
  g <- str_replace(reg_names, region_pattern, "")
  # also strip an extra underscore if left behind
  g <- str_replace(g, "_$", "")
  g
}

# clustering order from a numeric matrix
cluster_order <- function(mat, margin = 1) {
  # margin=1 order rows, margin=2 order cols
  if (margin == 1) {
    d <- dist(mat, method = dist_method)
    hc <- hclust(d, method = clust_method)
    return(hc$labels[hc$order])
  } else {
    d <- dist(t(mat), method = dist_method)
    hc <- hclust(d, method = clust_method)
    return(hc$labels[hc$order])
  }
}

# -------------------------
# Read matrices
# -------------------------
expr_df <- fread(expr_file, data.table = FALSE, check.names = FALSE)
rss_df  <- fread(rss_file,  data.table = FALSE, check.names = FALSE)

stopifnot(ncol(expr_df) >= 2, ncol(rss_df) >= 2)

rownames(expr_df) <- expr_df[[1]]; expr_df <- expr_df[, -1, drop = FALSE]
rownames(rss_df)  <- rss_df[[1]];  rss_df  <- rss_df[,  -1, drop = FALSE]

expr_mat0 <- as.matrix(expr_df); storage.mode(expr_mat0) <- "numeric"
rss_mat0  <- as.matrix(rss_df);  storage.mode(rss_mat0)  <- "numeric"

# align cell types first
common_ct <- intersect(rownames(expr_mat0), rownames(rss_mat0))
if (length(common_ct) < 2) stop("Too few overlapping cell types between expr and RSS.")

expr_mat0 <- expr_mat0[common_ct, , drop = FALSE]
rss_mat0  <- rss_mat0[ common_ct, , drop = FALSE]

# -------------------------
# Parse mapping: regulon -> gene
# -------------------------
regulons <- colnames(rss_mat0)
gene_of_reg <- regulon_to_gene(regulons)

map_dt <- data.frame(
  regulon = regulons,
  gene    = gene_of_reg,
  stringsAsFactors = FALSE
)

# intersect genes with expression columns
common_gene <- intersect(colnames(expr_mat0), unique(map_dt$gene))
if (length(common_gene) < 2) {
  # show examples to help debug
  msg <- paste0(
    "No/too few overlaps after mapping regulon->gene.\n",
    "Example expr genes: ", paste(head(colnames(expr_mat0), 10), collapse=", "), "\n",
    "Example regulons:   ", paste(head(regulons, 10), collapse=", "), "\n",
    "Example mapped gene:", paste(head(gene_of_reg, 10), collapse=", "), "\n"
  )
  stop(msg)
}

# Keep only regulons whose gene exists in expr
keep_reg <- map_dt$regulon[map_dt$gene %in% common_gene]
rss_mat  <- rss_mat0[, keep_reg, drop = FALSE]
map_dt   <- map_dt[map_dt$regulon %in% keep_reg, , drop = FALSE]

# Keep only expr genes present
expr_mat <- expr_mat0[, common_gene, drop = FALSE]

# -------------------------
# Compute expression z-score
# -------------------------
expr_z_gene <- expr_mat

if (zscore_by == "rbp") {
  expr_z_gene <- apply(expr_mat, 2, zscore_vec)
  expr_z_gene <- as.matrix(expr_z_gene)
  rownames(expr_z_gene) <- rownames(expr_mat)
} else {
  stop("Only zscore_by='rbp' implemented here (recommended).")
}

# optional clip
if (!is.null(tile_limits) && length(tile_limits) == 2) {
  expr_z_gene[expr_z_gene < tile_limits[1]] <- tile_limits[1]
  expr_z_gene[expr_z_gene > tile_limits[2]] <- tile_limits[2]
}

# ============================================================
# (A) MAIN FIG: collapse regulon RSS -> RBP-level RSS
# ============================================================
# Build celltype × gene RSS by collapsing all regulons per gene
genes_main <- common_gene

rss_gene <- sapply(genes_main, function(g) {
  regs <- map_dt$regulon[map_dt$gene == g]
  if (length(regs) == 0) return(rep(NA_real_, nrow(rss_mat)))
  x <- rss_mat[, regs, drop = FALSE]
  apply(x, 1, collapse_vec, fun = collapse_fun)
})
rss_gene <- as.matrix(rss_gene)
rownames(rss_gene) <- rownames(rss_mat)

# Select top RBPs for main
score_main <- if (select_mode_main == "maxRSS") {
  apply(rss_gene, 2, max, na.rm = TRUE)
} else {
  apply(rss_gene, 2, var, na.rm = TRUE)
}
keep_gene_main <- names(sort(score_main, decreasing = TRUE))[seq_len(min(top_n_rbps_main, length(score_main)))]

rss_gene_main   <- rss_gene[, keep_gene_main, drop = FALSE]
expr_z_main     <- expr_z_gene[, keep_gene_main, drop = FALSE]

# cluster order
ct_order_main <- rownames(expr_z_main)
gene_order_main <- colnames(expr_z_main)

if (cluster_celltypes) ct_order_main <- cluster_order(expr_z_main, margin = 1)
if (cluster_features)  gene_order_main <- cluster_order(expr_z_main, margin = 2)

expr_z_main <- expr_z_main[ct_order_main, gene_order_main, drop = FALSE]
rss_gene_main <- rss_gene_main[ct_order_main, gene_order_main, drop = FALSE]

# long df
df_main <- as.data.frame(as.table(expr_z_main), stringsAsFactors = FALSE) %>%
  rename(cell_type = Var1, feature = Var2, expr_z = Freq) %>%
  left_join(
    as.data.frame(as.table(rss_gene_main), stringsAsFactors = FALSE) %>%
      rename(cell_type = Var1, feature = Var2, RSS = Freq),
    by = c("cell_type", "feature")
  )

df_main$cell_type <- factor(df_main$cell_type, levels = ct_order_main)
df_main$feature   <- factor(df_main$feature,   levels = gene_order_main)

p_main <- ggplot(df_main, aes(x =feature, y = cell_type)) +
  geom_tile(aes(fill = expr_z)) +
  geom_point(aes(size = RSS), color = "black", alpha = 0.85) +
  scale_fill_gradient2(
    low = "#2b8cbe", mid = "white", high = "#f03b20", midpoint = 0,
    name = "z-score\nRBP expr"
  ) +
  scale_size_continuous(range = dot_size_range, name = paste0("RSS (", collapse_fun, ")")) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9)
  )

print(p_main)

ggsave(out_main_pdf, p_main, width = 7.5, height = 10, useDingbats = FALSE)
message("Saved main: ", out_main_pdf)

# ============================================================
# (B) SUPP FIG: keep all regulons (RBP_region)
# ============================================================
# For each regulon column, use expression z-score of its gene as the tile fill
# Create expr_z_regulon matrix with same dim as rss_mat
# rows = cell types, cols = regulons
expr_z_reg <- matrix(NA_real_, nrow = nrow(rss_mat), ncol = ncol(rss_mat),
                     dimnames = list(rownames(rss_mat), colnames(rss_mat)))

for (j in seq_len(ncol(rss_mat))) {
  reg <- colnames(rss_mat)[j]
  g   <- map_dt$gene[match(reg, map_dt$regulon)]
  expr_z_reg[, j] <- expr_z_gene[, g]
}

# Select top regulons for supplementary
score_supp <- if (select_mode_supp == "maxRSS") {
  apply(rss_mat, 2, max, na.rm = TRUE)
} else {
  apply(rss_mat, 2, var, na.rm = TRUE)
}
keep_reg_supp <- names(sort(score_supp, decreasing = TRUE))[seq_len(min(top_n_regulons_supp, length(score_supp)))]

rss_supp   <- rss_mat[, keep_reg_supp, drop = FALSE]
exprz_supp <- expr_z_reg[, keep_reg_supp, drop = FALSE]

# cluster order (can use exprz or rss; here use exprz)
ct_order_supp <- rownames(exprz_supp)
reg_order_supp <- colnames(exprz_supp)

if (cluster_celltypes) ct_order_supp <- cluster_order(exprz_supp, margin = 1)
if (cluster_features)  reg_order_supp <- cluster_order(exprz_supp, margin = 2)

exprz_supp <- exprz_supp[ct_order_supp, reg_order_supp, drop = FALSE]
rss_supp   <- rss_supp[ct_order_supp,  reg_order_supp, drop = FALSE]

df_supp <- as.data.frame(as.table(exprz_supp), stringsAsFactors = FALSE) %>%
  rename(cell_type = Var1, feature = Var2, expr_z = Freq) %>%
  left_join(
    as.data.frame(as.table(rss_supp), stringsAsFactors = FALSE) %>%
      rename(cell_type = Var1, feature = Var2, RSS = Freq),
    by = c("cell_type", "feature")
  )

df_supp$cell_type <- factor(df_supp$cell_type, levels = ct_order_supp)
df_supp$feature   <- factor(df_supp$feature,   levels = reg_order_supp)

p_supp <- ggplot(df_supp, aes(x = feature, y =cell_type )) +
  geom_tile(aes(fill = expr_z)) +
  geom_point(aes(size = RSS), color = "black", alpha = 0.85) +
  scale_fill_gradient2(
    low = "#2b8cbe", mid = "white", high = "#f03b20", midpoint = 0,
    name = "z-score\nRBP expr"
  ) +
  scale_size_continuous(range = dot_size_range, name = "RSS") +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 6)
  )


print(p_supp)
ggsave(out_supp_pdf, p_supp, width = 8.5, height = 12, useDingbats = FALSE)
message("Saved supp: ", out_supp_pdf)



