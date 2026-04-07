#2025-10-15


# ====== packages ======
suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(stringr)
})

# ====== Read in data ======
df <- read_csv("100_summary_10_Blood_cell_traits.csv")   # or read_tsv("table.tsv")
# Required columns: Regulon, cell_type, Trait, Regions, Dataset,
#           p_RAS, p_RGS, p_TRS, q_RAS, q_RGS, q_TRS,
#           and z_TRS (if absent, TRS can be used instead)
stopifnot(all(c("Regulon","cell_type","Trait",
                "p_RAS","p_RGS","p_TRS","q_RAS","q_RGS","q_TRS") %in% names(df)))

# ====== Select the Dataset to plot (e.g., only 6.5K; remove this line to combine both datasets) ======
DATASET_TO_USE <- "6.5K"   # can change to "8.7K"
if ("Dataset" %in% names(df)) df <- dplyr::filter(df, Dataset == DATASET_TO_USE)

# ====== Value column: prefer z_TRS; if absent, use TRS and compute row-wise Z-score inline ======
value_col <- if ("z_TRS" %in% names(df)) "z_TRS" else if ("TRS" %in% names(df)) "TRS" else NA
stopifnot(!is.na(value_col))

# -- If there are duplicate (Regulon, cell_type, Trait) entries (e.g., multiple regions / RBPs),
#    take the mean for values; take the min for p/q (most conservative significance)
agg <- df %>%
  group_by(Regulon, cell_type, Trait) %>%
  summarise(
    value  = mean(.data[[value_col]], na.rm = TRUE),
    p_RAS  = min(p_RAS, na.rm = TRUE),
    p_RGS  = min(p_RGS, na.rm = TRUE),
    p_TRS  = min(p_TRS, na.rm = TRUE),
    q_RAS  = min(q_RAS, na.rm = TRUE),
    q_RGS  = min(q_RGS, na.rm = TRUE),
    q_TRS  = min(q_TRS, na.rm = TRUE),
    .groups = "drop"
  )

# ====== Assign significance labels: priority ** > * > # ======
sign_table <- agg %>%
  rowwise() %>%
  mutate(
    star = case_when(
      # 1) All three q < 0.05 -> ** (highest priority)
      q_RAS < 0.05 & q_RGS < 0.05 & q_TRS < 0.05 ~ "**",
      # 2) All three q <= 0.1 but not all < 0.05 (i.e., not the condition above) -> *
      max(q_RAS, q_RGS, q_TRS, na.rm = TRUE) <= 0.1 &
        !(q_RAS < 0.05 & q_RGS < 0.05 & q_TRS < 0.05) ~ "*",
      # 3) Does not meet the above, but all three p < 0.05 -> #
      p_RAS < 0.05 & p_RGS < 0.05 & p_TRS < 0.05 ~ "#",
      TRUE ~ ""
    )
  ) %>%
  ungroup()

# ====== Wide table: columns = cell_type | Trait ======
wide <- sign_table %>%
  mutate(colkey = paste0(cell_type, " | ", Trait)) %>%
  select(Regulon, colkey, value, star)

# Value matrix
mat <- wide %>%
  select(Regulon, colkey, value) %>%
  pivot_wider(names_from = colkey, values_from = value) %>%
  arrange(Regulon) %>%
  column_to_rownames("Regulon") %>%
  as.matrix()

# Label matrix
lab <- wide %>%
  select(Regulon, colkey, star) %>%
  pivot_wider(names_from = colkey, values_from = star) %>%
  arrange(Regulon) %>%
  column_to_rownames("Regulon") %>%
  as.matrix()

# If z_TRS is not available and TRS is used, apply row-wise Z-score to avoid scale differences
if (value_col != "z_TRS") {
  mat <- t(scale(t(mat)))
}

# ====== Column order and faceting: split panels by cell_type, sort by Trait within each panel ======
# ====== Safe row-wise Z-score function (sd=0 -> all zeros; preserves NA) ======
row_zsafe <- function(m) {
  t(apply(m, 1, function(v) {
    mu  <- mean(v, na.rm = TRUE)
    sdv <- sd(v,   na.rm = TRUE)
    if (!is.finite(sdv) || sdv == 0) return(rep(0, length(v)))
    (v - mu) / sdv
  }))
}

col_meta <- tibble(col = colnames(mat)) %>%
  tidyr::separate(col, into = c("cell_type","Trait"),
                  sep = " \\| ", remove = FALSE)

celltype_levels <- c("Monocytes","T_cells")

col_meta <- col_meta %>%
  mutate(cell_type = factor(cell_type, levels = celltype_levels))


# ====== Clean up: remove rows/columns that are entirely NA ======
keep_rows <- rowSums(is.finite(mat)) > 0
keep_cols <- colSums(is.finite(mat)) > 0
mat <- mat[keep_rows, keep_cols, drop = FALSE]
lab <- lab[keep_rows, keep_cols, drop = FALSE]
col_meta <- col_meta[keep_cols, , drop = FALSE]

# ====== (If using TRS instead of z_TRS) Safe row-wise standardization ======
if (value_col != "z_TRS") {
  mat <- row_zsafe(mat)
}

# ====== Recalculate color scale (avoid ±Inf) =====
rng <- range(mat[is.finite(mat)], na.rm = TRUE)
pal <- colorRamp2(c(rng[1], 0, rng[2]),
                  c("#2166AC", "#F7F7F7", "#B2182B"))

# ====== Column order & faceting ======
# (If you already have col_meta/ord, you can use them directly here; otherwise same as before)
ord <- order(col_meta$cell_type, col_meta$Trait)
mat <- mat[, ord, drop = FALSE]
lab <- lab[, ord, drop = FALSE]
col_meta <- col_meta[ord, , drop = FALSE]

# Top annotation (if top_anno was not defined earlier)
celltype_cols <- c(
  "Monocytes"   = "#E6C229",
  "T_cells"      = "#1B9E77"
)
top_anno <- ComplexHeatmap::HeatmapAnnotation(
  `Cell type` = col_meta$cell_type,
  col = list(`Cell type` = celltype_cols),
  annotation_name_side = "left"
)

# ====== Custom clustering: temporarily replace NA with 0 for hclust to avoid dist() errors ======
cluster_na_ok <- function(m) {
  m2 <- m
  m2[!is.finite(m2)] <- 0
  hclust(dist(m2))
}

# ====== Plot (with significance labels inside cells) =====
ht <- Heatmap(
  mat,
  name = if (value_col == "z_TRS") "z_TRS" else "TRS (row-z)",
  col = pal, na_col = "grey90",
  # Use custom clustering function to avoid NA errors; set FALSE if clustering is not desired
  cluster_rows = cluster_na_ok,
  cluster_columns = T,
  column_split = col_meta$cell_type,
  column_gap   = unit(2.5, "mm"),
  top_annotation = top_anno,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 6),
  rect_gp = gpar(col = "white"),
  heatmap_legend_param = list(title = "TRS", at = seq(-2, 2, 1)),
  cell_fun = function(j, i, x, y, w, h, fill) {
    s <- lab[i, j]
    if (!is.na(s) && nzchar(s)) {
      val <- mat[i, j]
      # Automatically choose text color based on background brightness
      bg  <- pal(ifelse(is.finite(val), val, 0))
      rgb <- col2rgb(bg)
      lum <- (0.2126*rgb[1] + 0.7152*rgb[2] + 0.0722*rgb[3]) / 255
      txt_col <- ifelse(lum < 0.5, "white", "black")
      grid.text(s, x, y, gp = gpar(col = txt_col, fontsize = 5, fontface = "bold"))
    }
  }
)

draw(ht)


##---- Final version
# 1) Use a symmetric range (98th percentile of actual distribution to avoid extreme values darkening colors)
L <- quantile(abs(mat[is.finite(mat)]), 0.98, na.rm = TRUE)
L <- max(L, 1.5)  # lower bound to avoid too narrow a range

# 2) Use brighter red/blue at endpoints with pure white in the middle
pal <- circlize::colorRamp2(
  c(-L, 0, L),
  c("#5DADE2", "#FFFFFF", "#E74C3C")  # bright blue - white - bright red
)

# install.packages("scico")
library(scico)
L <- quantile(abs(mat[is.finite(mat)]), 0.98, na.rm = TRUE)
pal <- circlize::colorRamp2(
  seq(-L, L, length.out = 11),
  scico(11, palette = "vik")  # bright and perceptually uniform
)

ht <- Heatmap(
  mat,
  name = if (value_col == "z_TRS") "z_TRS" else "TRS (row-z)",
  col = pal, na_col = "grey90",
  cluster_rows    = cluster_na_ok,   # set FALSE if clustering is not desired
  cluster_columns = cluster_na_ok,   # columns may also contain NA; use safe clustering
  column_split = col_meta$cell_type,
  column_gap   = unit(2.5, "mm"),
  top_annotation = top_anno,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 6),
  rect_gp = gpar(col = "white"),
  heatmap_legend_param = list(title = "TRS", at = seq(-2, 2, 1)),
  cell_fun = function(j, i, x, y, w, h, fill) {
    s <- lab[i, j]
    if (!is.na(s) && nzchar(s)) {
      grid.text(s, x, y, gp = gpar(col = "white", fontsize = 6, fontface = "bold"))
    }
  }
)

draw(ht)
