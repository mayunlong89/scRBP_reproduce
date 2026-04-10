# --- 2025-12-12
# Spearman correlation heatmap (lower triangle) across 8 autoimmune disease TRS scores

# ============================================================
# Load required packages
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
})

# ============================================================
# Step 1: Read TRS data
# ============================================================
infile <- "01_heatmap_allTRS_only.csv"
dat <- readr::read_csv(infile, show_col_types = FALSE)

# ============================================================
# Step 2: Extract 8 disease TRS columns and aggregate by ID
# ============================================================
# Disease order (determines row/column order in the correlation heatmap)
trait_order <- c("IBD_TRS", "UC_TRS", "CD_TRS",
                 "MS_TRS", "RA_TRS", "T1D_TRS",
                 "PBC_TRS", "SLE_TRS")

trs_cols <- trait_order
stopifnot(all(trs_cols %in% names(dat)))

# Average TRS values per ID (in case of duplicate IDs)
mat_by_id <- dat |>
  dplyr::select(ID, dplyr::all_of(trs_cols)) |>
  dplyr::group_by(ID) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(trs_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Numeric matrix: rows = IDs, columns = 8 diseases
mat <- mat_by_id |>
  dplyr::select(-ID) |>
  as.data.frame()

# ============================================================
# Step 3: Compute pairwise correlation matrix (Spearman)
# ============================================================
cor_mat <- cor(
  mat,
  use    = "pairwise.complete.obs",
  method = "spearman"   # change to "pearson" if needed
)

# Enforce the desired row/column order
cor_mat <- cor_mat[trait_order, trait_order]

# ============================================================
# Step 4: Reshape into lower-triangle data frame for ggplot
# ============================================================
cor_df <- reshape2::melt(
  cor_mat,
  varnames   = c("Var1", "Var2"),
  value.name = "rho"
)

cor_df <- cor_df |>
  dplyr::mutate(
    Var1 = factor(Var1, levels = trait_order),
    Var2 = factor(Var2, levels = trait_order)
  )

# Keep lower triangle only (Var1 index > Var2 index)
cor_df_lower <- cor_df |>
  dplyr::filter(as.numeric(Var1) > as.numeric(Var2))

# ============================================================
# Step 5: Plot lower-triangle correlation heatmap with values
# ============================================================
p <- ggplot(cor_df_lower, aes(x = Var1, y = Var2, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0.3,
    limits   = c(0.3, 0.9),
    name     = "Spearman \u03c1"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "red"),
    axis.text.y = element_text(color = "red"),
    axis.title  = element_blank(),
    panel.grid  = element_blank()
  )

print(p)

# ============================================================
# Cell type color palette (for use in related figures)
# ============================================================
cols_fixed <- c(
  "Cycling"        = "#D6A18C",
  "B"              = "#E69F00",
  "Platelets"      = "#E6D86B",
  "mDC"            = "#C3D5DE",
  "pDC"            = "#A02D28",
  "ILC"            = "#6CC7BD",
  "NK"             = "#225EA8",
  "MAIT"           = "#B12A90",
  "T reg"          = "#6B7FA8",
  "T g/d"          = "#F564E3",
  "T CD8+"         = "#6A51A3",
  "T CD4+"         = "#FB6A4A",
  "Monocyte CD14"  = "#4A774D",
  "Monocyte CD16"  = "#8DC890"
)
