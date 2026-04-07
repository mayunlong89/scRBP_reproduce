##_-----correlation---for 10 blood cell traits
##_-----correlation---for 10 blood cell traits
##_-----correlation---for 10 blood cell traits
# ====== packages ======
suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
#!/usr/bin/env Rscript
# Read the "wide table" (containing ID, Regulon, cell_type, RBP, Regions, and multiple TRS_xxx columns),
# compute pairwise correlations among the 10 phenotypes (TRS_*),
# and draw a lower-triangle heatmap with values (similar to the example figure).
suppressPackageStartupMessages({
  library(tidyverse)
  library(corrplot)
})
# ========= 1) Read in data =========
# Change to your file name
infile <- "100_summary_10_Blood_cell_traits_correlation.csv"
dat <- readr::read_csv(infile, show_col_types = FALSE)
# ========= 2) Extract TRS_* columns, align by ID (if duplicate IDs exist, take the mean) =========
trs_cols <- grep("^TRS_", names(dat), value = TRUE)
stopifnot(length(trs_cols) >= 2)
# One row per ID (take mean in case of duplicate IDs)
mat_by_id <- dat |>
  select(ID, all_of(trs_cols)) |>
  group_by(ID) |>
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
# Convert to numeric matrix (rows = IDs, columns = 10 phenotypes)
mat <- mat_by_id |> select(-ID) |> as.data.frame()
# ========= 3) Correlation matrix (Spearman, more robust; can also use "pearson") =========
cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "spearman")
# Also create a ggplot version (optional, for publication-style figure)
library(reshape2)
melted <- reshape2::melt(cor_mat, varnames = c("Var1","Var2"), value.name = "rho") |>
  mutate(keep = as.numeric(Var1) > as.numeric(Var2)) |>
  filter(keep)
p <- ggplot(melted, aes(Var1, Var2, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0.5, limits = c(0.49,0.8), name = "Spearman ρ") +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "red"),
        axis.text.y = element_text(color = "red"),
        axis.title = element_blank(),
        panel.grid = element_blank())
p
ggsave("trs_correlation_heatmap_ggplot.png", p, width = 12, height = 9, dpi = 240)
# ======== Optional: export correlation matrix to CSV ========
readr::write_csv(as.data.frame(cor_mat) |>
                   mutate(Phenotype = rownames(cor_mat)) |>
                   relocate(Phenotype),
                 "trs_correlation_matrix.csv")
