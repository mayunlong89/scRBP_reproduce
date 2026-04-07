# ====== Required packages ======
suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
})

# ====== Read in data ======
# If TSV, use read_tsv; if CSV, use read_csv
# df <- read_tsv("trs_table.tsv")   # ← TSV
# df <- read_csv("trs_table.csv")   # ← CSV
# You can also directly assign an existing data.frame to df
df <- read.csv("00_4celltype_heatmap.csv")   # Change to your file name as needed

# Check for required columns
need_cols <- c("Regulon","cell_type","TRS_6.5K","TRS_8.7K")
stopifnot(all(need_cols %in% colnames(df)))

# ====== Reshape into "cell_type × (6.5K, 8.7K)" columns ======
long <- df %>%
  select(Regulon, cell_type, TRS_6.5K, TRS_8.7K) %>%
  pivot_longer(cols = starts_with("TRS_"),
               names_to = "Dataset", values_to = "TRS") %>%
  mutate(Dataset = recode(Dataset, `TRS_6.5K`="6.5K", `TRS_8.7K`="8.7K"))

# Create composite column names as "cell_type | Dataset"
wide <- long %>%
  mutate(colkey = paste0(cell_type, " | ", Dataset)) %>%
  select(Regulon, colkey, TRS) %>%
  # If duplicate (Regulon, colkey) entries exist, take the mean
  group_by(Regulon, colkey) %>%
  summarise(TRS = mean(TRS, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = colkey, values_from = TRS)

# Convert to matrix
mat <- wide %>% column_to_rownames("Regulon") %>% as.matrix()

# ====== Column order: group by cell_type, with 6.5K before 8.7K within each group ======
col_meta <- tibble(col = colnames(mat)) %>%
  separate(col, into = c("cell_type","Dataset"), sep = " \\| ", remove = FALSE) %>%
  mutate(Dataset = factor(Dataset, levels = c("6.5K","8.7K")))
ord <- order(col_meta$cell_type, col_meta$Dataset)
mat <- mat[, ord, drop = FALSE]
col_meta <- col_meta[ord, , drop = FALSE]

# ====== (Optional) Row-wise Z-score normalization and value capping ======
row_z  <- TRUE          # Set to FALSE to skip Z-score
cap_at <- 3             # Cap Z-scores to [-3, 3]; set to Inf to disable capping

if (row_z) {
  mat <- t(scale(t(mat)))            # Row-wise Z-score
  if (is.finite(cap_at)) {
    mat[mat >  cap_at] <-  cap_at
    mat[mat < -cap_at] <- -cap_at
  }
}

# ====== Colors and annotations ======
# Continuous color palette (negative - midpoint - positive)
pal <- colorRamp2(c(min(mat, na.rm=TRUE), 0, max(mat, na.rm=TRUE)),
                  c("#2166AC", "#F7F7F7", "#B2182B"))

top_anno <- HeatmapAnnotation(
  Dataset = col_meta$Dataset,
  col = list(Dataset = c("6.5K" = "#1f77b4", "8.7K" = "#ff7f0e"))
)

# ====== Plot ======
ht <- Heatmap(
  mat,
  name = if (row_z) "TRS (row-z)" else "TRS",
  col = pal,
  na_col = "grey90",
  top_annotation = top_anno,
  column_split = col_meta$cell_type,       # Facet by cell type
  column_title = "Cell type | Dataset",
  row_title = "Regulon",
  cluster_rows = TRUE,
  cluster_columns = FALSE,                  # Display columns in the specified order
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9)
)

# Display on screen
draw(ht)

# ====== Save to file (optional) ======
# PNG
png("regulon_TRS_heatmap.png", width = 2400, height = 1500, res = 300)
draw(ht)
dev.off()

# PDF (vector graphics, suitable for publication)
pdf("regulon_TRS_heatmap.pdf", width = 14, height = 9)
draw(ht)
dev.off()



# ====== Column order: first group (Monocytes vs T cells), then sort by cell type and dataset within each group ======
library(stringr)

col_meta <- tibble(col = colnames(mat)) %>%
  separate(col, into = c("cell_type","Dataset"), sep = " \\| ", remove = FALSE) %>%
  mutate(
    Dataset  = factor(Dataset, levels = c("6.5K","8.7K")),
    # Define two major groups: Monocytes and T cells
    Group    = factor(case_when(
      str_detect(cell_type, "Monocytes") ~ "Monocytes",
      str_detect(cell_type, "T cells")   ~ "T cells",
      TRUE ~ "Other"
    ), levels = c("Monocytes","T cells","Other")),
    # Cell type order within each major group (adjust as needed)
    CellType = factor(cell_type,
                      levels = c("CD14+ Monocytes", "FCGR3A+ Monocytes",
                                 "CD4+ T cells", "CD8+ T cells"))
  )

ord <- order(col_meta$Group, col_meta$CellType, col_meta$Dataset)
mat <- mat[, ord, drop = FALSE]
col_meta <- col_meta[ord, , drop = FALSE]

# ====== Colors and annotations (unchanged, or continue using your own as needed) =====
pal <- colorRamp2(c(min(mat, na.rm=TRUE), 0, max(mat, na.rm=TRUE)),
                  c("#2166AC", "#F7F7F7", "#B2182B"))

celltype_cols <- c(
  "CD14+ Monocytes"   = "#E6C229",
  "FCGR3A+ Monocytes" = "#B88A0E",
  "CD4+ T cells"      = "#4DAF4A",
  "CD8+ T cells"      = "#1B9E77"
)

top_anno <- HeatmapAnnotation(
  Dataset  = col_meta$Dataset,
  `Cell type` = col_meta$cell_type,
  col = list(
    Dataset  = c("6.5K" = "#1f77b4", "8.7K" = "#ff7f0e"),
    `Cell type` = celltype_cols
  )
)

# ====== Plot: split by Group (Monocytes / T cells), with whitespace between blocks ======
ht <- Heatmap(
  mat,
  name = if (row_z) "TRS (row-z)" else "TRS",
  col = pal, na_col = "grey90",
  top_annotation   = top_anno,
  column_split     = col_meta$Group,            # Split into Monocytes and T cells
  column_gap       = unit(2.5, "mm"),           # Gap between the two blocks
  column_title     = "Cell type | Dataset",
  row_title        = "Regulon",
  cluster_rows     = TRUE,
  cluster_columns  = FALSE,                     # Use custom order
  row_names_gp     = gpar(fontsize = 8),
  column_names_gp  = gpar(fontsize = 9)
)

draw(ht)
