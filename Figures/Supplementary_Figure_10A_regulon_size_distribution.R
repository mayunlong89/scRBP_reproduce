#2025-10-13
setwd("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/02_scRBP_trs_10blood_cell_traits_100times_v2/00_6.5K_regulons_gmt")
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
# ===== 1) Tidy input (reusing existing result) =====
# If your result already has columns named Regulon / Regulon_size, you can skip this step
result <- read.csv("summary_all_4regions_regulon_gene_counts_min1genes.csv")
regulon_data <- result %>%
  rename(Regulon = Regulon, Regulon_size = TargetGeneCount)
# Optional: if +1 was added previously, you can remove it here (keep your original logic if not needed)
# regulon_data <- regulon_data %>% mutate(Regulon_size = pmax(Regulon_size - 1, 0))
# ===== 2) Plot parameters =====
threshold <- 10
topN_to_label <- 12               # Only label the top N largest by size (to avoid overcrowding)
use_log10_x  <- FALSE             # Set to TRUE if log10 x-axis is needed
# Count of those >= threshold (for annotation on the plot)
n_ge  <- sum(regulon_data$Regulon_size >= threshold)
pct_ge <- mean(regulon_data$Regulon_size >= threshold)
# Select points to label (top N among those above the threshold)
lab_dat <- regulon_data %>%
  filter(Regulon_size >= threshold) %>%
  slice_max(Regulon_size, n = topN_to_label, with_ties = FALSE)
# ===== 3) Plot (histogram + density) =====
p <- ggplot(regulon_data, aes(x = Regulon_size)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.35) +
  geom_density(linewidth = 1.0) +
  geom_rug(sides = "b", alpha = 0.25) +
  geom_vline(xintercept = threshold, color = "red", linetype = "dashed") +
  annotate(
    "label",
    x = threshold, y = Inf, vjust = 1.2, hjust = 0,
    label = paste0("≥ ", threshold, ": ", n_ge, " (", percent(pct_ge), ")"),
    color = "red"
  ) +
  labs(
    title = "Regulon size distribution",
    x = "Regulon size (# targets)",
    y = "Density"
  ) +
  theme_minimal(base_size = 12)
# Mark regulons >= threshold near the x-axis, and label top N with ggrepel
p <- p +
  geom_point(
    data = regulon_data %>% filter(Regulon_size >= threshold),
    aes(x = Regulon_size, y = 0),
    inherit.aes = FALSE, size = 1.6, color = "red", alpha = 0.7
  ) +
  ggrepel::geom_text_repel(
    data = lab_dat,
    aes(x = Regulon_size, y = 0, label = Regulon),
    inherit.aes = FALSE,
    nudge_y = 0.03,    # Nudge labels slightly upward
    min.segment.length = 0,
    max.overlaps = Inf,
    size = 3.2
  )
# Optional: log10 x-axis (clearer for long-tail distributions)
if (use_log10_x) {
  p <- p + scale_x_continuous(trans = "log10", breaks = c(1, 3, 10, 30, 100, 300))
}
print(p)
# Save if needed
ggsave("regulon_size_density.png", p, width = 8.5, height = 5.2, dpi = 300)

##-------- Calculate the median, max, and min
# Retain regulons with 10 or more genes
regulon_data2 <- regulon_data[which(regulon_data$Regulon_size>9),]
summary(regulon_data2)
length(regulon_data2$Regulon_size)
cols_fixed <- c(
  "RBC"            = "#D73027",
  "Platelets"      = "#F1A340",
  "Cycling"        = "#7F7F7F",
  "Baso/Eos"       = "#F0E442",
  "HPC"            = "#1B9E77",
  "Plasma"         = "#2C7FB8",
  "B"              = "#1A9850",
  "DC"             = "#1F9A8A",
  "ILC"            = "#6CC7BD",
  "NK"             = "#225EA8",
  "MAIT"           = "#B12A90",
  "T reg"          = "#7CAE00",
  "T g/d"          = "#F564E3",
  "T CD8+"         = "#6A51A3",
  "T CD4+"         = "#FB6A4A",
  "Monocyte CD14"  = "#E69F00",  # orange-yellow
  "Monocyte CD16"  = "#8C564B"   # brown
)
cols_fixed <- c(
  "Cycling"        = "#D6A18C",
  "B"              = "#E69F00",
  "Platelets"      = "#E6D86B",
  "mDC"            = "#C3D5DE",
  "pDC"           = "#A02D28",
  "ILC"            = "#6CC7BD",
  "NK"             = "#225EA8",
  "MAIT"           = "#B12A90",
  "T reg"          = "#6B7FA8",
  "T g/d"          = "#F564E3",
  "T CD8+"         = "#6A51A3",
  "T CD4+"         = "#FB6A4A",
  "Monocyte CD14"  = "#4A774D",  # orange-yellow
  "Monocyte CD16"  = "#8DC890"   # brown
)
# Named vector with fixed order (recommended to copy directly)
my_colors <- c(
  "#9E3B30","#D6A18C","#C07E6B","#A63D36","#B45D4A","#8E2C28","#C5794B","#A02D28",
  "#D96CA0","#9429A4","#7F7F99","#6BAF4E","#C2B78C","#3F6843","#4A774D","#588B55",
  "#6AAE6B","#8DC890","#31693E","#E6D86B","#F3D447","#E2A957","#6699CC","#286BB0",
  "#1F4489","#6B7FA8","#C3D5DE"
)
names(my_colors) <- levels(obj$annotation_broad2)[1:length(my_colors)]
p <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "annotation_broad2",
  cols = my_colors,
  label = TRUE, repel = TRUE, label.size = 4, pt.size = 0.35, shuffle = TRUE
) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 16))
ggsave("PBMC_umap_annotation_broad2_custom_palette.png",
       p, device = ragg::agg_png, width = 8, height = 8, units = "in", dpi = 600, bg = "white")
