# --- 2026-01-07
# Regulon size distribution and PBMC UMAP with custom color palettes

# ============================================================
# Load required packages
# ============================================================
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)

# ============================================================
# Part 1: Regulon size distribution (fetal brain, 4 mRNA regions)
# ============================================================

# Read regulon gene count summary
result <- read.csv("summary_all_4regions_regulon_gene_counts_min1genes.csv")

regulon_data <- result %>%
  rename(Regulon = Regulon, Regulon_size = TargetGeneCount)

# Plot parameters
threshold     <- 10
topN_to_label <- 12       # Label only the top N largest regulons to avoid clutter
use_log10_x   <- FALSE    # Set TRUE for log10 x-axis (useful for long-tailed distributions)

# Count regulons at or above the threshold
n_ge   <- sum(regulon_data$Regulon_size >= threshold)
pct_ge <- mean(regulon_data$Regulon_size >= threshold)

# Select top regulons above threshold for labeling
lab_dat <- regulon_data %>%
  filter(Regulon_size >= threshold) %>%
  slice_max(Regulon_size, n = topN_to_label, with_ties = FALSE)

# Histogram + density + rug plot with threshold line
p <- ggplot(regulon_data, aes(x = Regulon_size)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.35) +
  geom_density(linewidth = 1.0) +
  geom_rug(sides = "b", alpha = 0.25) +
  geom_vline(xintercept = threshold, color = "red", linetype = "dashed") +
  annotate(
    "label",
    x = threshold, y = Inf, vjust = 1.2, hjust = 0,
    label = paste0("\u2265 ", threshold, ": ", n_ge, " (", percent(pct_ge), ")"),
    color = "red"
  ) +
  # Highlight regulons above threshold along x-axis
  geom_point(
    data = regulon_data %>% filter(Regulon_size >= threshold),
    aes(x = Regulon_size, y = 0),
    inherit.aes = FALSE, size = 1.6, color = "red", alpha = 0.7
  ) +
  # Label top regulons with repelled text
  ggrepel::geom_text_repel(
    data = lab_dat,
    aes(x = Regulon_size, y = 0, label = Regulon),
    inherit.aes = FALSE,
    nudge_y = 0.03,
    min.segment.length = 0,
    max.overlaps = Inf,
    size = 3.2
  ) +
  labs(
    title = "Regulon size distribution",
    x     = "Regulon size (# targets)",
    y     = "Density"
  ) +
  theme_minimal(base_size = 12)

# Optional: log10 x-axis for long-tailed distributions
if (use_log10_x) {
  p <- p + scale_x_continuous(trans = "log10", breaks = c(1, 3, 10, 30, 100, 300))
}

print(p)
ggsave("regulon_size_density.png", p, width = 8.5, height = 5.2, dpi = 300)

# Summary statistics for regulons with >= 10 genes
regulon_data2 <- regulon_data %>% filter(Regulon_size >= 10)
summary(regulon_data2)
length(regulon_data2$Regulon_size)

# ============================================================
# Part 2: PBMC UMAP with custom cell type color palette
# ============================================================

# Cell type color palette (14 broad PBMC cell types)
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

# Extended 27-color palette for fine-grained annotations (annotation_broad2)
my_colors <- c(
  "#9E3B30", "#D6A18C", "#C07E6B", "#A63D36", "#B45D4A", "#8E2C28", "#C5794B", "#A02D28",
  "#D96CA0", "#9429A4", "#7F7F99", "#6BAF4E", "#C2B78C", "#3F6843", "#4A774D", "#588B55",
  "#6AAE6B", "#8DC890", "#31693E", "#E6D86B", "#F3D447", "#E2A957", "#6699CC", "#286BB0",
  "#1F4489", "#6B7FA8", "#C3D5DE"
)
names(my_colors) <- levels(obj$annotation_broad2)[1:length(my_colors)]

p <- DimPlot(
  obj,
  reduction  = "umap",
  group.by   = "annotation_broad2",
  cols       = my_colors,
  label      = TRUE,
  repel      = TRUE,
  label.size = 4,
  pt.size    = 0.35,
  shuffle    = TRUE
) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 16)
  )

ggsave("PBMC_umap_annotation_broad2_custom_palette.png",
       p, device = ragg::agg_png, width = 8, height = 8, units = "in", dpi = 600, bg = "white")
