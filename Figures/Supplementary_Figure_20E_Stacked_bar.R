# --- Composition of developmental stages by cell types and subtypes
# Reproduces the stacked bar plots with stage dot annotations below x-axis

# ============================================================
# Load required packages
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(patchwork)
})

# ============================================================
# Load Seurat object (assumes GW_group, Type.v2, and subtype columns exist)
# ============================================================
# seu <- readRDS("your_seurat_object.rds")

# ============================================================
# Color palettes
# ============================================================

# Cell type colors (matching the figure)
celltype_cols <- c(
  "Div"               = "#F4A259",
  "Endothelial"       = "#E15759",
  "Excitatory Neuron" = "#8FBF88",
  "Inhibitory Neuron" = "#3A7D44",
  "IPC"               = "#66C2A5",
  "Microglia"         = "#76B7B2",
  "Newborn Neuron"    = "#A5D6A7",
  "Oligodendrocyte"   = "#59A14F",
  "OPC"               = "#4E79A7",
  "OPC.div"           = "#A0CBE8",
  "Pericyte"          = "#B07AA1",
  "Red blood cells"   = "#D37295",
  "RG"                = "#C79ECF",
  "RG.div"            = "#E8C8E0"
)

# Developmental stage colors (green gradient + postnatal red)
stage_cols <- c(
  "1: GW 6-10"        = "#1A5C1A",
  "2: GW 10-12"       = "#2E8B2E",
  "3: GW 12-15"       = "#6BBF6B",
  "4: GW 15-18"       = "#9CCB9C",
  "5: GW 18-21"       = "#B8DAB8",
  "6: GW 21-26"       = "#D0E7D0",
  "7: GW 26-41"       = "#E6F2E6",
  "8: 8 mo postnatal" = "#D7191C"
)

# ============================================================
# Helper function: stacked bar + stage dots
# ============================================================
plot_composition <- function(seu, group_col, stage_col = "GW_group",
                             color_pal, title = "") {

  # Build per-stage composition table
  meta <- seu@meta.data %>%
    dplyr::select(stage = !!sym(stage_col), group = !!sym(group_col)) %>%
    dplyr::filter(!is.na(stage), !is.na(group))

  comp <- meta %>%
    dplyr::count(stage, group) %>%
    dplyr::group_by(stage) %>%
    dplyr::mutate(pct = n / sum(n) * 100) %>%
    dplyr::ungroup()

  # Ensure stage ordering
  comp$stage <- factor(comp$stage, levels = levels(seu@meta.data[[stage_col]]))

  # Main stacked bar plot
  p_bar <- ggplot(comp, aes(x = stage, y = pct, fill = group)) +
    geom_col(width = 0.85, color = NA) +
    scale_fill_manual(values = color_pal, name = NULL) +
    scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    labs(
      title = title,
      x     = NULL,
      y     = "% of cells in developmental stage"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x       = element_blank(),
      axis.ticks.x       = element_blank(),
      plot.title         = element_text(size = 11, face = "bold", hjust = 0.5),
      legend.text        = element_text(size = 7),
      legend.key.size    = unit(3, "mm"),
      legend.position    = "bottom",
      legend.box         = "horizontal",
      plot.margin        = margin(5, 5, 0, 5)
    ) +
    guides(fill = guide_legend(ncol = 3, override.aes = list(size = 2)))

  # Stage dot annotation below x-axis
  stage_levels <- levels(comp$stage)
  dot_df <- data.frame(
    stage = factor(stage_levels, levels = stage_levels),
    y     = 0
  )

  p_dots <- ggplot(dot_df, aes(x = stage, y = y, color = stage)) +
    geom_point(size = 4) +
    scale_color_manual(values = stage_cols, guide = "none") +
    labs(x = "Developmental stage") +
    theme_void() +
    theme(
      axis.title.x = element_text(size = 10, margin = margin(t = 4)),
      plot.margin  = margin(0, 5, 5, 5)
    )

  # Combine bar + dots vertically
  p_bar / p_dots + plot_layout(heights = c(10, 1))
}

# ============================================================
# Plot 1: Composition by cell types
# ============================================================
p_celltype <- plot_composition(
  seu,
  group_col = "Type.v2",        # Your main cell type column
  stage_col = "GW_group",
  color_pal = celltype_cols,
  title     = "Composition of developmental stages by cell types"
)

# ============================================================
# Plot 2: Composition by cell subtypes
# ============================================================
# Generate a large palette for many subtypes (adjust if you have a curated one)
subtype_levels <- sort(unique(as.character(seu@meta.data$cell_subtype)))
n_sub <- length(subtype_levels)

# Use a combination of palettes to cover all subtypes
sub_pal <- setNames(
  colorRampPalette(
    c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
      "#A65628", "#F781BF", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
      "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02",
      "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#4E79A7",
      "#F28E2B", "#59A14F", "#B07AA1", "#76B7B2", "#EDC948", "#FF9DA7")
  )(n_sub),
  subtype_levels
)

p_subtype <- plot_composition(
  seu,
  group_col = "cell_subtype",   # Your cell subtype column name
  stage_col = "GW_group",
  color_pal = sub_pal,
  title     = "Composition of developmental stages by cell subtypes"
)

# ============================================================
# Combine side by side with shared stage legend
# ============================================================
# Stage legend (separate for the right side)
stage_legend_df <- data.frame(
  stage = factor(names(stage_cols), levels = names(stage_cols)),
  y = seq_along(stage_cols)
)

p_stage_legend <- ggplot(stage_legend_df, aes(x = 1, y = y, color = stage)) +
  geom_point(size = 4) +
  scale_color_manual(values = stage_cols, name = "Developmental stages") +
  theme_void() +
  theme(legend.text = element_text(size = 9))

# Combine everything
p_combined <- (p_celltype | p_subtype) +
  plot_layout(widths = c(1, 1.3))

# Save outputs
ggsave("composition_by_celltype.pdf", p_celltype, width = 7, height = 7)
ggsave("composition_by_subtype.pdf",  p_subtype,  width = 9, height = 8)
ggsave("composition_combined.pdf",    p_combined,  width = 16, height = 7)

ggsave("composition_combined.tiff", p_combined,
       width = 16, height = 7, units = "in", dpi = 600, compression = "lzw")
