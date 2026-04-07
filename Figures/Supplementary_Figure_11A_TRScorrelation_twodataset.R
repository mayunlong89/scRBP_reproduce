##--------Figure--3--related-Supplementary Figure S11A---correlation

library(tidyverse)
library(stringr)
library(scales)

# 1) Read in and tidy data (columns required: ID, cell_type, RBP, TRS_6.5K, TRS_8.7K)
df <- read.csv("Figure_2_Summary_monocyte_count_correlation_TRS.csv", check.names = FALSE, stringsAsFactors = FALSE)

cor.test(df$TRS_6.5K,df$TRS_8.7K)


# Extract region from ID
df <- df %>%
  mutate(
    region    = str_extract(ID, "(5UTR|3UTR|CDS|Introns?)$"),
    region    = recode(region, "Intron" = "Introns"),
    region    = factor(region, levels = c("5UTR","3UTR","CDS","Introns")),
    cell_type = factor(cell_type, levels = c("Monocytes","T_cells"))
  ) %>%
  drop_na(`TRS_6.5K`, `TRS_8.7K`)  # Keep only paired observations

# 2) Compute correlation for each cell type, for in-plot annotation
spearman_p <- function(x, y) suppressWarnings(cor.test(x, y, method = "spearman")$p.value)
stats_ct <- df %>%
  group_by(cell_type) %>%
  summarise(
    n    = n(),
    rho  = suppressWarnings(cor(`TRS_6.5K`, `TRS_8.7K`, method = "spearman")),
    pval = spearman_p(`TRS_6.5K`, `TRS_8.7K`),
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(pval, "BH"),
         label = sprintf("%s: ρ=%.2f, n=%d, FDR=%.3g",
                         as.character(cell_type), rho, n, FDR))

# 3) Colors/shapes (consistent with the red/teal two-color scheme used in the UMAP; adjust as needed)
cols_ct <- c(Monocytes = "#A84B4B", T_cells = "#53A6A3")
shapes  <- c("5UTR"=21, "3UTR"=22, "CDS"=24, "Introns"=23)  # 21/22/23/24 allow fill

# 4) Combined main plot: color = cell type, shape = region; one fitted line per cell type
p_all <- ggplot(df, aes(`TRS_6.5K`, `TRS_8.7K`)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(fill = cell_type, shape = region),
             size = 2.2, alpha = 0.9, color = "grey25", stroke = 0.35) +
  stat_smooth(aes(color = cell_type), method = "lm", se = FALSE, linewidth = 0.9) +
  scale_fill_manual(values = cols_ct, name = "Cell type") +
  scale_color_manual(values = cols_ct, guide = "none") +
  scale_shape_manual(values = shapes, name = "mRNA region") +
  coord_equal() +
  labs(
    title = "TRS correlation between 6.5K and 8.7K",
    subtitle = "Color = cell type; shape = mRNA region; dashed = y = x",
    x = "TRS (6.5K)", y = "TRS (8.7K)"
  ) +
  # Annotate rho and n for each cell type inside the plot
  geom_label(
    data = stats_ct,
    aes(x = -Inf, y = Inf, label = label, fill = cell_type),
    inherit.aes = FALSE,
    hjust = -0.02, vjust = 1.1, color = "white",
    alpha = 0.9, label.size = 0
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(hjust = 0),
    legend.position = "right"
  )

ggsave("Figure_2_TRS_correlation_combined_main.pdf", p_all, width = 8.5, height = 6)
p_all




# Side-by-side facets: one panel per cell type
p_facets <- ggplot(df, aes(`TRS_6.5K`, `TRS_8.7K`)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(fill = cell_type, shape = region),
             size = 2.2, alpha = 0.9, color = "grey25", stroke = 0.35) +
  stat_smooth(method = "lm", se = FALSE, linewidth = 0.9,
              aes(color = cell_type)) +
  scale_fill_manual(values = cols_ct, guide = "none") +
  scale_color_manual(values = cols_ct, guide = "none") +
  scale_shape_manual(values = shapes, name = "mRNA region") +
  coord_equal() +
  labs(x = "TRS (6.5K)", y = "TRS (8.7K)") +
  facet_wrap(~cell_type, nrow = 1, scales = "fixed",
             labeller = labeller(cell_type = c(Monocytes = "Monocytes",
                                               T_cells   = "T cells"))) +
  # Annotate rho, n, and FDR in the corner of each panel
  geom_label(data = stats_ct,
             aes(x = -Inf, y = Inf, label = label, fill = cell_type),
             inherit.aes = FALSE, hjust = -0.02, vjust = 1.1,
             color = "white", alpha = 0.9, label.size = 0) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

ggsave("Figure_2_TRS_correlation_by_celltype_facets.pdf", p_facets,
       width = 7.0, height = 3.6)
p_facets
