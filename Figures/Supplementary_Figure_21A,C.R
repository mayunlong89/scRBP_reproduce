
# --- 2026-01-14
# Developmental regulon visualization: violin, scatter, and top-N bar plots
# by trend cluster and overall Fisher test significance

# ============================================================
# Load required packages
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# ============================================================
# 0) Input / Output
# ============================================================
fisher_csv   <- "DE_overall_stage_Fisher.csv"
clusters_csv <- "RegulonTrends_global_v3.clusters.csv"
out_prefix   <- "devRegulon_trendCluster"

# ============================================================
# 1) Read data and merge Fisher results with cluster assignments
# ============================================================
de0 <- fread(fisher_csv, data.table = FALSE, check.names = TRUE)
cl0 <- fread(clusters_csv, data.table = FALSE, check.names = TRUE)

# Report proportion of significant devRegulons
n_sig <- sum(de0$FDR_overall < 0.05, na.rm = TRUE)
message("[INFO] Significant devRegulons (FDR < 0.05): ", n_sig, " / ", nrow(de0),
        " (", round(100 * n_sig / nrow(de0), 1), "%)")

# Validate required columns
needed_de <- c("regulon", "n_contrasts", "FDR_overall", "avg_logFC")
needed_cl <- c("regulon", "cluster")
stopifnot(all(needed_de %in% colnames(de0)), all(needed_cl %in% colnames(cl0)))

de_df <- de0 %>%
  inner_join(cl0, by = "regulon") %>%
  mutate(
    cluster     = as.integer(cluster),
    FDR_overall = as.numeric(FDR_overall),
    avg_logFC   = as.numeric(avg_logFC),
    n_contrasts = as.numeric(n_contrasts)
  ) %>%
  filter(is.finite(FDR_overall), FDR_overall > 0, is.finite(avg_logFC), is.finite(n_contrasts))

if (nrow(de_df) == 0) stop("After merging, de_df has 0 rows. Check regulon names match.")

# ============================================================
# 2) Color palette for 10 trend clusters
# ============================================================
# Display order matching the trend heatmap (top -> bottom)
cluster_order_in_fig <- c(3, 6, 9, 4, 1, 8, 10, 7, 5, 2)

cluster_pal10 <- c(
  "#2C7BB6",  # blue
  "#D7191C",  # red
  "#FDAE61",  # orange
  "#1A9641",  # green
  "#ABD9E9",  # light blue
  "#F46D43",  # orange-red
  "#7B3294",  # purple
  "#4D4D4D",  # dark gray
  "#66C2A5",  # teal-green
  "#E6AB02"   # mustard
)

# Map colors to cluster IDs in display order
pal_map <- setNames(cluster_pal10, as.character(cluster_order_in_fig))

# Warn if any clusters in data are missing from the palette
miss <- setdiff(sort(unique(de_df$cluster)), cluster_order_in_fig)
if (length(miss) > 0) {
  warning("Clusters in data but not in palette: ", paste(miss, collapse = ", "))
}

# Define sorted cluster levels (used for factor ordering in plots)
klev <- sort(unique(de_df$cluster))

# ============================================================
# 3) Figure 1: Violin plot of significance by cluster
# ============================================================
plot_df1 <- de_df %>%
  mutate(
    cluster   = factor(cluster, levels = klev),
    neglogFDR = -log10(FDR_overall)
  )

p1 <- ggplot(plot_df1, aes(x = cluster, y = neglogFDR, fill = cluster)) +
  geom_violin(width = 0.95, trim = TRUE, color = NA, alpha = 0.85) +
  geom_boxplot(width = 0.18, outlier.shape = NA, color = "grey15", linewidth = 0.35) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0),
    size = 0.55, alpha = 0.25, color = "grey10"
  ) +
  scale_fill_manual(values = pal_map, guide = "none") +
  labs(
    title = "Significance distribution by regulon trend cluster",
    x     = "Trend cluster",
    y     = expression(-log[10]("FDR (overall)"))
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

ggsave(paste0(out_prefix, ".violin.pdf"), p1, width = 7.6, height = 4.8, useDingbats = FALSE)
ggsave(paste0(out_prefix, ".violin.png"), p1, width = 7.6, height = 4.8, dpi = 300)

# ============================================================
# 4) Figure 2: Scatter plot (avg logFC vs -log10 FDR)
# ============================================================
plot_df2 <- de_df %>%
  mutate(
    cluster   = factor(cluster, levels = klev),
    neglogFDR = -log10(FDR_overall),
    direction = factor(ifelse(avg_logFC >= 0, "Up", "Down"), levels = c("Down", "Up"))
  )

thr  <- -log10(0.05)
topN_label <- 12

lab_df <- plot_df2 %>% arrange(FDR_overall) %>% slice_head(n = topN_label)

p2 <- ggplot(plot_df2, aes(x = avg_logFC, y = neglogFDR)) +
  geom_hline(yintercept = thr, linetype = "dashed", linewidth = 0.6, color = "grey35") +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.45, color = "grey75") +
  geom_point(aes(color = cluster, shape = direction, size = n_contrasts), alpha = 0.88) +
  scale_color_manual(values = pal_map, name = "Trend cluster") +
  scale_shape_manual(values = c(16, 17), name = "Direction") +
  scale_size_continuous(range = c(1.6, 5.5), breaks = c(1, 3, 5, 7),
                        name = "# significant\ncontrasts") +
  ggrepel::geom_text_repel(
    data = lab_df, aes(label = regulon),
    size = 3.4, min.segment.length = 0, box.padding = 0.35,
    point.padding = 0.25, max.overlaps = Inf,
    segment.color = "grey55", segment.size = 0.35
  ) +
  labs(
    title = "Development-associated regulons (overall Fisher test)",
    x     = "Average logFC across stage contrasts",
    y     = expression(-log[10]("FDR (overall)"))
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 18, hjust = 0),
    axis.title    = element_text(face = "bold"),
    legend.title  = element_text(face = "bold"),
    legend.box    = "vertical",
    legend.position = "right"
  )

ggsave(paste0(out_prefix, ".scatter.pdf"), p2, width = 7.4, height = 5.8, useDingbats = FALSE)
ggsave(paste0(out_prefix, ".scatter.png"), p2, width = 7.4, height = 5.8, dpi = 300)

# ============================================================
# 5) Figure 3: Top-N bar plot by overall Fisher FDR
# ============================================================
de_df <- de_df %>%
  mutate(
    neglog10FDR = -log10(pmax(FDR_overall, 1e-300)),
    direction   = factor(ifelse(avg_logFC >= 0, "Up", "Down"), levels = c("Down", "Up"))
  )

topN_bar <- 50

top_df <- de_df %>%
  arrange(FDR_overall) %>%
  slice_head(n = topN_bar) %>%
  mutate(regulon = factor(regulon, levels = rev(regulon)))

dir_pal <- c("Down" = "#66C2A5", "Up" = "#D7191C")

p_top <- ggplot(top_df, aes(x = regulon, y = neglog10FDR)) +
  geom_col(aes(fill = direction), color = "black", linewidth = 0.25, width = 0.85) +
  scale_fill_manual(values = dir_pal) +
  coord_flip() +
  labs(
    x     = NULL,
    y     = expression(-log[10]("FDR_overall")),
    fill  = NULL,
    title = paste0("Top ", topN_bar, " devRegulons by overall FDR")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title      = element_text(face = "bold"),
    axis.text.y     = element_text(size = 9)
  )

ggsave(paste0(out_prefix, ".top", topN_bar, ".bar.pdf"),
       p_top, limitsize = FALSE,
       width = 7.2, height = 0.28 * topN_bar + 2.0, useDingbats = FALSE)

cat("DONE:\n",
    paste0(out_prefix, ".violin.pdf\n"),
    paste0(out_prefix, ".scatter.pdf\n"),
    paste0(out_prefix, ".top", topN_bar, ".bar.pdf\n"))
