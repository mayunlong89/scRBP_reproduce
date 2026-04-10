
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
})

# -------------------------
# Input
# -------------------------

trend_file    <- "RegulonTrends_global_v3.trend01_matrix.csv"   # regulon x g1..g80 (0-1 scaled)
cluster_file  <- "RegulonTrends_global_v3.clusters.csv"         # columns: regulon, cluster
out_pdf     <- "RegulonTrends_cluster_curves.pdf"

# grids / stage bins

# stage settings (你现在是 Stage1~Stage9)
stage_labels <- paste0("Stage", 1:9)

# -------------------------
# Read
# -------------------------
trend <- fread(trend_file, data.table = FALSE, check.names = FALSE)
# trend 第一列通常是 regulon 名（可能列名是 "" 或 "V1"）
if (!("regulon" %in% colnames(trend))) {
  colnames(trend)[1] <- "regulon"
}
cl <- fread(cluster_file, data.table = FALSE)
stopifnot(all(c("regulon","cluster") %in% colnames(cl)))

# 只保留在 cluster 表里的 regulon
trend <- trend %>% inner_join(cl, by = "regulon")

# 找到 grid 列（g1,g2,...）
grid_cols <- grep("^g\\d+$", colnames(trend), value = TRUE)
stopifnot(length(grid_cols) >= 10)

n_grid <- length(grid_cols)

# long format
df_long <- trend %>%
  select(regulon, cluster, all_of(grid_cols)) %>%
  pivot_longer(cols = all_of(grid_cols), names_to = "grid", values_to = "y") %>%
  mutate(
    grid_i = as.integer(sub("^g", "", grid)),
    cluster = as.factor(cluster)
  )

# -------------------------
# Make stage bins on a continuous grid
# 80 grids -> 9 stages => 每个 stage 占大约 80/9 个 grid
# 用等分切分（你之前也是“等分 intervals”思路）
# -------------------------
# breaks: 0..n_grid 切成 9 段
stage_breaks <- round(seq(0, n_grid, length.out = length(stage_labels) + 1))
stage_breaks[1] <- 0
stage_breaks[length(stage_breaks)] <- n_grid

# stage id for each grid_i
grid2stage <- function(i) {
  # 返回 1..9
  findInterval(i, vec = stage_breaks, rightmost.closed = TRUE)
}
df_long <- df_long %>%
  mutate(stage_id = pmax(1, pmin(length(stage_labels), grid2stage(grid_i))),
         stage = factor(stage_labels[stage_id], levels = stage_labels))

# stage vertical lines at boundaries (1..n_grid)
vlines <- data.frame(x = stage_breaks[-c(1, length(stage_breaks))])

# stage label positions (每段中心)
stage_centers <- floor((stage_breaks[-1] + stage_breaks[-length(stage_breaks)] + 1) / 2)
stage_axis_df <- data.frame(
  stage = factor(stage_labels, levels = stage_labels),
  x = stage_centers
)

# -------------------------
# Cluster summary: mean ± sd per grid
# -------------------------
summ <- df_long %>%
  group_by(cluster, grid_i) %>%
  summarize(
    mean = mean(y, na.rm = TRUE),
    sd   = sd(y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = pmax(0, mean - sd),
    ymax = pmin(1, mean + sd)
  )

# -------------------------
# Nice cluster colors (接近 Cell 那种“每个 trend 一种色”)
# 8 clusters 给 8 个清晰颜色（你也可以换）
# -------------------------
cluster_levels <- sort(unique(as.integer(as.character(df_long$cluster))))
cluster_levels <- as.character(cluster_levels)

cluster_pal <- c(
  "#2C7BB6", # blue
  "#D7191C", # red
  "#FDAE61", # orange
  "#1A9641", # green
  "#ABD9E9", # light blue
  "#F46D43", # orange-red
  "#7B3294", # purple
  "#4D4D4D", # dark gray
  "#66A61E", # olive green (new)
  "#E6AB02"  # mustard (new)
)

names(cluster_pal) <- cluster_levels[seq_len(min(length(cluster_levels), length(cluster_pal)))]

# 如果你的 cluster 编号不是 1..8，补齐映射
missing_cols <- setdiff(levels(df_long$cluster), names(cluster_pal))
if (length(missing_cols) > 0) {
  # 给缺失的随便补一圈灰色系
  cluster_pal[missing_cols] <- rep("#4D4D4D", length(missing_cols))
}

# -------------------------
# Build plot (Cell-like)
# -------------------------
p <- ggplot() +
  # faint member curves
  geom_line(
    data = df_long,
    aes(x = grid_i, y = y, group = regulon),
    color = "grey75", linewidth = 0.25, alpha = 0.35
  ) +
  # ribbon = mean ± sd (colored)
  geom_ribbon(
    data = summ,
    aes(x = grid_i, ymin = ymin, ymax = ymax, fill = cluster),
    alpha = 0.35, color = NA
  ) +
  # mean line (colored + thicker)
  geom_line(
    data = summ,
    aes(x = grid_i, y = mean, color = cluster),
    linewidth = 1.2, lineend = "round"
  ) +
  # stage split lines
  geom_vline(
    data = vlines,
    aes(xintercept = x),
    linewidth = 0.35, color = "grey65"
  ) +
  # facet by cluster
  facet_wrap(~ cluster, ncol = 1, scales = "fixed") +
  scale_x_continuous(
    breaks = stage_axis_df$x,
    labels = as.character(stage_axis_df$stage),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.02, 0.02)) +
  scale_fill_manual(values = cluster_pal) +
  scale_color_manual(values = cluster_pal) +
  guides(fill = "none", color = "none") +
  labs(
    x = NULL,
    y = "Normalized regulon activity",
    title = "Regulon trend clusters across development"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 0, vjust = 0.8),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.25, "lines")
  )

ggsave(out_pdf, p, width = 10.5, height = 7.5, useDingbats = FALSE)
cat("Wrote:", out_pdf, "\n")

