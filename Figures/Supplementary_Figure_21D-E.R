
#2026-01-16

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# =========================
# Inputs
# =========================
rss_file   <- "../04_scRBP_ras_ct_fetalbrain/10X_fetal_brain_590Kcells_ras_ct_all4regions.symbol_main_celltype_rss.csv"
trend_file <- "RegulonTrends_global_v3.trend01_matrix.csv"   # regulon x g1..g80
cluster_file <- "RegulonTrends_global_v3.clusters.csv"
out_dir    <- "TrendCluster_vs_Celltype_Method2_weightedTrends"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 1) Read trend matrix
# =========================
trend <- fread(trend_file, data.table=FALSE, check.names=FALSE)
if (!("regulon" %in% colnames(trend))) colnames(trend)[1] <- "regulon"
grid_cols <- grep("^g\\d+$", colnames(trend), value=TRUE)
stopifnot(length(grid_cols) >= 10)

trend_long <- trend %>%
  select(regulon, all_of(grid_cols)) %>%
  pivot_longer(all_of(grid_cols), names_to="grid", values_to="trend01") %>%
  mutate(grid_i = as.integer(sub("^g","",grid)),
         trend01 = as.numeric(trend01))

# =========================
# 2) Read RSS
# =========================
rss <- fread(rss_file, data.table=FALSE, check.names=FALSE)
colnames(rss)[1] <- "cell_type"
rss_long <- rss %>%
  pivot_longer(-cell_type, names_to="regulon", values_to="RSS") %>%
  mutate(RSS = as.numeric(RSS))

# =========================
# 3) Read clusters
# =========================
cl <- fread(cluster_file, data.table=FALSE)
stopifnot(all(c("regulon","cluster") %in% colnames(cl)))
cl$cluster <- as.integer(cl$cluster)

# overlap
common_regs <- Reduce(intersect, list(unique(trend_long$regulon),
                                      unique(rss_long$regulon),
                                      unique(cl$regulon)))
if (length(common_regs) < 10) stop("Too few overlap regulons among trend/RSS/cluster.")
trend_long <- trend_long %>% filter(regulon %in% common_regs)
rss_long   <- rss_long   %>% filter(regulon %in% common_regs)
cl         <- cl         %>% filter(regulon %in% common_regs)

# =========================
# 4) Build weighted celltype trend curve
# =========================
# normalize weights within each cell_type so sum(w)=1 (prevents celltype scale differences)
w <- rss_long %>%
  group_by(cell_type) %>%
  mutate(w = RSS / sum(RSS, na.rm=TRUE)) %>%
  ungroup() %>%
  select(cell_type, regulon, w)

celltype_trend <- trend_long %>%
  inner_join(w, by="regulon") %>%
  group_by(cell_type, grid_i) %>%
  summarise(weighted_trend = sum(w * trend01, na.rm=TRUE), .groups="drop")

fwrite(celltype_trend, file.path(out_dir, "celltype_weighted_trend.curve.csv"))

p <- ggplot(celltype_trend, aes(x=grid_i, y=weighted_trend, color=cell_type)) +
  geom_line(linewidth=1.0) +
  theme_classic(base_size = 12) +
  labs(x="Developmental grid", y="Weighted trend (RSS-weighted)",
       title="Celltype-level weighted regulon trend")
ggsave(file.path(out_dir, "celltype_weighted_trend.curves.pdf"), p, width=10, height=5)

# =========================
# 5) Cluster contribution per celltype (which trend cluster drives a celltype curve?)
# =========================
trend_cluster_long <- trend_long %>%
  inner_join(cl, by="regulon")

celltype_cluster_curve <- trend_cluster_long %>%
  inner_join(w, by="regulon") %>%
  group_by(cell_type, cluster, grid_i) %>%
  summarise(weighted_trend = sum(w * trend01, na.rm=TRUE), .groups="drop")

fwrite(celltype_cluster_curve, file.path(out_dir, "celltype_weighted_trend.byCluster.csv"))

# quick visualization: facet by celltype, color by cluster
p2 <- ggplot(celltype_cluster_curve,
             aes(x=grid_i, y=weighted_trend, color=factor(cluster))) +
  geom_line(linewidth=0.9) +
  facet_wrap(~ cell_type, scales="free_y") +
  theme_classic(base_size = 11) +
  labs(x="Developmental grid", y="Weighted trend",
       color="Trend cluster",
       title="Cluster-resolved contributions to each celltype weighted trend")
ggsave(file.path(out_dir, "celltype_weighted_trend.byCluster.pdf"), p2, width=12, height=8)

cat("DONE. Outputs in:", out_dir, "\n")






# celltype_cluster_curve 需要至少包含：
# grid_i, cluster, cell_type, weighted_trend

library(ggplot2)
library(dplyr)

celltype_cluster_curve <- celltype_cluster_curve %>%
  mutate(
    cluster   = factor(cluster),
    cell_type = factor(cell_type)
  )

p_trend_by_cluster <- ggplot(
  celltype_cluster_curve,
  aes(x = grid_i, y = weighted_trend, group = cell_type, color = cell_type)
) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  facet_wrap(~ cluster, ncol = 2, scales = "free_y") +
  theme_classic(base_size = 11) +
  labs(
    x = "Developmental grid",
    y = "Cluster contribution (RSS-weighted)",
    color = "Cell type",
    title = "Cell-type trajectories within each regulon trend cluster"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "cluster_contribution.by_cluster.pdf"),
       p_trend_by_cluster, width = 12, height = 8, useDingbats = FALSE)

cat("DONE. Outputs in:", out_dir, "\n")






library(ggplot2)
library(dplyr)

topN <- 6  # 每个 cluster 保留前 N 个 cell types（按 AUC/mean 贡献）

# 计算每个 cluster 内：每个 cell type 的总体贡献（用 mean 或 AUC 都行）
rank_tbl <- celltype_cluster_curve %>%
  group_by(cluster, cell_type) %>%
  summarize(score = mean(weighted_trend, na.rm = TRUE), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(order_by = score, n = topN, with_ties = FALSE) %>%
  ungroup()

df_top <- celltype_cluster_curve %>%
  inner_join(rank_tbl, by = c("cluster", "cell_type")) %>%
  mutate(
    cluster   = factor(cluster),
    cell_type = factor(cell_type)
  )

p_top <- ggplot(
  df_top,
  aes(x = grid_i, y = weighted_trend, group = cell_type, color = cell_type)
) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~ cluster, ncol = 2, scales = "fixed") +
  theme_classic(base_size = 11) +
  labs(
    x = "Developmental grid",
    y = "Cluster contribution (RSS-weighted)",
    color = "Top cell types",
    title = paste0("Top cell-type trajectories within each trend cluster (Top ", topN, ")")
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(out_dir, paste0("cluster_contribution.by_cluster.top", topN, ".pdf")),
       p_top, width = 11, height = 7.5, useDingbats = FALSE)

cat("DONE. Outputs in:", out_dir, "\n")




library(ggplot2)
library(dplyr)

celltype_cluster_curve <- celltype_cluster_curve %>%
  mutate(
    cluster   = factor(cluster),
    cell_type = factor(cell_type)
  )

p_trend_by_cluster <- ggplot(
  celltype_cluster_curve,
  aes(x = grid_i, y = weighted_trend, group = cell_type, color = cell_type)
) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  facet_wrap(~ cluster, ncol = 2, scales = "free_y") +   # <- 关键改动
  theme_classic(base_size = 11) +
  labs(
    x = "Developmental grid",
    y = "Cluster contribution (RSS-weighted)",
    color = "Cell type",
    title = "Cell-type trajectories within each regulon trend cluster"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "cluster_contribution.by_cluster.freeY.pdf"),
       p_trend_by_cluster, width = 12, height = 8, useDingbats = FALSE)

cat("DONE. Outputs in:", out_dir, "\n")



library(ggplot2)
library(dplyr)

# 每个 cluster 内按 1%~99% 分位数裁剪
df_plot <- celltype_cluster_curve %>%
  mutate(cluster = factor(cluster), cell_type = factor(cell_type)) %>%
  group_by(cluster) %>%
  mutate(
    lo = quantile(weighted_trend, 0.01, na.rm = TRUE),
    hi = quantile(weighted_trend, 0.99, na.rm = TRUE),
    weighted_trend_clip = pmin(pmax(weighted_trend, lo), hi)
  ) %>%
  ungroup()

p_trend_by_cluster <- ggplot(
  df_plot,
  aes(x = grid_i, y = weighted_trend_clip, group = cell_type, color = cell_type)
) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  facet_wrap(~ cluster, ncol = 2, scales = "free_y") +
  theme_classic(base_size = 11) +
  labs(
    x = "Developmental grid",
    y = "Cluster contribution (RSS-weighted, clipped)",
    color = "Cell type",
    title = "Cell-type trajectories within each regulon trend cluster"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "cluster_contribution.by_cluster.freeY.clipped.pdf"),
       p_trend_by_cluster, width = 12, height = 8, useDingbats = FALSE)

cat("DONE. Outputs in:", out_dir, "\n")




