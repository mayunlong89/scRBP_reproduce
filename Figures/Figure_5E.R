
#2026-01-16

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})

# =========================
# Inputs
# =========================
rss_file     <- "../04_scRBP_ras_ct_fetalbrain/10X_fetal_brain_590Kcells_ras_ct_all4regions.symbol_main_celltype_rss.csv"
cluster_file <- "RegulonTrends_global_v3.clusters.csv"
out_dir      <- "TrendCluster_vs_Celltype_Method1"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 1) Read RSS (celltype x regulon)
# =========================
rss <- fread(rss_file, data.table = FALSE, check.names = FALSE)
stopifnot(ncol(rss) >= 2)
colnames(rss)[1] <- "cell_type"

# long
rss_long <- rss %>%
  pivot_longer(-cell_type, names_to = "regulon", values_to = "RSS") %>%
  mutate(RSS = as.numeric(RSS))

# =========================
# 2) Read regulon -> cluster
# =========================
cl <- fread(cluster_file, data.table = FALSE)
stopifnot(all(c("regulon","cluster") %in% colnames(cl)))
cl$cluster <- as.integer(cl$cluster)

# keep overlap regulons
common_regs <- intersect(unique(rss_long$regulon), unique(cl$regulon))
if (length(common_regs) < 10) stop("Too few overlap regulons between RSS and clusters.")
rss_long <- rss_long %>% filter(regulon %in% common_regs)
cl       <- cl       %>% filter(regulon %in% common_regs)

# merge
df <- rss_long %>% inner_join(cl, by="regulon")

# =========================
# 3) Winner celltype per regulon (max RSS)
# =========================
winner <- df %>%
  group_by(regulon) %>%
  slice_max(order_by = RSS, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(regulon, winner_celltype = cell_type, winner_RSS = RSS) %>%
  inner_join(cl, by="regulon")

fwrite(winner, file.path(out_dir, "regulon_winner_celltype.csv"))

# cluster-wise winner composition
comp <- winner %>%
  count(cluster, winner_celltype, name="n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

fwrite(comp, file.path(out_dir, "cluster_winner_celltype.composition.csv"))

p1 <- ggplot(comp, aes(x=factor(cluster), y=prop, fill=winner_celltype)) +
  geom_col(width=0.85) +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  labs(x="Trend cluster", y="Winner celltype proportion") +
  theme_classic(base_size = 12)
ggsave(file.path(out_dir, "cluster_winner_celltype.stackedbar.pdf"), p1, width=9, height=5)

# =========================
# 4) Cluster x celltype mean RSS + zscore enrichment
# =========================
mat_mean <- df %>%
  group_by(cluster, cell_type) %>%
  summarise(mean_RSS = mean(RSS, na.rm=TRUE),
            median_RSS = median(RSS, na.rm=TRUE),
            .groups="drop")

# z-score within each celltype across clusters (avoid "neuron regulons too many" bias)
mat_z <- mat_mean %>%
  group_by(cell_type) %>%
  mutate(z = (mean_RSS - mean(mean_RSS)) / sd(mean_RSS)) %>%
  ungroup()

fwrite(mat_mean, file.path(out_dir, "cluster_celltype.meanRSS.csv"))
fwrite(mat_z,    file.path(out_dir, "cluster_celltype.zscore.csv"))

# wide for heatmap
mean_wide <- mat_mean %>%
  dplyr::select(cluster, cell_type, mean_RSS) %>%
  pivot_wider(names_from = cell_type, values_from = mean_RSS) %>%
  arrange(cluster)

z_wide <- mat_z %>%
  dplyr::select(cluster, cell_type, z) %>%
  pivot_wider(names_from = cell_type, values_from = z) %>%
  arrange(cluster)

mean_mat <- as.matrix(mean_wide[,-1, drop=FALSE]); rownames(mean_mat) <- paste0("C", mean_wide$cluster)
z_mat    <- as.matrix(z_wide[,-1, drop=FALSE]);    rownames(z_mat)    <- paste0("C", z_wide$cluster)

pdf(file.path(out_dir, "heatmap.cluster_celltype.meanRSS.pdf"), width=10, height=5.5)
pheatmap(mean_mat, cluster_rows=T, cluster_cols=TRUE, border_color=NA,
         main="Mean RSS per trend cluster")
dev.off()

pdf(file.path(out_dir, "heatmap.cluster_celltype.zscore.pdf"), width=10, height=5.5)
pheatmap(z_mat, cluster_rows=T, cluster_cols=TRUE, border_color=NA,
         main="Celltype-enrichment z-score (within celltype across clusters)")
dev.off()



# ---- z-score heatmap (diverging; centered at 0) ----
z_lim <- 3  # 常用截断，避免极端值把颜色拉爆
z_plot <- pmax(pmin(z_mat, z_lim), -z_lim)

bk <- seq(-z_lim, z_lim, length.out = 101)
col_fun <- colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(100)  # 蓝-白-红

pdf(file.path(out_dir, "heatmap.cluster_celltype.zscore_min_max.pdf"), width=10, height=5.5)
pheatmap(
  z_plot,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  color = col_fun,
  breaks = bk,
  main = "Cell-type enrichment across trend clusters (z-score)",
  fontsize_row = 10,
  fontsize_col = 10
)
dev.off()





# =========================
# 5) Top celltypes per cluster (table)
# =========================
top_ct <- mat_mean %>%
  group_by(cluster) %>%
  arrange(desc(mean_RSS)) %>%
  slice_head(n=5) %>%
  ungroup()

fwrite(top_ct, file.path(out_dir, "Top5_celltypes_per_cluster.byMeanRSS.csv"))

cat("DONE. Outputs in:", out_dir, "\n")
