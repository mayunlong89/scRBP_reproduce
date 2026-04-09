
#2025-06-10

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(RColorBrewer)

# === Step 1: 批量读取所有 CSV ===
files <- list.files(path = "./", pattern = "*.csv", full.names = TRUE)

# 提取方法和区域信息
extract_info <- function(filename) {
  parts <- strsplit(basename(filename), "_")[[1]]
  list(
    method = gsub(".csv", "", parts[8]),
    cell_line = parts[4],
    region = parts[9]
  )
}

all_data <- lapply(files, function(file) {
  df <- read_csv(file)
  info <- extract_info(file)
  df$Method <- info$method
  df$CellLine <- info$cell_line
  df$Region <- info$region
  return(df)
}) %>% bind_rows()


# === Step 1.1: 每个方法 × 区域 × RBP 先保留 F1 最优的 motif（去除冗余 motif）===
best_motifs_per_method_region <- all_data %>%
  group_by(RBP, CellLine, Method, Region) %>%
  slice_max(order_by = F1, n = 1, with_ties = FALSE) %>%
  ungroup()

# === Step 2: 获取每个 RBP + CellLine 在 clustered 方法中 F1 最佳区域 ===
best_clustered_regions <- best_motifs_per_method_region %>%
  filter(Method == "clustered") %>%
  group_by(RBP, CellLine) %>%
  slice_max(order_by = F1, n = 1, with_ties = FALSE) %>%
  select(RBP, CellLine, Region) %>%
  rename(BestRegion = Region)

# === Step 3: 所有方法都按 clustered 的最佳区域筛选 ===
best_data <- best_motifs_per_method_region %>%
  inner_join(best_clustered_regions, by = c("RBP", "CellLine")) %>%
  filter(Region == BestRegion)

# === Step 4: 转为长格式准备画图 ===
plot_data <- best_data %>%
  select(Method, CellLine, F1, Precision, Recall) %>%
  pivot_longer(cols = c("F1", "Precision", "Recall"),
               names_to = "Metric", values_to = "Score")

# 设置统一颜色和方法顺序
plot_data$Method <- factor(plot_data$Method,
                           levels = c("clustered", "singleton", "archetype", "FIMO", "HOMER2"))
method_colors <- brewer.pal(5, "Set2")

# === Step 5: 画图 ===
# Boxplot without dots
p1 <- ggplot(plot_data, aes(x = Method, y = Score, fill = Method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  facet_wrap(~Metric, scales = "free_y", nrow = 3) +
  scale_fill_manual(values = method_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  labs(title = "Motif Method Comparison (Each RBP Uses Clustered's Best Region)",
       x = "Motif Scanning Method", y = "Score")

pdf("motif_method_comparison_boxplot_clustered_best_region.pdf", width = 10, height = 8)
print(p1)
dev.off() 

# Boxplot with jitter
p2 <- ggplot(plot_data, aes(x = Method, y = Score, fill = Method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.4, color = "black") +
  facet_wrap(~Metric, scales = "free_y", nrow = 3) +
  scale_fill_manual(values = method_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  labs(title = "Motif Method Comparison (Each RBP Uses Clustered's Best Region)",
       x = "Motif Scanning Method", y = "Score")

pdf("motif_method_comparison_boxplot_with_dots_classic_clustered_best_region.pdf", width = 10, height = 8)
print(p2)
dev.off()

