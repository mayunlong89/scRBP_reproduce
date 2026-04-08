library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(ggplot2)
library(pheatmap)

set.seed(123)

# ==== STEP 1: 读取文件 ====
files <- list.files(pattern = "*2K_grn_network_.*csv")

networks <- lapply(files, function(file) {
  df <- read_csv(file, show_col_types = FALSE)
  df %>%
    group_by(RBP) %>%
    summarise(targets = list(unique(Gene)), .groups = "drop") %>%
    deframe()
})
names(networks) <- gsub("_2K_grn_network_IM0005_seed42\\.csv", "", files)

# ==== STEP 2: 构建 background genes ====
background_genes <- networks %>%
  purrr::map(~ unlist(.)) %>%
  unlist() %>%
  unique()

# ==== STEP 3: 定义函数 ====
jaccard <- function(set1, set2) {
  inter <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  if (union == 0) return(NA) else return(inter / union)
}

sample_null_modules <- function(size, n = 100) {
  replicate(n, sample(background_genes, size, replace = FALSE), simplify = FALSE)
}

# ==== STEP 4: 主循环 ====
strategies <- setdiff(names(networks), c("empiricalP0.01", "empiricalP0.05"))
results <- list()

for (method in strategies) {
  message("Processing: ", method)
  method_net <- networks[[method]]
  
  for (rbp in names(method_net)) {
    real_targets <- method_net[[rbp]]
    size <- length(real_targets)
    if (size == 0) next
    
    emp_targets_01 <- networks[["empiricalP0.01"]][[rbp]]
    emp_targets_05 <- networks[["empiricalP0.05"]][[rbp]]
    
    real_jaccard_01 <- jaccard(real_targets, emp_targets_01)
    real_jaccard_05 <- jaccard(real_targets, emp_targets_05)
    
    null_modules <- sample_null_modules(size, n = 100)
    null_jaccard_01 <- sapply(null_modules, function(nm) jaccard(nm, emp_targets_01))
    null_jaccard_05 <- sapply(null_modules, function(nm) jaccard(nm, emp_targets_05))
    
    # 加鲁棒性防止 Inf
    z_01 <- (real_jaccard_01 - mean(null_jaccard_01, na.rm = TRUE)) / max(sd(null_jaccard_01, na.rm = TRUE), 1e-6)
    p_01 <- (sum(null_jaccard_01 >= real_jaccard_01, na.rm = TRUE) + 1) / (length(null_jaccard_01) + 1)
    
    z_05 <- (real_jaccard_05 - mean(null_jaccard_05, na.rm = TRUE)) / max(sd(null_jaccard_05, na.rm = TRUE), 1e-6)
    p_05 <- (sum(null_jaccard_05 >= real_jaccard_05, na.rm = TRUE) + 1) / (length(null_jaccard_05) + 1)
    
    results[[length(results) + 1]] <- tibble(
      Method = method,
      RBP = rbp,
      Real_Jaccard_01 = real_jaccard_01,
      Zscore_P0.01 = z_01,
      EmpiricalP_P0.01 = p_01,
      Real_Jaccard_05 = real_jaccard_05,
      Zscore_P0.05 = z_05,
      EmpiricalP_P0.05 = p_05,
      ModuleSize = size
    )
  }
}

result_df <- bind_rows(results)
write_csv(result_df, "rbp_module_jaccard_zscore_summary.csv")

                              
                              
# ==== STEP 4.5: cluster 信息 ====                              
 # 方法分组信息
method_cluster <- tibble::tibble(
  Method = c(
    "top500perRBP", "percent75", "top25perTarget", "top1000perRBP", "percent50", "top50perTarget",
    "top300perRBP", "top200perRBP", "percent90",
    "top10perTarget", "top5perTarget", "percent95",
    "top100perRBP", "top3perTarget",
    "top50perRBP"
  ),
  Cluster = c(
    rep("Cluster 1", 6),
    rep("Cluster 2", 3),
    rep("Cluster 3", 3),
    rep("Cluster 4", 2),
    "Cluster 5"
  )
)

# 方法顺序
method_order <- method_cluster %>%
  arrange(Cluster) %>%
  pull(Method)

                              
                              
                              
                              
# ==== STEP 5: Boxplot ====
library(forcats)
# 调整顺序（上 → 下与 heatmap 保持一致）
method_order_reversed <- rev(method_order)

long_df <- result_df %>%
  pivot_longer(cols = c(Zscore_P0.01, Zscore_P0.05),
               names_to = "Metric", values_to = "Zscore") %>%
  left_join(method_cluster, by = "Method") %>%
  mutate(Method = factor(Method, levels = method_order_reversed))  # 修正顺序

# P < 0.01 boxplot
p1 <- ggplot(filter(long_df, Metric == "Zscore_P0.01"), aes(x = Method, y = Zscore, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.8, color = "black") +
  coord_flip() +
  theme_bw() +
  labs(title = "Z-score of RBP modules vs empiricalP0.01",
       y = "Z-score", x = "Method") +
  theme(axis.text.y = element_text(size = 9), legend.position = "right")

ggsave("rbp_module_zscore_boxplot_P0.01.pdf", p1, width = 9, height = 6)

# P < 0.05 boxplot
p2 <- ggplot(filter(long_df, Metric == "Zscore_P0.05"), aes(x = Method, y = Zscore, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.8, color = "black") +
  coord_flip() +
  theme_bw() +
  labs(title = "Z-score of RBP modules vs empiricalP0.05",
       y = "Z-score", x = "Method") +
  theme(axis.text.y = element_text(size = 9), legend.position = "right")

ggsave("rbp_module_zscore_boxplot_P0.05.pdf", p2, width = 9, height = 6)
           
                              
                              
# ==== STEP 6: Heatmap ====
# 平均 z-score，确保包含所有 methods
mean_z <- result_df %>%
  group_by(Method) %>%
  summarise(
    mean_Z_P0.01 = mean(Zscore_P0.01, na.rm = TRUE),
    mean_Z_P0.05 = mean(Zscore_P0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Method %in% method_order) %>%
  column_to_rownames("Method") %>%
  .[method_order, ]  # 按 cluster 顺序排序

# 转换矩阵并替换 Inf/NA
mat <- as.matrix(mean_z)
mat[!is.finite(mat)] <- NA
mat[is.na(mat)] <- max(mat, na.rm = TRUE)

# 绘图
pdf("module_construction_heatmap_meanZscore_byCluster.pdf", width = 6, height = 10)
pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white","#F8D6A1","orange", "red"))(100),
  main = "Mean Jaccard Z-score per Method (Clustered)",
  angle_col = 45
)
dev.off()
