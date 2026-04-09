# 1. 解析命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11) {
  stop("Usage: Rscript script.R <real1> <fake1> ... <real5> <fake5> <metric> <output_plot>")
}

real_files <- args[1:5]
fake_files <- args[6:10]
metric <- args[11]
output_plot <- args[12]

# 检查 metric 合法性
if (!(metric %in% c("scRBP_Jaccard", "Z_Scaled_Jaccard"))) {
  stop("Metric must be 'scRBP_Jaccard' or 'Z_Scaled_Jaccard'")
}

# 2. 读取并处理 5 组数据
library(dplyr)
combined_list <- list()

method_names <- c("Method1", "Method2", "Method3", "Method4", "Method5")

for (i in 1:5) {
  real_df <- read.csv(real_files[i])
  fake_df <- read.csv(fake_files[i])
  
  real_df$Group <- "Real"
  fake_df$Group <- "Fake"
  
  # 随机抽取 fake 样本
  set.seed(123)
  fake_sample <- fake_df[sample(nrow(fake_df), nrow(real_df)), ]
  
  # 合并，并添加 Method 列
  combined <- rbind(real_df, fake_sample) %>%
    mutate(Method = method_names[i])
  
  combined_list[[i]] <- combined
}

# 3. 汇总所有方法
final_df <- bind_rows(combined_list)

# 4. 绘图
library(ggplot2)

p <- ggplot(final_df, aes(x = Group, y = .data[[metric]], fill = Group)) +
  geom_boxplot() +
  facet_wrap(~Method, scales = "free_y") +
  theme_minimal() +
  labs(title = paste("Comparison of", metric, "across 5 Methods"), y = metric) +
  scale_fill_manual(values = c("Real" = "skyblue", "Fake" = "orange"))

ggsave(output_plot, plot = p, width = 10, height = 6)

