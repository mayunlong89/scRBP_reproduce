# 1. 解析命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript script.R <real_file> <fake_file> <metric> <output_plot> <output_pval>")
}

real_file <- args[1]
fake_file <- args[2]
metric <- args[3]  # "scRBP_Jaccard" or "Z_Scaled_Jaccard"
output_plot <- args[4]
output_pval <- args[5]

# 2. 读取数据
real_df <- read.csv(real_file)
fake_df <- read.csv(fake_file)

# 检查 metric 合法性
if (!(metric %in% c("scRBP_Jaccard", "Z_Scaled_Jaccard"))) {
  stop("Metric must be 'scRBP_Jaccard' or 'Z_Scaled_Jaccard'")
}

# 3. 计算真实均值
real_mean <- mean(real_df[[metric]])

# 4. Bootstrap
set.seed(123)
bootstrap_means <- replicate(1000, {
  fake_sample <- fake_df[sample(nrow(fake_df), nrow(real_df)), ]
  mean(fake_sample[[metric]])
})

# 5. 绘图
library(ggplot2)

boot_df <- data.frame(Mean = bootstrap_means)

p <- ggplot(boot_df, aes(x = Mean)) +
  geom_histogram(bins = 30, fill = "grey", color = "black") +
  geom_vline(xintercept = real_mean, color = "red", linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = paste("Bootstrap Null Distribution of Fake Means (", metric, ")", sep = ""),
       x = paste("Mean", metric), y = "Frequency",
       subtitle = paste("Red dashed line = Real Mean =", round(real_mean, 3)))

ggsave(output_plot, plot = p, width = 6, height = 4)

# 6. 计算 p 值 (双尾)
p_value <- mean(abs(bootstrap_means) >= abs(real_mean))
write(paste("P-value:", p_value), file = output_pval)

