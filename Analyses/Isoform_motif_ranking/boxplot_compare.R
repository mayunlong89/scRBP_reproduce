
# 1. 解析命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript script.R <real_file> <fake_file> <metric> <output_plot> <plot_title>")
}

real_file <- args[1]
fake_file <- args[2]
metric <- args[3]  # "scRBP_Jaccard" or "Z_Scaled_Jaccard"
output_plot <- args[4]
plot_title <- args[5]

# 2. 读取数据
real_df <- read.csv(real_file)
fake_df <- read.csv(fake_file)

# 检查 metric 合法性
if (!(metric %in% c("scRBP_Jaccard", "Z_Scaled_Jaccard"))) {
  stop("Metric must be 'scRBP_Jaccard' or 'Z_Scaled_Jaccard'")
}

# 3. 添加标签
real_df$Group <- "Real"
fake_df$Group <- "Fake"

# 4. 从 fake_df 随机抽取与 real_df 相同行数
set.seed(123)
#fake_sample <- fake_df[sample(nrow(fake_df), nrow(real_df)), ]
fake_sample <- fake_df
# 5. 合并
combined_df <- rbind(real_df, fake_sample)

# 6. 绘图
library(ggplot2)

p <- ggplot(combined_df, aes(x = Group, y = .data[[metric]], fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = plot_title, y = metric) +
  scale_fill_manual(values = c("Real" = "skyblue", "Fake" = "orange"))

ggsave(output_plot, plot = p, width = 6, height = 4)




