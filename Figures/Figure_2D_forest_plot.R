
#2025-09-04


library(tidyverse)
library(ggplot2)
library(scales)

df <- read.csv("Figure_2_monocytes_forest_GWAS_SNP_RBP_all_SNPs.csv", check.names = FALSE, stringsAsFactors = FALSE)
df <- read.csv("Figure_2_monocytes_forest_GWAS_SNP_RBP_keep_topOneSNP.csv", check.names = FALSE, stringsAsFactors = FALSE)
df$sample_size <- as.numeric(gsub(",", "", df$sample_size))
df$p_value <- pmax(as.numeric(df$p_value), .Machine$double.xmin)  # 防下溢
df$region    <- factor(df$region, levels = c("5UTR","3UTR","Introns"))
df$Cell_type <- factor(df$Cell_type, levels = c("Monocytes","T_cells"))


# ---- 假设 df 已经读入（含 RBP, SNP, beta, p_value, sample_size, Cell_type, region）----

# 计算区间 & 排序：每个 Cell_type 内按 beta 升序（更负在底部）
df1 <- df %>%
  mutate(
    sample_size = as.numeric(gsub(",", "", sample_size)),
    p_value = pmax(as.numeric(p_value), .Machine$double.xmin),
    z   = qnorm(1 - p_value/2),
    se  = abs(beta) / z,
    lo95 = beta - 1.96*se,
    hi95 = beta + 1.96*se,
    mlogp = -log10(p_value),
    label = paste0(RBP, " (", SNP, ")")
  ) %>%
  group_by(Cell_type) %>%
  arrange(beta, .by_group = TRUE) %>%                     # 按 beta 排序
  mutate(y = factor(label, levels = unique(label))) %>%   # 排序写进 y 轴
  ungroup()

# 颜色
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)

# x 轴留白
xl <- quantile(df1$lo95, 0.02, na.rm = TRUE)
xr <- quantile(df1$hi95, 0.98, na.rm = TRUE)
rng <- xr - xl
xlim_use <- c(xl - 0.05*rng, xr + 0.05*rng)

p <- ggplot(df1, aes(x = beta, y = y)) +
  facet_grid(Cell_type ~ ., scales = "free_y", space = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo95, xmax = hi95),
                 height = 0.22, color = "grey60", linewidth = 0.6) +
  # —— 不要 region shape：统一用 shape=21，颜色映射只保留 -log10(P)，大小映射样本量 ——
  geom_point(aes(fill = mlogp, size = sample_size),
             shape = 21, color = "grey30", stroke = 0.3) +
  scale_fill_gradientn(
    colours = pal, name = "-log10(P)",
    guide = guide_colorbar(nbin = 256, raster = TRUE,
                           barheight = unit(35, "mm"), barwidth = unit(4, "mm"))
  ) +
  scale_size_continuous(name = "Sample size", range = c(1.6, 3.6), labels = comma) +
  coord_cartesian(xlim = xlim_use, clip = "off") +
  labs(x = "Effect size (Beta)", y = NULL,
       title = "RBP regulons associated with monocyte count") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9, color = "grey20"),
    legend.position = "right",
    plot.margin = margin(5, 10, 5, 5)
  )

# 高度跟行数走，避免拥挤
h <- max(6, 0.28 * nrow(df1))
ggsave("Figure_2_monocytes_forest_GWAS_SNP_RBP_keep_topOneSNP.pdf", p, width = 9, height = h)
p
