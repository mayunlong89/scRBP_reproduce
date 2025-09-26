
#----2025-09-08

library(tidyverse)
library(stringr)
library(scales)

# 1) 读入你的汇总表（ID, cell_type, RBP, TRS_6.5K, TRS_8.7K）
df <- read.csv("02_Summary_monocyte_count_correlation_TRS.csv", check.names = FALSE, stringsAsFactors = FALSE)

# 2) 从 ID 解析 region（末尾的 5UTR/3UTR/CDS/Introns）
df <- df %>%
  mutate(
    region = str_extract(ID, "(5UTR|3UTR|CDS|Introns?)$"),
    region = recode(region, "Intron" = "Introns"),
    region = factor(region, levels = c("5UTR","3UTR","CDS","Introns")),
    cell_type = factor(cell_type, levels = c("Monocytes","T_cells"))
  )

# 3) 计算各面板的 Spearman 相关
spearman_p <- function(x, y) suppressWarnings(cor.test(x, y, method = "spearman")$p.value)

stats <- df %>%
  group_by(cell_type, region) %>%
  summarise(
    n    = n(),
    rho  = suppressWarnings(cor(`TRS_6.5K`, `TRS_8.7K`, method = "spearman",
                                use = "pairwise.complete.obs")),
    pval = spearman_p(`TRS_6.5K`, `TRS_8.7K`),
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(pval, "BH"),
         label = sprintf("ρ = %.2f\nn = %d\nFDR = %.3g", rho, n, FDR))

# 4) Facet 散点相关图
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)  # 你的配色

p_scatter <- ggplot(df, aes(`TRS_6.5K`, `TRS_8.7K`)) +
  geom_point(alpha = 0.8, size = 1.8, color = "#6b6b6b", fill = "#bdbdbd", shape = 21, stroke = 0.2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = "grey40") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
  facet_grid(cell_type ~ region, scales = "free") +
  # 在每个面板左上角标注 ρ / n / FDR
  geom_text(data = stats, inherit.aes = FALSE,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.1, size = 3.2, color = "grey20") +
  labs(
    title = "TRS correlation between 6.5K and 8.7K",
    x = "TRS (6.5K)", y = "TRS (8.7K)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0)
  )

# 保存散点图（高度按面板数自动）
ggsave("TRS_correlation_scatter.pdf", p_scatter, width = 10, height = 6)

p_scatter



##--------Figure--2---Figure 2D---correlation

library(tidyverse)
library(stringr)
library(scales)

# 1) 读取并整理（列需含：ID, cell_type, RBP, TRS_6.5K, TRS_8.7K）
df <- read.csv("02_Summary_monocyte_count_correlation_TRS.csv", check.names = FALSE, stringsAsFactors = FALSE)

cor.test(df$TRS_6.5K,df$TRS_8.7K)


# 从 ID 提取 region
df <- df %>%
  mutate(
    region    = str_extract(ID, "(5UTR|3UTR|CDS|Introns?)$"),
    region    = recode(region, "Intron" = "Introns"),
    region    = factor(region, levels = c("5UTR","3UTR","CDS","Introns")),
    cell_type = factor(cell_type, levels = c("Monocytes","T_cells"))
  ) %>%
  drop_na(`TRS_6.5K`, `TRS_8.7K`)  # 仅保留成对存在的数据

# 2) 计算每个 cell type 的相关，用于图内标注
spearman_p <- function(x, y) suppressWarnings(cor.test(x, y, method = "spearman")$p.value)
stats_ct <- df %>%
  group_by(cell_type) %>%
  summarise(
    n    = n(),
    rho  = suppressWarnings(cor(`TRS_6.5K`, `TRS_8.7K`, method = "spearman")),
    pval = spearman_p(`TRS_6.5K`, `TRS_8.7K`),
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(pval, "BH"),
         label = sprintf("%s: ρ=%.2f, n=%d, FDR=%.3g",
                         as.character(cell_type), rho, n, FDR))

# 3) 颜色/形状（与 UMAP 的红/青两色保持一致，可按需微调）
cols_ct <- c(Monocytes = "#A84B4B", T_cells = "#53A6A3")
shapes  <- c("5UTR"=21, "3UTR"=22, "CDS"=24, "Introns"=23)  # 21/22/23/24 允许填充

# 4) 合成主图：颜色=cell type，形状=region；每个 cell type 一条拟合线
p_all <- ggplot(df, aes(`TRS_6.5K`, `TRS_8.7K`)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(fill = cell_type, shape = region),
             size = 2.2, alpha = 0.9, color = "grey25", stroke = 0.35) +
  stat_smooth(aes(color = cell_type), method = "lm", se = FALSE, linewidth = 0.9) +
  scale_fill_manual(values = cols_ct, name = "Cell type") +
  scale_color_manual(values = cols_ct, guide = "none") +
  scale_shape_manual(values = shapes, name = "mRNA region") +
  coord_equal() +
  labs(
    title = "TRS correlation between 6.5K and 8.7K",
    subtitle = "Color = cell type; shape = mRNA region; dashed = y = x",
    x = "TRS (6.5K)", y = "TRS (8.7K)"
  ) +
  # 在图内标注两个 cell type 的 ρ 与 n
  geom_label(
    data = stats_ct,
    aes(x = -Inf, y = Inf, label = label, fill = cell_type),
    inherit.aes = FALSE,
    hjust = -0.02, vjust = 1.1, color = "white",
    alpha = 0.9, label.size = 0
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(hjust = 0),
    legend.position = "right"
  )

ggsave("TRS_correlation_combined_main.pdf", p_all, width = 8.5, height = 6)
p_all



library(ggplot2)

# 颜色/形状（与之前保持一致）
cols_ct <- c(Monocytes = "#A84B4B", T_cells = "#53A6A3")
shapes  <- c("5UTR"=21, "3UTR"=22, "CDS"=24, "Introns"=23)

# 全局相关用于角标（可选）
rho_all <- suppressWarnings(cor(df$`TRS_6.5K`, df$`TRS_8.7K`, method = "spearman"))
p_all   <- suppressWarnings(cor.test(df$`TRS_6.5K`, df$`TRS_8.7K`, method = "spearman")$p.value)
lab_all <- sprintf("ρ = %.2f, n = %d, FDR = %.3g", rho_all, nrow(df), p.adjust(p_all, "BH"))

p <- ggplot(df, aes(`TRS_6.5K`, `TRS_8.7K`)) +
  # 散点：颜色=cell type，形状=region
  geom_point(aes(fill = cell_type, shape = region),
             size = 2.2, alpha = 0.9, color = "grey25", stroke = 0.35) +
  # **唯一一条全局直线虚线（线性拟合）**
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              linetype = "dashed", color = "black", linewidth = 1) +
  scale_fill_manual(values = cols_ct, name = "Cell type") +
  scale_shape_manual(values = shapes,  name = "mRNA region") +
  coord_equal() +
  labs(
    title = "TRS correlation between 6.5K and 8.7K",
    subtitle = "Dashed line = global linear fit (all points)",
    x = "TRS (6.5K)", y = "TRS (8.7K)"
  ) +
  # 角标（可选）
  annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
           label = lab_all, size = 3.2, label.size = 0, fill = "grey95") +
  theme_classic(base_size = 12)

ggsave("TRS_correlation_one_dashed_line.pdf", p, width = 8.5, height = 6)
p

cols_ct <- c(Monocytes = "#A84B4B", T_cells = "#53A6A3")
shapes  <- c("5UTR"=21, "3UTR"=22, "CDS"=24, "Introns"=23)  # 21/22/23/24 才会显示 fill

p <- ggplot(df, aes(`TRS_6.5K`, `TRS_8.7K`)) +
  geom_point(aes(fill = cell_type, shape = region),
             size = 2.2, alpha = 0.9, color = "grey25", stroke = 0.35) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              linetype = "dashed", color = "black", linewidth = 1) +
  scale_shape_manual(values = shapes, name = "mRNA region") +
  scale_fill_manual(values = cols_ct, name = "Cell type") +
  guides(
    fill = guide_legend(                     # 让图例按“填充色”显示
      override.aes = list(shape = 21,        # 实心圆作为图例样式
                          size  = 4,
                          color = "grey25",  # 与点的描边一致
                          alpha = 1)
    )
  ) +
  coord_equal() +
  labs(title = "TRS correlation between 6.5K and 8.7K",
       subtitle = "Dashed line = global linear fit (all points)",
       x = "TRS (6.5K)", y = "TRS (8.7K)") +
  theme_classic(base_size = 12)
p



ct   <- cor.test(df$TRS_6.5K, df$TRS_8.7K, method = "spearman")
rho  <- unname(ct$estimate)
npts <- sum(complete.cases(df$TRS_6.5K, df$TRS_8.7K))

p_txt <- format.pval(ct$p.value, digits = 3, eps = 1e-3000)  # 关键：避免显示为 0

lab_txt <- sprintf("rho = %.2f\nn = %d\nP = %s", rho, npts, p_txt)

p + annotate("label",
             x = min(df$TRS_6.5K, na.rm = TRUE) + 5,
             y = max(df$TRS_8.7K, na.rm = TRUE) - 5,
             label = lab_txt,
             size = 4, label.size = 0, fill = "grey95")



p <- ggplot(df, aes(`TRS_6.5K`, `TRS_8.7K`)) +
  geom_point(aes(fill = cell_type, shape = region),
             size = 2.2, alpha = 0.9, color = "grey25", stroke = 0.35) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              linetype = "dashed", color = "black", linewidth = 1) +
  ## 参考虚线：TRS = 3
  geom_vline(xintercept = 3, linetype = "dashed", color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "steelblue", linewidth = 1) +
  scale_shape_manual(values = shapes, name = "mRNA region") +
  scale_fill_manual(values = cols_ct, name = "Cell type") +
  guides(fill = guide_legend(override.aes = list(
    shape = 21, size = 4, color = "grey25", alpha = 1))) +
  coord_equal() +
  labs(title = "TRS correlation between 6.5K and 8.7K",
       subtitle = "Dashed line = global linear fit (all points)",
       x = "TRS (6.5K)", y = "TRS (8.7K)") +
  theme_classic(base_size = 12)

p


library(ggplot2)
library(scales)

# 把 cut 左侧按 factor 倍压缩（保持单调，适合本图）
compress_left_trans <- function(cut = -25, factor = 4){
  trans_new(
    paste0("compress_left_", cut, "_", factor),
    transform = function(x) ifelse(is.na(x), NA_real_,
                                   ifelse(x < cut, cut + (x - cut)/factor, x)),
    inverse   = function(x) ifelse(is.na(x), NA_real_,
                                   ifelse(x < cut, cut + (x - cut)*factor, x))
  )
}

p2 <- p +
  coord_trans(
    x = compress_left_trans(cut = -25, factor = 4),
    y = compress_left_trans(cut = -25, factor = 4)
  ) +
  scale_x_continuous(breaks = c(-80,-60,-40,-25,-10,-5,0,5,10), minor_breaks = NULL) +
  scale_y_continuous(breaks = c(-80,-60,-40,-25,-10,-5,0,5,10), minor_breaks = NULL)
# 如需 TRS=3 的参考虚线:
# + geom_vline(xintercept = 3, linetype = "22", color = "grey40") +
#   geom_hline(yintercept = 3, linetype = "22", color = "grey40")

p2


library(scales)

tfun <- scales::pseudo_log_trans(sigma = 12)  # 可调：8~20，越大越“线性”
p2 <- p +
  coord_trans(x = tfun, y = tfun) +
  scale_x_continuous(
    breaks = c(-80, -60, -40, -25, -10, -5, 0, 5),
    labels = label_number()
  ) +
  scale_y_continuous(
    breaks = c(-80, -60, -40, -25, -10, -5, 0, 5),
    labels = label_number()
  )
p2





##--------Figure--2---Supplementary_Figure 10B-C---correlation--Lymphocyte count

library(tidyverse)
library(stringr)
library(scales)

# 1) 读取并整理（列需含：ID, cell_type, RBP, TRS_6.5K, TRS_8.7K）
df <- read.csv("02_Summary_monocyte_lymphocyte_count_correlation_TRS.csv", check.names = FALSE, stringsAsFactors = FALSE)

cor.test(df$TRS_6.5K_lymp_count,df$TRS_8.7K_lymp_count)


# 从 ID 提取 region
df <- df %>%
  mutate(
    region    = str_extract(ID, "(5UTR|3UTR|CDS|Introns?)$"),
    region    = recode(region, "Intron" = "Introns"),
    region    = factor(region, levels = c("5UTR","3UTR","CDS","Introns")),
    cell_type = factor(cell_type, levels = c("Monocytes","T_cells"))
  ) %>%
  drop_na(`TRS_6.5K_lymp_count`, `TRS_8.7K_lymp_count`)  # 仅保留成对存在的数据

head(df)


# 2) 计算每个 cell type 的相关，用于图内标注
spearman_p <- function(x, y) suppressWarnings(cor.test(x, y, method = "spearman")$p.value)
stats_ct <- df %>%
  group_by(cell_type) %>%
  summarise(
    n    = n(),
    rho  = suppressWarnings(cor(`TRS_6.5K_lymp_count`, `TRS_8.7K_lymp_count`, method = "spearman")),
    pval = spearman_p(`TRS_6.5K_lymp_count`, `TRS_8.7K_lymp_count`),
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(pval, "BH"),
         label = sprintf("%s: ρ=%.2f, n=%d, FDR=%.3g",
                         as.character(cell_type), rho, n, FDR))

# 3) 颜色/形状（与 UMAP 的红/青两色保持一致，可按需微调）
cols_ct <- c(Monocytes = "#A84B4B", T_cells = "#53A6A3")
shapes  <- c("5UTR"=21, "3UTR"=22, "CDS"=24, "Introns"=23)  # 21/22/23/24 允许填充

# 4) 合成主图：颜色=cell type，形状=region；每个 cell type 一条拟合线
p_all <- ggplot(df, aes(`TRS_6.5K_lymp_count`, `TRS_8.7K_lymp_count`)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(fill = cell_type, shape = region),
             size = 2.2, alpha = 0.9, color = "grey25", stroke = 0.35) +
  stat_smooth(aes(color = cell_type), method = "lm", se = FALSE, linewidth = 0.9) +
  scale_fill_manual(values = cols_ct, name = "Cell type") +
  scale_color_manual(values = cols_ct, guide = "none") +
  scale_shape_manual(values = shapes, name = "mRNA region") +
  coord_equal() +
  labs(
    title = "TRS correlation between 6.5K_lymp_count and 8.7K_lymp_count",
    subtitle = "Color = cell type; shape = mRNA region; dashed = y = x",
    x = "TRS (6.5K)", y = "TRS (8.7K)"
  ) +
  # 在图内标注两个 cell type 的 ρ 与 n
  geom_label(
    data = stats_ct,
    aes(x = -Inf, y = Inf, label = label, fill = cell_type),
    inherit.aes = FALSE,
    hjust = -0.02, vjust = 1.1, color = "white",
    alpha = 0.9, label.size = 0
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(hjust = 0),
    legend.position = "right"
  )

ggsave("TRS_correlation_combined_main_lymp_count.pdf", p_all, width = 8.5, height = 6)
p_all



library(ggplot2)

# 颜色/形状（与之前保持一致）
cols_ct <- c(Monocytes = "#A84B4B", T_cells = "#53A6A3")
shapes  <- c("5UTR"=21, "3UTR"=22, "CDS"=24, "Introns"=23)

# 全局相关用于角标（可选）
rho_all <- suppressWarnings(cor(df$`TRS_6.5K_lymp_count`, df$`TRS_8.7K_lymp_count`, method = "spearman"))
p_all   <- suppressWarnings(cor.test(df$`TRS_6.5K_lymp_count`, df$`TRS_8.7K_lymp_count`, method = "spearman")$p.value)
lab_all <- sprintf("ρ = %.2f, n = %d, FDR = %.3g", rho_all, nrow(df), p.adjust(p_all, "BH"))

p <- ggplot(df, aes(`TRS_6.5K_lymp_count`, `TRS_8.7K_lymp_count`)) +
  # 散点：颜色=cell type，形状=region
  geom_point(aes(fill = cell_type, shape = region),
             size = 2.2, alpha = 0.9, color = "grey25", stroke = 0.35) +
  # **唯一一条全局直线虚线（线性拟合）**
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              linetype = "dashed", color = "black", linewidth = 1) +
  scale_fill_manual(values = cols_ct, name = "Cell type") +
  scale_shape_manual(values = shapes,  name = "mRNA region") +
  coord_equal() +
  labs(
    title = "TRS correlation between 6.5K_lymp_count and 8.7K_lymp_count",
    subtitle = "Dashed line = global linear fit (all points)",
    x = "TRS (6.5K)", y = "TRS (8.7K)"
  ) +
  # 角标（可选）
  annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
           label = lab_all, size = 3.2, label.size = 0, fill = "grey95") +
  theme_classic(base_size = 12)

ggsave("TRS_correlation_one_dashed_line_lymp_count.pdf", p, width = 8.5, height = 6)
p

cols_ct <- c(Monocytes = "#A84B4B", T_cells = "#53A6A3")
shapes  <- c("5UTR"=21, "3UTR"=22, "CDS"=24, "Introns"=23)  # 21/22/23/24 才会显示 fill

p <- ggplot(df, aes(`TRS_6.5K_lymp_count`, `TRS_8.7K_lymp_count`)) +
  geom_point(aes(fill = cell_type, shape = region),
             size = 2.2, alpha = 0.9, color = "grey25", stroke = 0.35) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              linetype = "dashed", color = "black", linewidth = 1) +
  scale_shape_manual(values = shapes, name = "mRNA region") +
  scale_fill_manual(values = cols_ct, name = "Cell type") +
  guides(
    fill = guide_legend(                     # 让图例按“填充色”显示
      override.aes = list(shape = 21,        # 实心圆作为图例样式
                          size  = 4,
                          color = "grey25",  # 与点的描边一致
                          alpha = 1)
    )
  ) +
  coord_equal() +
  labs(title = "TRS correlation between 6.5K_lymp_count and 8.7K_lymp_count",
       subtitle = "Dashed line = global linear fit (all points)",
       x = "TRS (6.5K)", y = "TRS (8.7K)") +
  theme_classic(base_size = 12)
p



ct   <- cor.test(df$TRS_6.5K_lymp_count, df$TRS_8.7K_lymp_count, method = "spearman")
rho  <- unname(ct$estimate)
npts <- sum(complete.cases(df$TRS_6.5K_lymp_count, df$TRS_8.7K_lymp_count))

p_txt <- format.pval(ct$p.value, digits = 3, eps = 1e-3000)  # 关键：避免显示为 0

lab_txt <- sprintf("rho = %.2f\nn = %d\nP = %s", rho, npts, p_txt)

p + annotate("label",
             x = min(df$TRS_6.5K_lymp_count, na.rm = TRUE) + 5,
             y = max(df$TRS_8.7K_lymp_count, na.rm = TRUE) - 5,
             label = lab_txt,
             size = 4, label.size = 0, fill = "grey95")



p <- ggplot(df, aes(`TRS_6.5K_lymp_count`, `TRS_8.7K_lymp_count`)) +
  geom_point(aes(fill = cell_type, shape = region),
             size = 2.2, alpha = 0.9, color = "grey25", stroke = 0.35) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              linetype = "dashed", color = "black", linewidth = 1) +
  ## 参考虚线：TRS = 3
  geom_vline(xintercept = 3, linetype = "dashed", color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "steelblue", linewidth = 1) +
  scale_shape_manual(values = shapes, name = "mRNA region") +
  scale_fill_manual(values = cols_ct, name = "Cell type") +
  guides(fill = guide_legend(override.aes = list(
    shape = 21, size = 4, color = "grey25", alpha = 1))) +
  coord_equal() +
  labs(title = "TRS correlation between 6.5K_lymp_count and 8.7K_lymp_count",
       subtitle = "Dashed line = global linear fit (all points)",
       x = "TRS (6.5K)", y = "TRS (8.7K)") +
  theme_classic(base_size = 12)

p


library(ggplot2)
library(scales)

# 把 cut 左侧按 factor 倍压缩（保持单调，适合本图）
compress_left_trans <- function(cut = -25, factor = 4){
  trans_new(
    paste0("compress_left_", cut, "_", factor),
    transform = function(x) ifelse(is.na(x), NA_real_,
                                   ifelse(x < cut, cut + (x - cut)/factor, x)),
    inverse   = function(x) ifelse(is.na(x), NA_real_,
                                   ifelse(x < cut, cut + (x - cut)*factor, x))
  )
}

p2 <- p +
  coord_trans(
    x = compress_left_trans(cut = -25, factor = 4),
    y = compress_left_trans(cut = -25, factor = 4)
  ) +
  scale_x_continuous(breaks = c(-80,-60,-40,-25,-10,-5,0,5,10), minor_breaks = NULL) +
  scale_y_continuous(breaks = c(-80,-60,-40,-25,-10,-5,0,5,10), minor_breaks = NULL)
# 如需 TRS=3 的参考虚线:
# + geom_vline(xintercept = 3, linetype = "22", color = "grey40") +
#   geom_hline(yintercept = 3, linetype = "22", color = "grey40")

p2


library(scales)

tfun <- scales::pseudo_log_trans(sigma = 12)  # 可调：8~20，越大越“线性”
p2 <- p +
  coord_trans(x = tfun, y = tfun) +
  scale_x_continuous(
    breaks = c(-80, -60, -40, -25, -10, -5, 0, 5),
    labels = label_number()
  ) +
  scale_y_continuous(
    breaks = c(-80, -60, -40, -25, -10, -5, 0, 5),
    labels = label_number()
  )
p2





##--------Figure--2---Supplementary_Figure 10B-C---correlation--Lymphocyte percent 

library(tidyverse)
library(stringr)
library(scales)

# 1) 读取并整理（列需含：ID, cell_type, RBP, TRS_6.5K, TRS_8.7K）
df <- read.csv("02_Summary_monocyte_lymphocyte_count_correlation_TRS.csv", check.names = FALSE, stringsAsFactors = FALSE)

cor.test(df$TRS_6.5K_lymp_percent,df$TRS_8.7K_lymp_percent)


# 从 ID 提取 region
df <- df %>%
  mutate(
    region    = str_extract(ID, "(5UTR|3UTR|CDS|Introns?)$"),
    region    = recode(region, "Intron" = "Introns"),
    region    = factor(region, levels = c("5UTR","3UTR","CDS","Introns")),
    cell_type = factor(cell_type, levels = c("Monocytes","T_cells"))
  ) %>%
  drop_na(`TRS_6.5K_lymp_percent`, `TRS_8.7K_lymp_percent`)  # 仅保留成对存在的数据

head(df)


# 2) 计算每个 cell type 的相关，用于图内标注
spearman_p <- function(x, y) suppressWarnings(cor.test(x, y, method = "spearman")$p.value)
stats_ct <- df %>%
  group_by(cell_type) %>%
  summarise(
    n    = n(),
    rho  = suppressWarnings(cor(`TRS_6.5K_lymp_percent`, `TRS_8.7K_lymp_percent`, method = "spearman")),
    pval = spearman_p(`TRS_6.5K_lymp_percent`, `TRS_8.7K_lymp_percent`),
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(pval, "BH"),
         label = sprintf("%s: ρ=%.2f, n=%d, FDR=%.3g",
                         as.character(cell_type), rho, n, FDR))

# 3) 颜色/形状（与 UMAP 的红/青两色保持一致，可按需微调）
cols_ct <- c(Monocytes = "#A84B4B", T_cells = "#53A6A3")
shapes  <- c("5UTR"=21, "3UTR"=22, "CDS"=24, "Introns"=23)  # 21/22/23/24 允许填充

# 4) 合成主图：颜色=cell type，形状=region；每个 cell type 一条拟合线
p_all <- ggplot(df, aes(`TRS_6.5K_lymp_percent`, `TRS_8.7K_lymp_percent`)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(fill = cell_type, shape = region),
             size = 2.2, alpha = 0.9, color = "grey25", stroke = 0.35) +
  stat_smooth(aes(color = cell_type), method = "lm", se = FALSE, linewidth = 0.9) +
  scale_fill_manual(values = cols_ct, name = "Cell type") +
  scale_color_manual(values = cols_ct, guide = "none") +
  scale_shape_manual(values = shapes, name = "mRNA region") +
  coord_equal() +
  labs(
    title = "TRS correlation between 6.5K_lymp_percent and 8.7K_lymp_percent",
    subtitle = "Color = cell type; shape = mRNA region; dashed = y = x",
    x = "TRS (6.5K)", y = "TRS (8.7K)"
  ) +
  # 在图内标注两个 cell type 的 ρ 与 n
  geom_label(
    data = stats_ct,
    aes(x = -Inf, y = Inf, label = label, fill = cell_type),
    inherit.aes = FALSE,
    hjust = -0.02, vjust = 1.1, color = "white",
    alpha = 0.9, label.size = 0
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(hjust = 0),
    legend.position = "right"
  )

ggsave("TRS_correlation_combined_main_lymp_percent.pdf", p_all, width = 8.5, height = 6)
p_all



library(ggplot2)

# 颜色/形状（与之前保持一致）
cols_ct <- c(Monocytes = "#A84B4B", T_cells = "#53A6A3")
shapes  <- c("5UTR"=21, "3UTR"=22, "CDS"=24, "Introns"=23)

# 全局相关用于角标（可选）
rho_all <- suppressWarnings(cor(df$`TRS_6.5K_lymp_percent`, df$`TRS_8.7K_lymp_percent`, method = "spearman"))
p_all   <- suppressWarnings(cor.test(df$`TRS_6.5K_lymp_percent`, df$`TRS_8.7K_lymp_percent`, method = "spearman")$p.value)
lab_all <- sprintf("ρ = %.2f, n = %d, FDR = %.3g", rho_all, nrow(df), p.adjust(p_all, "BH"))

p <- ggplot(df, aes(`TRS_6.5K_lymp_percent`, `TRS_8.7K_lymp_percent`)) +
  # 散点：颜色=cell type，形状=region
  geom_point(aes(fill = cell_type, shape = region),
             size = 2.2, alpha = 0.9, color = "grey25", stroke = 0.35) +
  # **唯一一条全局直线虚线（线性拟合）**
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              linetype = "dashed", color = "black", linewidth = 1) +
  scale_fill_manual(values = cols_ct, name = "Cell type") +
  scale_shape_manual(values = shapes,  name = "mRNA region") +
  coord_equal() +
  labs(
    title = "TRS correlation between 6.5K_lymp_percent and 8.7K_lymp_percent",
    subtitle = "Dashed line = global linear fit (all points)",
    x = "TRS (6.5K)", y = "TRS (8.7K)"
  ) +
  # 角标（可选）
  annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
           label = lab_all, size = 3.2, label.size = 0, fill = "grey95") +
  theme_classic(base_size = 12)

ggsave("TRS_correlation_one_dashed_line_lymp_percent.pdf", p, width = 8.5, height = 6)
p

cols_ct <- c(Monocytes = "#A84B4B", T_cells = "#53A6A3")
shapes  <- c("5UTR"=21, "3UTR"=22, "CDS"=24, "Introns"=23)  # 21/22/23/24 才会显示 fill

p <- ggplot(df, aes(`TRS_6.5K_lymp_count`, `TRS_8.7K_lymp_count`)) +
  geom_point(aes(fill = cell_type, shape = region),
             size = 2.2, alpha = 0.9, color = "grey25", stroke = 0.35) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              linetype = "dashed", color = "black", linewidth = 1) +
  scale_shape_manual(values = shapes, name = "mRNA region") +
  scale_fill_manual(values = cols_ct, name = "Cell type") +
  guides(
    fill = guide_legend(                     # 让图例按“填充色”显示
      override.aes = list(shape = 21,        # 实心圆作为图例样式
                          size  = 4,
                          color = "grey25",  # 与点的描边一致
                          alpha = 1)
    )
  ) +
  coord_equal() +
  labs(title = "TRS correlation between 6.5K_lymp_percent and 8.7K_lymp_percent",
       subtitle = "Dashed line = global linear fit (all points)",
       x = "TRS (6.5K)", y = "TRS (8.7K)") +
  theme_classic(base_size = 12)
p



ct   <- cor.test(df$TRS_6.5K_lymp_percent, df$TRS_8.7K_lymp_percent, method = "spearman")
rho  <- unname(ct$estimate)
npts <- sum(complete.cases(df$TRS_6.5K_lymp_percent, df$TRS_8.7K_lymp_percent))

p_txt <- format.pval(ct$p.value, digits = 3, eps = 1e-3000)  # 关键：避免显示为 0

lab_txt <- sprintf("rho = %.2f\nn = %d\nP = %s", rho, npts, p_txt)

p + annotate("label",
             x = min(df$TRS_6.5K_lymp_percent, na.rm = TRUE) + 5,
             y = max(df$TRS_8.7K_lymp_percent, na.rm = TRUE) - 5,
             label = lab_txt,
             size = 4, label.size = 0, fill = "grey95")



p <- ggplot(df, aes(`TRS_6.5K_lymp_percent`, `TRS_8.7K_lymp_percent`)) +
  geom_point(aes(fill = cell_type, shape = region),
             size = 2.2, alpha = 0.9, color = "grey25", stroke = 0.35) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              linetype = "dashed", color = "black", linewidth = 1) +
  ## 参考虚线：TRS = 3
  geom_vline(xintercept = 3, linetype = "dashed", color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "steelblue", linewidth = 1) +
  scale_shape_manual(values = shapes, name = "mRNA region") +
  scale_fill_manual(values = cols_ct, name = "Cell type") +
  guides(fill = guide_legend(override.aes = list(
    shape = 21, size = 4, color = "grey25", alpha = 1))) +
  coord_equal() +
  labs(title = "TRS correlation between 6.5K_lymp_percent and 8.7K_lymp_percent",
       subtitle = "Dashed line = global linear fit (all points)",
       x = "TRS (6.5K)", y = "TRS (8.7K)") +
  theme_classic(base_size = 12)

p


library(ggplot2)
library(scales)

# 把 cut 左侧按 factor 倍压缩（保持单调，适合本图）
compress_left_trans <- function(cut = -25, factor = 4){
  trans_new(
    paste0("compress_left_", cut, "_", factor),
    transform = function(x) ifelse(is.na(x), NA_real_,
                                   ifelse(x < cut, cut + (x - cut)/factor, x)),
    inverse   = function(x) ifelse(is.na(x), NA_real_,
                                   ifelse(x < cut, cut + (x - cut)*factor, x))
  )
}

p2 <- p +
  coord_trans(
    x = compress_left_trans(cut = -25, factor = 4),
    y = compress_left_trans(cut = -25, factor = 4)
  ) +
  scale_x_continuous(breaks = c(-80,-60,-40,-25,-10,-5,0,5,10), minor_breaks = NULL) +
  scale_y_continuous(breaks = c(-80,-60,-40,-25,-10,-5,0,5,10), minor_breaks = NULL)
# 如需 TRS=3 的参考虚线:
# + geom_vline(xintercept = 3, linetype = "22", color = "grey40") +
#   geom_hline(yintercept = 3, linetype = "22", color = "grey40")

p2


library(scales)

tfun <- scales::pseudo_log_trans(sigma = 12)  # 可调：8~20，越大越“线性”
p2 <- p +
  coord_trans(x = tfun, y = tfun) +
  scale_x_continuous(
    breaks = c(-80, -60, -40, -25, -10, -5, 0, 5),
    labels = label_number()
  ) +
  scale_y_continuous(
    breaks = c(-80, -60, -40, -25, -10, -5, 0, 5),
    labels = label_number()
  )
p2
