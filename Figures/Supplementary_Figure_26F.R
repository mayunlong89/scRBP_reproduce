
#2026-01-28
setwd("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/03_scRBP_trs_fetal_brain_15times/01_Final_scRBP_TRS_fetal_brain_50times")
library(tidyverse)


#1) 读入数据
df <- read.csv("01_heatmap_allTRS_only_8psychiatric.csv", stringsAsFactors = FALSE)

# 确保列名一致（按你截图：Cell_types / SCZ_TRS / BIP_TRS）
stopifnot(all(c("Cell_types", "SCZ_TRS", "BIP_TRS") %in% colnames(df)))


#2) 按 cell type 计算 correlation（Pearson / Spearman 都给你）
#2.1 计算每个 cell type 的相关系数 + p 值 + CI（推荐）

library(dplyr)

# 一个小函数：对每个 cell type 做 cor.test
cor_one_group <- function(x, y, method = "pearson"){
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if(n < 3){
    return(tibble(n = n, rho = NA_real_, p = NA_real_,
                  ci_low = NA_real_, ci_high = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x, y, method = method))
  tibble(
    n = n,
    rho = unname(ct$estimate),
    p = ct$p.value,
    ci_low = if(!is.null(ct$conf.int)) ct$conf.int[1] else NA_real_,
    ci_high= if(!is.null(ct$conf.int)) ct$conf.int[2] else NA_real_
  )
}

# 你可以选 pearson 或 spearman
method_use <- "pearson"   # 或 "spearman"

cor_by_ct <- df %>%
  group_by(Cell_types) %>%
  summarise(
    tmp = list(cor_one_group(SCZ_TRS, BIP_TRS, method = method_use)),
    .groups = "drop"
  ) %>%
  tidyr::unnest(tmp) %>%
  mutate(
    padj = p.adjust(p, method = "BH"),
    pair = "SCZ_TRS vs BIP_TRS"
  )

cor_by_ct

#议你设一个最小样本数阈值，比如每个 cell type 至少 10 个点才画：
min_n <- 10
cor_by_ct2 <- cor_by_ct %>% filter(n >= min_n)


#3) 画 heatmap（ggplot2：行=Cell_types，列=pair）
#3.1 基础版：相关系数热图 + 数值标注

library(ggplot2)
library(scales)

# 如果你想固定 cell type 的显示顺序（不想就注释掉）
celltype_order <- cor_by_ct2 %>%
  arrange(desc(rho)) %>%
  pull(Cell_types)

cor_by_ct2 <- cor_by_ct2 %>%
  mutate(Cell_types = factor(Cell_types, levels = celltype_order))

p_heat <- ggplot(cor_by_ct2, aes(x = Cell_types, y = pair, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(is.na(rho), "", sprintf("%.2f", rho))), size = 3.2) +
  scale_fill_gradient2(
    low  = "#2166AC",
    mid  = "white",
    high = "#B2182B",
    midpoint = 0.5,                # 中心点（相关性中位数）
    limits   = c(0.5, 0.8),        # 视你实际结果可再微调
    name     = "Spearman \u03c1"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

p_heat 

p_heat_h <- p_heat + coord_flip()
p_heat_h










cor_by_ct2 <- cor_by_ct2 %>%
  mutate(sig = case_when(
    is.na(padj) ~ "",
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    TRUE ~ ""
  ),
  label = ifelse(is.na(rho), "", sprintf("%.2f%s\nn=%d", rho, sig, n))
  )

p_heat2 <- ggplot(cor_by_ct2, aes(x = pair, y = Cell_types, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 3.0, lineheight = 0.95) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1), oob = squish,
    name = paste0(method_use, " r")
  ) +
  theme_bw(base_size = 12) +
  theme(axis.title = element_blank(), panel.grid = element_blank())

p_heat2
