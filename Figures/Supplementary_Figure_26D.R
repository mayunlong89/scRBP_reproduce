



#2026-01-28
###----Figure S16a
setwd("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/03_scRBP_trs_fetal_brain_15times/01_Final_scRBP_TRS_fetal_brain_50times")

library(pheatmap)

# ====== packages ======
suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
})

## 1) 读入数据 -------------------------------------------------------------
# 改成你的文件名
infile <- "01_heatmap_allTRS_only_8psychiatric.csv"
dat <- readr::read_csv(infile, show_col_types = FALSE)

dat <- dat[,-2]

# 把 ID 设成行名（可选）
# rownames(dat) <- dat$ID

## 2) 选出 8 个疾病的 TRS 列，并按 ID 聚合（万一有重复 ID 就取均值） ----
# 想要的疾病及顺序（也将决定相关性热图的行列顺序）
trait_order <- c("SCZ_TRS", "BIP_TRS",  "ASD_TRS", "ADHD_TRS",
                 "MDD_TRS", "AN_TRS",
                 "TS_TRS","OCD_TRS")

trs_cols <- trait_order
stopifnot(all(trs_cols %in% names(dat)))

mat_by_id <- dat |>
  dplyr::select(ID, dplyr::all_of(trs_cols)) |>
  dplyr::group_by(ID) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(trs_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# 数值矩阵：行为 ID，列为 8 个疾病
mat <- mat_by_id |>
  dplyr::select(-ID) |>
  as.data.frame()

## 3) 计算相关矩阵（Spearman / Pearson 都行） -----------------------------
cor_mat <- cor(
  mat,
  use    = "pairwise.complete.obs",
  method = "spearman"   # 如想改 Pearson 就写 "pearson"
)

# 确保列顺序是我们想要的顺序
cor_mat <- cor_mat[trait_order, trait_order]

## 4) 整理成下三角数据框，用 ggplot 画图 ---------------------------------
cor_df <- reshape2::melt(
  cor_mat,
  varnames  = c("Var1", "Var2"),
  value.name = "rho"
)

# 把 Var1/Var2 设为 factor，保证顺序一致
cor_df <- cor_df |>
  dplyr::mutate(
    Var1 = factor(Var1, levels = trait_order),
    Var2 = factor(Var2, levels = trait_order)
  )

# 只保留下三角（Var1 在 Var2 的“下方”）
cor_df_lower <- cor_df |>
  dplyr::filter(as.numeric(Var1) > as.numeric(Var2))

## 5) 画图（下三角相关性热图 + 数字） -------------------------------------
p <- ggplot(cor_df_lower, aes(x = Var1, y = Var2, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
  # 颜色范围可以按需要调，这里类似你之前代码的设定
  scale_fill_gradient2(
    low  = "#2166AC",
    mid  = "white",
    high = "#B2182B",
    midpoint = 0.5,                # 中心点（相关性中位数）
    limits   = c(0.5, 0.8),        # 视你实际结果可再微调
    name     = "Spearman \u03c1"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, color = "red"
    ),
    axis.text.y = element_text(color = "red"),
    axis.title  = element_blank(),
    panel.grid  = element_blank()
  )

p


