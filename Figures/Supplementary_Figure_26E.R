

##---Figure S16b----- hclust-------------------------------------------------
##---Figure S16b----- hclust-------------------------------------------------
##---Figure S16b----- hclust-------------------------------------------------
# ====== packages ======
suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
})

## 1) 读入数据 -------------------------------------------------------------
infile <- "01_heatmap_allTRS_only_8psychiatric.csv"
dat <- readr::read_csv(infile, show_col_types = FALSE)

## 2) 选出 8 个疾病 TRS 列，并按 ID 聚合（重复 ID 取均值） -------------------
trs_cols <- c(
  "SCZ_TRS","BIP_TRS","AN_TRS","TS_TRS",
  "ASD_TRS","ADHD_TRS","MDD_TRS","OCD_TRS"
)
stopifnot(all(trs_cols %in% names(dat)))

mat_by_id <- dat |>
  dplyr::select(ID, dplyr::all_of(trs_cols)) |>
  dplyr::group_by(ID) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(trs_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

mat <- mat_by_id |>
  dplyr::select(-ID) |>
  as.data.frame()

## 3) 计算相关矩阵（Spearman） --------------------------------------------
cor_mat <- cor(
  mat,
  use    = "pairwise.complete.obs",
  method = "spearman"
)

## 4) 自动聚类得到 trait 顺序 ---------------------------------------------
# 距离：1 - correlation（相关越高距离越近）
dist_mat <- as.dist(1 - cor_mat)

# 层次聚类方法可选: "average"(UPGMA), "complete", "ward.D2" 等
hc <- hclust(dist_mat, method = "average")

plot(hc, main = "Trait clustering (1 - Spearman rho)", xlab = "", sub = "")

