
#2026-01-06
rm(list = ls())

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(corrplot)
})

infile <- "./ldsc_brain.xlsx"
df <- read_excel(infile)

trait_order <- c("SCZ","BIP","AN","TS","ASD","ADHD","MDD","OCD")

# -------------------------
# 1) Clean p1/p2 and numeric columns
# -------------------------
dat <- df %>%
  dplyr::select(p1, p2, rg, p) %>%
  dplyr::mutate(
    p1 = str_trim(as.character(p1)),
    p2 = str_trim(as.character(p2)),
    # 去掉可能的引号或奇怪符号
    p1 = str_replace_all(p1, "['\"`]", ""),
    p2 = str_replace_all(p2, "['\"`]", ""),
    # 如果你的表型名里可能有后缀（比如 _TRS），可以在这里去掉：
    p1 = str_replace(p1, "_TRS$", ""),
    p2 = str_replace(p2, "_TRS$", ""),
    rg = suppressWarnings(as.numeric(rg)),
    p  = suppressWarnings(as.numeric(p))
  ) %>%
  dplyr::filter(p1 %in% trait_order, p2 %in% trait_order) %>%
  dplyr::distinct(p1, p2, .keep_all = TRUE) %>%
  dplyr::mutate(
    rg = pmax(pmin(rg, 1), -1)   # clip to [-1, 1]
  )

# -------------------------
# 2) Build 8x8 matrices directly (most robust)
# -------------------------
rg_mat <- matrix(NA_real_, nrow = length(trait_order), ncol = length(trait_order),
                 dimnames = list(trait_order, trait_order))
p_mat  <- matrix(NA_real_, nrow = length(trait_order), ncol = length(trait_order),
                 dimnames = list(trait_order, trait_order))

# 填充 (p1,p2)
for (i in seq_len(nrow(dat))) {
  a <- dat$p1[i]; b <- dat$p2[i]
  rg_mat[a, b] <- dat$rg[i]
  p_mat[a, b]  <- dat$p[i]
  # 同时填充对称位置 (p2,p1)，保证矩阵对称
  rg_mat[b, a] <- dat$rg[i]
  p_mat[b, a]  <- dat$p[i]
}

diag(rg_mat) <- 1
diag(p_mat)  <- 0

# 最后安全检查（避免 corrplot 坐标错误）
stopifnot(all(rownames(rg_mat) == trait_order),
          all(colnames(rg_mat) == trait_order),
          identical(dim(rg_mat), dim(p_mat)))

rg_mat[!is.finite(rg_mat)] <- NA_real_
p_mat[!is.finite(p_mat)]   <- NA_real_
# 可选：把缺失 p 设为 1（表示不显著），避免 insig 标记怪异
p_mat[is.na(p_mat)] <- 1

# -------------------------
# 3) Plot
# -------------------------
col_fun <- colorRampPalette(c("#004E71","white","#B1182D"))

pdf("./LDSC_brain_corrplot_diagCircle.pdf", width = 9, height = 9)

# 1) 上三角：圆圈（包含对角线圆圈）
corrplot(
  rg_mat,
  method = "circle",
  type = "upper",
  col = col_fun(200),
  tl.col = "black",
  tl.srt = 45,
  tl.pos = "lt",
  diag = TRUE,          # ✅ 对角线画大红圈
  cl.pos = "r",
  p.mat = p_mat,
  insig = "label_sig",
  sig.level = c(0.001, 0.01, 0.05),
  pch.cex = 0.9
)

# 2) 下三角：数字（但不画对角线，避免盖住圆圈）
corrplot(
  rg_mat,
  method = "number",
  type = "lower",
  add = TRUE,
  tl.pos = "n",
  cl.pos = "n",
  diag = FALSE,         # ✅ 关键：别在对角线上画 1.00
  number.cex = 0.9,
  col = "black"
)

dev.off()


