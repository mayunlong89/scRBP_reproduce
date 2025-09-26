
#---2025-09-05

# ---- 1) 读入数据（两种方式二选一）----

## 方式A：若你有一个Excel/CSV，列名=集合名，列里是条目（允许为空格/NA）
# path <- "your_file.xlsx"  # 或 "your_file.csv"
# if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
#   library(readxl)
#   df <- readxl::read_excel(path)
# } else {
#   library(readr)
#   df <- readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()))
# }

## 方式B：直接在 R 里定义（示例，替换为你的6列条目）
# 每个向量就是一列里的条目
Monocyte_count_6.5K <- c(
  "PABPC1_Introns","RBM47_Introns","XRN2_Introns","YBX3_Introns","HNRNPC_Introns",
  "ZNF148_Introns","CELF1_3UTR","DDX21_5UTR","DDX3X_5UTR","TRA2A_5UTR","RYBP_3UTR",
  "ILF3_3UTR","TAF15_3UTR","NIPBL_5UTR","FAM120A_5UTR","QKI_5UTR","IGF2BP1_5UTR",
  "PUM2_5UTR","EIF4A1_5UTR","AKAP1_5UTR","LIN28A_5UTR","WTAP_5UTR","RPRD2_5UTR",
  "HNRNPH2_5UTR","YTHDF3_5UTR","DDX21_3UTR","DDX3X_3UTR","TRA2A_3UTR"
)
Monocyte_count_8.7K <- c(
  "PABPC1_Introns","RBM47_Introns","XRN2_Introns","YBX3_Introns","HNRNPC_Introns",
  "ZNF148_Introns","CELF1_3UTR","DDX21_5UTR","DDX3X_5UTR","TRA2A_5UTR","RYBP_3UTR",
  "ILF3_3UTR","TAF15_3UTR","NIPBL_5UTR","FAM120A_5UTR","QKI_5UTR","IGF2BP1_5UTR",
  "PUM2_5UTR","EIF4A1_5UTR","AKAP1_5UTR","LIN28A_5UTR","WTAP_5UTR","RPRD2_5UTR",
  "HNRNPH2_5UTR","YTHDF3_5UTR","DDX21_3UTR","DDX3X_3UTR","TRA2A_3UTR"
)
lymp_count_6.5K   <- c("EEF2_5UTR","TAF15_5UTR","NIPBL_5UTR","DDX3X_5UTR","QKI_5UTR",
                       "LARP7_Introns","SRSF1_Introns","RTCA_Introns","DEK_Introns","RBM47_Introns")
lymp_count_8.7K   <- c("EEF2_5UTR","TAF15_5UTR","NIPBL_5UTR","DDX3X_5UTR","QKI_5UTR",
                       "LARP7_Introns","SRSF1_Introns","RTCA_Introns","DEK_Introns","RBM47_Introns")
lymp_percent_6.5K <- c("NIPBL_5UTR","NXF1_5UTR","PUM1_5UTR","TNRC6C_5UTR","HNRNPC_Introns",
                       "HNRNPC_5UTR","FMR1_5UTR","QKI_5UTR","DEK_Introns")
lymp_percent_8.7K <- c("NIPBL_5UTR","NXF1_5UTR","PUM1_5UTR","TNRC6C_5UTR","HNRNPC_5UTR",
                       "FMR1_5UTR","QKI_5UTR","DEK_Introns")

# 组装成 named list（如果用了方式A，就把 df 转为 sets）
sets <- list(
  Monocyte_count_6.5K = Monocyte_count_6.5K,
  Monocyte_count_8.7K = Monocyte_count_8.7K,
  lymp_count_6.5K     = lymp_count_6.5K,
  lymp_count_8.7K     = lymp_count_8.7K,
  lymp_percent_6.5K   = lymp_percent_6.5K,
  lymp_percent_8.7K   = lymp_percent_8.7K
)

# 若你用方式A（df），请用下面两行把每列转集合（自动去空、去重）：
# sets <- lapply(df, function(x) unique(na.omit(trimws(x[x != ""]))))
# names(sets) <- make.unique(names(sets))

# ---- 2) 计算 Jaccard 相似度矩阵 ----
jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(1)  # 两个空集时定义为1
  length(intersect(a, b)) / length(union(a, b))
}

n <- length(sets)
J <- matrix(0, n, n, dimnames = list(names(sets), names(sets)))
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    J[i, j] <- jaccard(sets[[i]], sets[[j]])
  }
}

# （可选）看看集合大小与交集数
sizes <- sapply(sets, length)
overlaps <- outer(names(sets), names(sets),
                  Vectorize(function(x,y) length(intersect(sets[[x]], sets[[y]]))))
dimnames(overlaps) <- list(names(sets), names(sets))
sizes; overlaps  # 查看

# ---- 3) 画热图 ----
library(pheatmap)
pheatmap(J,
         cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE, number_format = "%.2f",
         border_color = NA, legend = TRUE,
         main = "Pairwise Jaccard similarity")

# （可选）保存图片
# ggsave("jaccard_heatmap.png", width = 6, height = 5, dpi = 300)




library(RColorBrewer)
warm_cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)

pheatmap(J,
         color = warm_cols,
         breaks = seq(0, 1, length.out = 101),  # 0~1 等距映射
         cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE, number_format = "%.2f",
         border_color = NA, main = "Pairwise Jaccard similarity")



warm_cols <- colorRampPalette(c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272",
                                "#FB6A4A","#EF3B2C","#CB181D","#99000D"))(120)

pheatmap(J, color = warm_cols, breaks = seq(0,1,length.out=121),
         cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE, number_format = "%.2f",
         border_color = NA, main = "Pairwise Jaccard similarity")





# ---- Packages ----
library(pheatmap)
library(RColorBrewer)

# ---- 1) 定义6个集合（每个向量是一列的条目；可随时增删） ----
Monocyte_count_6.5K <- c(
  "PABPC1_Introns","RBM47_Introns","XRN2_Introns","YBX3_Introns","HNRNPC_Introns",
  "ZNF148_Introns","CELF1_3UTR","DDX21_5UTR","DDX3X_5UTR","TRA2A_5UTR","RYBP_3UTR",
  "ILF3_3UTR","TAF15_3UTR","NIPBL_5UTR","FAM120A_5UTR","QKI_5UTR","IGF2BP1_5UTR",
  "PUM2_5UTR","EIF4A1_5UTR","AKAP1_5UTR","LIN28A_5UTR","WTAP_5UTR","RPRD2_5UTR",
  "HNRNPH2_5UTR","YTHDF3_5UTR","DDX21_3UTR","DDX3X_3UTR","TRA2A_3UTR"
)

Monocyte_count_8.7K <- c(
  "PABPC1_Introns","RBM47_Introns","XRN2_Introns","YBX3_Introns","HNRNPC_Introns",
  "ZNF148_Introns","CELF1_3UTR","DDX21_5UTR","DDX3X_5UTR","TRA2A_5UTR","RYBP_3UTR",
  "ILF3_3UTR","TAF15_3UTR","NIPBL_5UTR","FAM120A_5UTR","QKI_5UTR","IGF2BP1_5UTR",
  "PUM2_5UTR","EIF4A1_5UTR","AKAP1_5UTR","LIN28A_5UTR","WTAP_5UTR","RPRD2_5UTR",
  "HNRNPH2_5UTR","YTHDF3_5UTR","DDX21_3UTR","DDX3X_3UTR","TRA2A_3UTR"
)

lymp_count_6.5K <- c(
  "EEF2_5UTR","TAF15_5UTR","NIPBL_5UTR","DDX3X_5UTR","QKI_5UTR",
  "LARP7_Introns","SRSF1_Introns","RTCA_Introns","DEK_Introns","RBM47_Introns",
  "HNRNPC_Introns",          # 新增
  "TAF15_3UTR",              # 新增
  "FMR1_5UTR","NIPBL_CDS","DDX3X_3UTR"
)

lymp_count_8.7K <- c(
  "EEF2_5UTR","TAF15_5UTR","NIPBL_5UTR","DDX3X_5UTR","QKI_5UTR",
  "LARP7_Introns","SRSF1_Introns","RTCA_Introns","DEK_Introns","RBM47_Introns",
  "HNRNPC_Introns",          # 新增
  "TAF15_3UTR",              # 新增
  "FMR1_5UTR","NIPBL_CDS","DDX3X_3UTR"
)

lymp_percent_6.5K <- c(
  "NIPBL_5UTR","NXF1_5UTR","PUM1_5UTR","TNRC6C_5UTR",
  "HNRNPC_Introns","HNRNPC_5UTR","FMR1_5UTR","QKI_5UTR","DEK_Introns",
  "EEF2_5UTR","SRSF1_CDS",   # 新增
  "DDX3X_3UTR","DDX3X_5UTR","RBM47_Introns","FMR1_3UTR"
)

lymp_percent_8.7K <- c(
  "NIPBL_5UTR","NXF1_5UTR","PUM1_5UTR","TNRC6C_5UTR",
  "HNRNPC_5UTR","FMR1_5UTR","QKI_5UTR","DEK_Introns",
  "EEF2_5UTR","SRSF1_CDS",   # 新增
  "DDX3X_3UTR","DDX3X_5UTR","RBM47_Introns","FMR1_3UTR"
)

# ---- 2) 组装为 named list，并做清洗（去重、去空格与空值） ----
sets <- list(
  Monocyte_count_6.5K = Monocyte_count_6.5K,
  Monocyte_count_8.7K = Monocyte_count_8.7K,
  lymp_count_6.5K     = lymp_count_6.5K,
  lymp_count_8.7K     = lymp_count_8.7K,
  lymp_percent_6.5K   = lymp_percent_6.5K,
  lymp_percent_8.7K   = lymp_percent_8.7K
)
clean_vec <- function(x) unique(na.omit(trimws(x[x != ""])))
sets <- lapply(sets, clean_vec)

# ---- 3) 计算两两 Jaccard 相似度、交集数与集合大小 ----
jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(1)
  length(intersect(a, b)) / length(union(a, b))
}

n <- length(sets)
J <- matrix(0, n, n, dimnames = list(names(sets), names(sets)))
for (i in seq_len(n)) for (j in seq_len(n)) J[i, j] <- jaccard(sets[[i]], sets[[j]])

sizes <- sapply(sets, length)
overlaps <- outer(names(sets), names(sets),
                  Vectorize(function(x,y) length(intersect(sets[[x]], sets[[y]]))))
dimnames(overlaps) <- list(names(sets), names(sets))

cat("Sizes:\n"); print(sizes)
cat("\nPairwise overlaps (counts):\n"); print(overlaps)
cat("\nJaccard matrix:\n"); print(round(J, 3))

# ---- 4) 纯暖色渐变 + 用 (1 - J) 作为聚类距离 ----
warm_cols <- colorRampPalette(c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272",
                                "#FB6A4A","#EF3B2C","#CB181D","#99000D"))(200)
warm_cols <- colorRampPalette(c("grey","#FFF5F0","#FEE0D2","#FCBBA1","#FCBBA1","#FC9272","#FC9272",
                                "#FB6A4A"))(200)
warm_cols <- colorRampPalette(c("#FFF5F0","#FEE0D2","#FCBBA1","#FCBBA1","#FC9272","#FC9272",
                                "#FB6A4A"))(200)

# 用 1-J 做层次聚类，便于和相似度一致地解释
d  <- as.dist(1 - J)
hc <- hclust(d, method = "average")

pheatmap(J,
         color = warm_cols,
         breaks = seq(0, 1, length.out = length(warm_cols) + 1),
         cluster_rows = hc, cluster_cols = hc,   # 使用上面的距离/树
         display_numbers =F, number_format = "%.2f",
         number_color = "black", fontsize_number = 9,
         border_color = NA, angle_col = 45,
         main = "Pairwise Jaccard similarity (warm palette)")



pheatmap(J,
         color = warm_cols,
         breaks = seq(0,1,length.out=length(warm_cols)+1),
         cluster_rows = hc, cluster_cols = hc,
         display_numbers = F, number_format = "%.2f",
         number_color = "black", fontsize_number = 9,
         border_color = "white",      # ← 每个格子白色边框
         angle_col = 45, main = "Pairwise Jaccard similarity")



# ---- 5) 可选：导出结果 ----
# write.csv(J,        "jaccard_matrix.csv")
# write.csv(overlaps, "overlap_counts.csv")
# png("jaccard_heatmap_warm.png", width = 1800, height = 1500, res = 200)
# pheatmap(J, color = warm_cols,
#          breaks = seq(0, 1, length.out = length(warm_cols) + 1),
#          cluster_rows = hc, cluster_cols = hc,
#          display_numbers = TRUE, number_format = "%.2f",
#          number_color = "black", fontsize_number = 9,
#          border_color = NA, angle_col = 45,
#          main = "Pairwise Jaccard similarity (warm palette)")
# dev.off()



# ---- 数据 ----
txt <- "
ZNF589   9.817870786
VPS53    8.186819467
UBXN11   43.12496486
TRAF3    10.58236226
TNS1     6.36713796
TNPO1    14.58526429
TFCP2    5.630969778
STXBP6   50
SNX1     50
SLC45A4  13.58479555
RBPJ     8.868702203
RAF1     11.3705904
RAD52    16.57228862
RAB10    9.020269569
POU2F2   13.55424016
POU2F1   20.25125346
OPRM1    8.549289122
NCBP3    9.786747948
KLF7     10.49227902
ITGAL    38.6644053
IREB2    9.338187314
IGF2R    8.488250289
FRS2     6.271808601
EVI5     50
E2F3     11.17652577
CELF1    18.00436481
BCL6     13.98046832
ATF7     44.43427021
ANKRD44  16.65619767
"

dat <- read.table(text = txt, header = FALSE,
                  col.names = c("Gene","Value"),
                  stringsAsFactors = FALSE)
dat <- dat[order(dat$Value, decreasing = TRUE), ]

# 矩阵（1 列），行名为基因
mat <- matrix(dat$Value, ncol = 1)
rownames(mat) <- dat$Gene
colnames(mat) <- "Score"

# ---- 颜色（纯暖色）----
library(pheatmap)
warm_cols <- colorRampPalette(c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272",
                                "#FB6A4A","#EF3B2C","#CB181D","#99000D"))(200)

# ---- 默认：log10 色阶（更均衡）----
mat_log <- mat   # +1 避免 log(0)
pheatmap(mat_log,
         color = warm_cols,
         breaks = seq(min(mat_log), max(mat_log), length.out = length(warm_cols)+1),
         cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = "white",           # 每格白线分隔
         display_numbers = FALSE,
         main = "Heatmap (log10 scale, warm palette)")

# ---- 可选：线性色阶（若不想压缩动态范围则用这段）----
# pheatmap(mat,
#          color = warm_cols,
#          breaks = seq(min(mat), max(mat), length.out = length(warm_cols)+1),
#          cluster_rows = FALSE, cluster_cols = FALSE,
#          border_color = "white",
#          display_numbers = FALSE,
#          main = "Heatmap (linear, warm palette)")

