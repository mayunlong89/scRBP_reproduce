
##-------Figure 2 ----dot plot-------2 cell types
##-------Figure 2 ----dot plot---monocyte count

library(tidyverse)
library(ggplot2)
library(grid)  # for unit()

# 1) 读入你的长表：Cell_type, RBP, TRS, Dataset, Region
df <- read.table("Figure_2_monocyte_heatmap.txt", sep = "\t", header = TRUE, check.names = FALSE)
# ---- Dot plot with robust RBP clustering (ready to run) ----

df$TRS <- as.numeric(df$TRS)

# 2) 因子顺序（可按需改）
df <- df %>%
  mutate(
    Cell_type = factor(Cell_type, levels = c("Monocytes","T_cells")),
    Region    = factor(Region,    levels = c("5UTR","CDS","Introns","3UTR")),
    Dataset   = factor(Dataset,   levels = c("8.7K","6.5K"))
  )


# 3) —— 稳健的 RBP 聚类顺序：列标准化 + NA=0 + 欧氏距离
cluster_rbps_euclid <- function(df) {
  mat_df <- df %>%
    mutate(key = paste(Cell_type, Region, Dataset, sep = "|")) %>%
    group_by(RBP, key) %>%
    summarise(TRS = mean(TRS, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = key, values_from = TRS)
  
  mat <- mat_df %>% column_to_rownames("RBP") %>% as.matrix()
  if (nrow(mat) < 2) return(rownames(mat))
  
  # 列中心化/标准化，避免某些列方差为0
  col_means <- colMeans(mat, na.rm = TRUE)
  col_sds   <- apply(mat, 2, sd, na.rm = TRUE)
  for (j in seq_len(ncol(mat))) {
    s <- col_sds[j]; m <- col_means[j]
    if (!is.finite(s) || s == 0) {
      mat[, j] <- 0
    } else {
      mat[, j] <- (mat[, j] - m) / s
    }
  }
  mat[!is.finite(mat)] <- 0  # NA -> 0（列中心化后等价于均值）
  
  d  <- stats::dist(mat, method = "euclidean")   # 显式用 stats::，避免方法冲突
  hc <- stats::hclust(d, method = "average")
  rownames(mat)[hc$order]
}

rbp_order <- cluster_rbps_euclid(df)

# 若只剩 1–2 个 RBP，回退为按全局均值排序
if (is.null(rbp_order) || length(rbp_order) < 2) {
  rbp_order <- df %>% group_by(RBP) %>%
    summarise(m = mean(TRS, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(m)) %>% pull(RBP)
}

# 4) 作图数据：大小=|TRS| 截顶；颜色=正/负；补全缺失组合保证网格完整
CAP <- 6  # winsorize 上限（可调 5~8）
dfp <- df %>%
  mutate(
    RBP      = factor(RBP, levels = rbp_order),
    sign     = ifelse(TRS >= 0, "pos", "neg"),
    TRS_size = pmin(abs(TRS), CAP)
  ) %>%
  tidyr::complete(Cell_type, Region, Dataset, RBP,
                  fill = list(TRS = NA_real_, TRS_size = NA_real_, sign = NA_character_))

rbp_order <- c("DDX21","DDX3X","TRA2A","FAM120A","QKI","EIF4A1","WTAP","HNRNPH2",
               "YTHDF3","RBM47","YBX3","PABPC1","XRN2","CELF1","RYBP","VIM","NIPBL",
               "RPRD2","TAF15","AKAP1","IGF2BP1","LIN28A","PUM2","ILF3","ZNF148")
               
rbp_order_full <- c(rbp_order, setdiff(sort(unique(df$RBP)), rbp_order))
dfp <- dfp %>% mutate(RBP = factor(RBP, levels = rbp_order_full))


# 5) 画图：纵向分面 = Cell_type → Region；面板内 y 轴为 Dataset
p <- ggplot(dfp, aes(x = RBP, y = Dataset)) +
  geom_point(aes(size = TRS_size, fill = sign),
             shape = 21, color = "white", stroke = 0.3, alpha = 0.9, na.rm = TRUE) +
  scale_x_discrete(limits = rbp_order_full, drop = FALSE) + 
  scale_fill_manual(values = c(neg = "#4C78A8", pos = "#E45756"),
                    na.value = "grey85", name = "Sign") +
  scale_size_area(name = "|TRS|", max_size = 6,
                  breaks = c(1, 3, 5), limits = c(0, CAP),
                  oob = scales::squish) +
  facet_grid(rows = vars(Cell_type, Region), cols = vars(), switch = "y") +
  labs(x = "RBP", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(size = 0.3, colour = "grey85"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    legend.position = "bottom",
    panel.spacing.y = unit(6, "pt")
  )

print(p)

