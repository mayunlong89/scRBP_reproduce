#2026-03-04

setwd("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/09_scRBP_infer_grn_adultbrain_fetal_brain_Lung")
# --------- 安装依赖（首次运行需要） ----------
# install.packages(c("ggtern","dplyr","purrr","stringr"))

library(ggtern)
library(dplyr)
library(purrr)
library(stringr)

# --------- 1) 读 GMT 为 list[name -> unique gene vec] ----------
read_gmt <- function(path, canonical_by = c("first","second"),
                     toupper_genes = TRUE, min_genes = 1) {
  canonical_by <- match.arg(canonical_by)
  lines <- readLines(path, warn = FALSE)
  out <- list()
  for (ln in lines) {
    if (!nzchar(ln)) next
    f <- strsplit(ln, "\t")[[1]]
    if (length(f) < 3) next
    # GMT: 第1列集合名, 第2列描述, 第3列起为基因
    nm_raw1 <- f[[1]]
    nm_raw2 <- f[[2]]
    nm <- if (canonical_by == "first") nm_raw1 else nm_raw2
    # 去掉常见后缀以得到 RBP 名（按需）
    nm <- sub("_regulon$", "", nm, ignore.case = TRUE)
    genes <- unique(f[-c(1,2)])
    if (toupper_genes) genes <- toupper(genes)
    genes <- genes[nzchar(genes)]
    if (length(genes) >= min_genes) out[[nm]] <- genes
  }
  out
}

# --------- 2) 分摊计分 + 并集归一化 的三元比例 ----------
frac_props <- function(g1, g2, g3) {
  U <- union(union(g1, g2), g3)
  if (length(U) == 0) return(c(PBMC=0, Fetal=0, Adult=0, n_union=0))
  m1 <- U %in% g1; m2 <- U %in% g2; m3 <- U %in% g3
  k  <- (m1 + m2 + m3)             # 每个基因出现于几份集合
  credit1 <- sum(m1 / k); credit2 <- sum(m2 / k); credit3 <- sum(m3 / k)
  tot <- credit1 + credit2 + credit3   # == length(U)
  c(PBMC = credit1/tot, Fetal = credit2/tot, Adult = credit3/tot, n_union = length(U))
}

# （可选）简单计数比例，不分摊、直接相加做分母
naive_props <- function(g1, g2, g3) {
  tot <- length(g1) + length(g2) + length(g3)
  if (tot == 0) return(c(PBMC=0, Fetal=0, Adult=0, n_union=0))
  U <- length(union(union(g1,g2), g3))
  c(PBMC = length(g1)/tot, Fetal = length(g2)/tot, Adult = length(g3)/tot, n_union = U)
}

# --------- 3) 主函数：读入、取交集、算比例、作图 ----------
ternary_from_gmts <- function(pbmc_gmt, fetal_gmt, adult_gmt,
                              canonical_by = "first",   # regulon 名来自第1列；若你需要用第2列就设 "second"
                              min_genes = 1,
                              method = c("fractional","naive"),
                              top_label_n = 0           # >0 时，给偏向度最大的前 n 个打标签
) {
  method <- match.arg(method)
  PB <- read_gmt(pbmc_gmt, canonical_by = canonical_by, min_genes = min_genes)
  FT <- read_gmt(fetal_gmt, canonical_by = canonical_by, min_genes = min_genes)
  AD <- read_gmt(adult_gmt, canonical_by = canonical_by, min_genes = min_genes)

  common_ids <- Reduce(intersect, list(names(PB), names(FT), names(AD)))
  if (length(common_ids) == 0) stop("三个 GMT 没有共同的 regulon 名。请检查命名是否一致。")

  calc_fun <- if (method == "fractional") frac_props else naive_props

  df <- map_dfr(common_ids, function(id){
    p <- calc_fun(PB[[id]], FT[[id]], AD[[id]])
    tibble(rbp = id, PBMC = p[1], Fetal = p[2], Adult = p[3], n_union = as.numeric(p[4]))
  })

  # 可选：给“最偏向”的若干个加标签
  if (top_label_n > 0) {
    df <- df %>%
      mutate(max_share = pmax(PBMC,Fetal,Adult)) %>%
      arrange(desc(max_share)) %>%
      mutate(label = if_else(row_number() <= top_label_n, rbp, ""))
  } else {
    df$label <- ""
  }

  p <- ggtern(df, aes(x = PBMC, y = Fetal, z = Adult)) +
    geom_point(aes(size = n_union), alpha = 0.85) +
    geom_text(aes(label = label), vjust = -0.4, size = 3, show.legend = FALSE) +
    labs(T = "Fetal", L = "PBMC", R = "Adult",
         title = paste0("Ternary of shared RBP-regulons (", method, ")"),
         size = "Union size") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

  list(data = df, plot = p)
}

# --------- 4) 运行示例 ----------
# 把路径替换为你的三个 GMT（同一个 mRNA 区域的一组）

##---5UTR
lung_gmt  <- "scRBP_prune_matched_results_min10genes_5UTR_lung.gmt"
fetal_gmt <- "scRBP_prune_matched_results_min10genes_5UTR_fetal.gmt"
adult_gmt <- "scRBP_prune_matched_results_min10genes.symbol_default_palmer_adult_100Kcells_5UTR.gmt"


##---CDS
lung_gmt  <- "scRBP_prune_matched_results_min10genes_CDS_lung.gmt"
fetal_gmt <- "scRBP_prune_matched_results_min10genes_CDS_fetal.gmt"
adult_gmt <- "scRBP_prune_matched_results_min10genes.symbol_default_palmer_adult_100Kcells_CDS.gmt"


##---Introns
lung_gmt  <- "scRBP_prune_matched_results_min10genes_Introns_lung.gmt"
fetal_gmt <- "scRBP_prune_matched_results_min10genes_Introns_fetal.gmt"
adult_gmt <- "scRBP_prune_matched_results_min10genes.symbol_default_palmer_adult_100Kcells_Introns.gmt"

##---3UTR
lung_gmt  <- "scRBP_prune_matched_results_min10genes_3UTR_lung.gmt"
fetal_gmt <- "scRBP_prune_matched_results_min10genes_3UTR_fetal.gmt"
adult_gmt <- "scRBP_prune_matched_results_min10genes.symbol_default_palmer_adult_100Kcells_3UTR.gmt"


##---All four regions
lung_gmt  <- "scRBP_allRegions_min10genes.symbol_default_4regions_lung.gmt"
fetal_gmt <- "scRBP_allRegions_min10genes.symbol_default_4regions_Fetal_brain.gmt"
adult_gmt <- "scRBP_allRegions_min10genes.symbol_default_4region_adult_brain.gmt"






res <- ternary_from_gmts(lung_gmt, fetal_gmt, adult_gmt,
                         canonical_by = "first",   # 若你的集合名在第2列（比如 AATF_regulon），设为 "second"
                         min_genes = 5,            # 例如仅保留 >=5 基因的 regulon
                         method = "fractional",    # 推荐：分摊计分
                         top_label_n = 10)         # 标注前 10 个最偏向的点

print(res$plot)
head(res$data)

data2 <- as.data.frame(res$data)


#write.csv(res$data, file="04_ternary_plot_5UTR.csv",quote=F,row.names = F)
#write.csv(res$data, file="04_ternary_plot_CDS.csv",quote=F,row.names = F)
#write.csv(res$data, file="04_ternary_plot_Introns.csv",quote=F,row.names = F)
#write.csv(res$data, file="04_ternary_plot_3UTR.csv",quote=F,row.names = F)

library(dplyr)

df <- res$data %>%   # 包含列: rbp, PBMC, Fetal, Adult, n_union
  mutate(
    winner    = c("PBMC","Fetal","Adult")[max.col(as.matrix(select(., PBMC, Fetal, Adult)))],
    max_share = pmax(PBMC,Fetal,Adult),
    class = case_when(
      max_share >= 0.60 ~ paste0(winner, "-dominant"),
      max_share >= 0.45 ~ paste0(winner, "-leaning"),
      TRUE              ~ "balanced"
    ),
    dist_center = sqrt( (PBMC-1/3)^2 + (Fetal-1/3)^2 + (Adult-1/3)^2 )
  )

# 汇总有多少是 Adult/Fetal/PBMC 主导或偏向
df %>% count(class, sort = TRUE)

# 若想看“脑总体” vs PBMC，合并 Fetal+Adult：
df %>% mutate(brain = Fetal + Adult) %>%
  summarize(median_brain = median(brain), mean_brain = mean(brain))

library(dplyr)

df %>%
  mutate(brain = Adult + Fetal,
         af_diff = abs(Adult - Fetal)) %>%
  summarize(
    median_brain = median(brain),
    mean_brain   = mean(brain),
    prop_brain_ge_2_3 = mean(brain >= 2/3),          # 脑主导的比例
    median_afdiff = median(af_diff),                  # 成人-胎脑差的中位数
    prop_close_brain = mean(brain >= 2/3 & af_diff <= 0.10),  # 脑主导且A≈F的比例（阈值可改）
    rho_AF = cor(Adult, Fetal, method = "spearman"),  # 成人 vs 胎脑相关
    p_AF   = cor.test(Adult, Fetal, method = "spearman")$p.value
  )


###_-------标柱---dominant regulons for each tissues
###_-------标柱---dominant regulons for each tissues
###_-------标柱---dominant regulons for each tissues

##---marked top 5 RBP-regulons
library(ggtern)
library(dplyr)
library(scales)
library(grid)

# 计算赢家与每个组织的Top5
df_top <- res$data %>%
  mutate(
    winner    = c("PBMC","Fetal","Adult")[max.col(as.matrix(select(., PBMC, Fetal, Adult)))],
    max_share = pmax(PBMC, Fetal, Adult)
  ) %>%
  group_by(winner) %>%
  arrange(desc(max_share), desc(n_union), .by_group = TRUE) %>%
  mutate(rank_within = row_number(),
         label = if_else(rank_within <= 5, rbp, "")) %>%
  ungroup()

#pal <- c(PBMC = "#2B6CB0", Fetal = "#2F855A", Adult = "#C53030")
pal <- c(PBMC = "#74add1", Fetal = "#ffbba4", Adult = "#84A8A6")
pal <- c(PBMC = "#108DDF", Fetal = "#F47771", Adult = "#6D9B91")

p <- ggtern(df_top, aes(x = PBMC, y = Fetal, z = Adult)) +
  geom_point(aes(color = winner, size = n_union), alpha = 0.9) +

  ## 仅给每个组织Top5打标签；用 geom_text + check_overlap 避免严重重叠
  geom_text(
    data = subset(df_top, label != ""),
    aes(label = label, color = winner),
    size = 3, check_overlap = F, show.legend = FALSE
  ) +

  scale_T_continuous(limits = c(0,1), breaks = seq(0,1,0.2), labels = seq(0,100,20)) +
  scale_L_continuous(limits = c(0,1), breaks = seq(0,1,0.2), labels = seq(0,100,20)) +
  scale_R_continuous(limits = c(0,1), breaks = seq(0,1,0.2), labels = seq(0,100,20)) +

  scale_color_manual(values = pal, name = "Dominant tissue") +
  scale_size_continuous(range = c(2,7), name = "Union size") +

  theme_bw() +
  theme(
    tern.panel.grid.major = element_line(color = "grey85"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),  # ← 改这里
    tern.axis.ticks = element_line(linewidth = 0.6),
    tern.axis.ticks.length.major = unit(6, "pt"),
    tern.axis.text.T = element_text(color = "black"),
    tern.axis.text.L = element_text(color = "black"),
    tern.axis.text.R = element_text(color = "black")
  ) +
  labs(title = "Ternary of shared RBP-regulons (fractional share)",
       T = "Fetal", L = "PBMC", R = "Adult") +
  geom_Lline(Lintercept = 1/3, linetype = 2, linewidth = 0.4, color = "red")  # PBMC = 1/3

print(p)





###_-------标柱---dominant regulons for each tissues
###_-------标柱---dominant regulons for each tissues
###_-------标柱---dominant regulons for each tissues

##---marked top 5 RBP-regulons
library(ggtern)
library(dplyr)
library(scales)
library(grid)

# 计算赢家与每个组织的Top5
df_top <- res$data %>%
  mutate(
    winner    = c("PBMC","Fetal","Adult")[max.col(as.matrix(select(., PBMC, Fetal, Adult)))],
    max_share = pmax(PBMC, Fetal, Adult)
  ) %>%
  group_by(winner) %>%
  arrange(desc(max_share), desc(n_union), .by_group = TRUE) %>%
  mutate(rank_within = row_number(),
         label = if_else(rank_within <= 5, rbp, "")) %>%
  ungroup()

#pal <- c(PBMC = "#2B6CB0", Fetal = "#2F855A", Adult = "#C53030")
#pal <- c(PBMC = "#74add1", Fetal = "#ffbba4", Adult = "#84A8A6")
pal <- c(PBMC = "#FFB000", Fetal = "#A6D854", Adult = "#0B6B3A")


p <- ggtern(df_top, aes(x = PBMC, y = Fetal, z = Adult)) +
  geom_point(aes(color = winner, size = n_union), alpha = 0.9) +
  
  ## 仅给每个组织Top5打标签；用 geom_text + check_overlap 避免严重重叠
  geom_text(
    data = subset(df_top, label != ""),
    aes(label = label, color = winner),
    size = 3, check_overlap = F, show.legend = FALSE
  ) +
  
  scale_T_continuous(limits = c(0,1), breaks = seq(0,1,0.2), labels = seq(0,100,20)) +
  scale_L_continuous(limits = c(0,1), breaks = seq(0,1,0.2), labels = seq(0,100,20)) +
  scale_R_continuous(limits = c(0,1), breaks = seq(0,1,0.2), labels = seq(0,100,20)) +
  
  scale_color_manual(values = pal, name = "Dominant tissue") +
  scale_size_continuous(range = c(2,7), name = "Union size") +
  
  theme_bw() +
  theme(
    tern.panel.grid.major = element_line(color = "grey85"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),  # ← 改这里
    tern.axis.ticks = element_line(linewidth = 0.6),
    tern.axis.ticks.length.major = unit(6, "pt"),
    tern.axis.text.T = element_text(color = "black"),
    tern.axis.text.L = element_text(color = "black"),
    tern.axis.text.R = element_text(color = "black")
  ) +
  labs(title = "Ternary of shared RBP-regulons (fractional share)",
       T = "Fetal", L = "PBMC", R = "Adult") +
  geom_Lline(Lintercept = 1/3, linetype = 2, linewidth = 0.4, color = "red")  # PBMC = 1/3

print(p)

