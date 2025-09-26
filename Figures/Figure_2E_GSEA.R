#2025-09-07

#monocyte count


setwd("/Users/mayunlong/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/02_scRBP_trs_10blood_cell_traits_100times/00-OpenTargets_monocyte_count/")

# =========================
# GSEA: OpenTargets CSV  ×  RBP regulons (GMT)
# =========================
# 依赖
#BiocManager::install("fgsea")

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(stringr)
  library(ggplot2)
})

# ---- 读取 OpenTargets CSV：给出基因列与分数字段 ----
# 会做：去空白、去 NA、同基因多行时取最大分数、按分数降序
read_opentargets_csv <- function_csv <- function(path, gene_col, score_col, to_upper = TRUE) {
  dt <- fread(path)
  stopifnot(gene_col %in% names(dt), score_col %in% names(dt))
  dt <- dt[, .(gene = get(gene_col), score = as.numeric(get(score_col)))]
  dt <- dt[!is.na(gene) & !is.na(score)]
  dt[, gene := str_trim(gene)]
  if (to_upper) dt[, gene := toupper(gene)]
  # 同基因出现多次：取最大（也可改成 mean/median）
  dt <- dt[, .(score = max(score)), by = gene][order(-score)]
  setNames(dt$score, dt$gene)
}

# ---- 读取 GMT（NAME\tDESC\tg1\tg2... 或 NAME\tg1\tg2... 都兼容）----
read_gmt <- function(path, to_upper = TRUE) {
  lines <- readLines(path, warn = FALSE)
  sets <- lapply(lines, function(ln) {
    f <- strsplit(ln, "\t")[[1]]
    f <- f[nzchar(f)]
    if (length(f) < 2) return(NULL)
    nm <- f[1]
    genes <- f[-1]
    # 如果第二列是描述，且后面仍有基因，把第二列空掉：
    # 经验规则：若第二列不全为字母/数字/下划线组合且后面>1基因，则看作描述列
    if (length(genes) > 1 && !grepl("^[A-Za-z0-9._-]+$", genes[1])) {
      genes <- genes[-1]
    }
    genes <- str_trim(genes)
    if (to_upper) genes <- toupper(genes)
    genes[genes != ""]
  })
  nms <- sapply(strsplit(lines, "\t"), `[`, 1)
  names(sets) <- nms
  sets <- Filter(Negate(is.null), sets)
  # 去重
  sets <- lapply(sets, function(v) unique(v[nzchar(v)]))
  sets
}

# ---- 主函数：运行 fgsea ----
run_gsea <- function(stats_named_vector, pathways_list,
                     minSize = 10, maxSize = 5000, nperm = 10000) {
  universe <- names(stats_named_vector)
  # 交集过滤 + 大小限制
  pathways <- lapply(pathways_list, function(g) intersect(g, universe))
  sizes <- sapply(pathways, length)
  keep <- sizes >= minSize & sizes <= maxSize
  pathways <- pathways[keep]
  if (length(pathways) == 0) stop("No pathways left after filtering by universe & size.")
  fg <- fgsea(pathways = pathways,
              stats = stats_named_vector,
              nperm = nperm,
              minSize = minSize,
              maxSize = maxSize)
  # 增加一些实用列
  fg$overlap <- vapply(fg$pathway,
                       function(p) length(intersect(pathways[[p]], universe)),
                       integer(1))
  fg$set_size <- vapply(fg$pathway, function(p) length(pathways[[p]]), integer(1))
  fg[order(padj, pval, -NES), c("pathway","NES","pval","padj","size","set_size","overlap","leadingEdge")]
}

# ---- 可视化单个 regulon 的富集曲线 ----
gsea_plot <- function(pathways_list, stats_named_vector, pathway_name, title = NULL) {
  stopifnot(pathway_name %in% names(pathways_list))
  p <- fgsea::plotEnrichment(pathways_list[[pathway_name]], stats_named_vector)
  if (!is.null(title)) p <- p + ggtitle(title)
  p + theme_minimal(base_size = 12)
}

# =========================  
# 使用示例（替换成你的文件路径与列名）
# =========================
# 1) OpenTargets CSV：需要告诉哪一列是基因/分数
#   常见：gene 符号列类似 "gene" / "geneSymbol"
#        分数列类似 "overallScore" / "association_score" / "mlog10p"/"z"
ot_stats <- function_csv("01_opentargts_monocyte.csv",
                         gene_col  = "symbol",
                         score_col = "gwasCredibleSets")         # 改成你的分数字段

# 2) Regulons GMT：包含多条 RBP_regulon
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_3UTR.gmt")
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_5UTR.gmt")
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_Introns.gmt")

# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 10, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)

  
# 5) 导出结果
data.table::fwrite(fg_res, "fgsea_RBP_regulons_vs_OpenTargets_3UTR.csv")
data.table::fwrite(fg_res, "fgsea_RBP_regulons_vs_OpenTargets_5UTR.csv")
data.table::fwrite(fg_res, "fgsea_RBP_regulons_vs_OpenTargets_Introns.csv")


# 找到 CELF1 regulon 这一行（名字可能是 "CELF1" 或 "CELF1_regulon"）
row <- fg_res[fg_res$pathway %in% c("CELF1", "CELF1_regulon"), ][1, ]

# 取出 leading edge 基因向量
le <- unlist(row$leadingEdge)
length(le); le[1:10]  


# 假设 fg_res 是 run_gsea(...) 的结果；ot_stats 是命名向量(已按降序)
row <- fg_res[fg_res$pathway %in% c("QKI","QKI_regulon"), ][1, ]
row[, c("pathway","NES","size","pval","padj")]

# 取 leading edge 并按 OpenTargets 排名从前到后列出
le <- unlist(row$leadingEdge)
ord <- order(rank(ot_stats))                               # fgsea 用的顺序
le_tab <- data.frame(
  gene = le,
  score = ot_stats[le],
  rank  = match(le, names(ot_stats))
)[order(le_tab$rank), ]
nrow(le_tab)                         # leading edge 基因数
min(le_tab$rank) / length(ot_stats)  # 覆盖起点（占总排名的比例）
max(le_tab$rank) / length(ot_stats)  # 覆盖终点
head(le_tab, 15)                     # 前 15 个 leading edge 基因


# 4) 绘制某个 regulon 的富集曲线（例如 YBX3_regulon）
# 4) 绘制某个 regulon 的富集曲线（例如 YBX3_regulon）
#--3UTR
#gsea_plot(regulons, ot_stats, "CELF1", title = "CELF1_regulon ~ OpenTargets ranks")
#gsea_plot(regulons, ot_stats, "SF3B4", title = "SF3B4_regulon ~ OpenTargets ranks")

#--5UTR
#gsea_plot(regulons, ot_stats, "QKI", title = "QKI_regulon ~ OpenTargets ranks")
#gsea_plot(regulons, ot_stats, "DDX3X", title = "DDX3X_regulon ~ OpenTargets ranks")


suppressPackageStartupMessages({
  library(fgsea)
  library(ggplot2)
})

# 将 data.frame / 向量统一为“命名数值向量”的工具
coerce_stats_vec <- function(stats, gene_col = NULL, score_col = NULL, to_upper = TRUE) {
  if (is.numeric(stats) && !is.null(names(stats))) {
    v <- stats
  } else if (is.data.frame(stats)) {
    stopifnot(!is.null(gene_col), !is.null(score_col))
    df <- stats[, c(gene_col, score_col)]
    colnames(df) <- c("gene","score")
    df <- df[!is.na(df$gene) & !is.na(df$score), ]
    if (to_upper) df$gene <- toupper(trimws(df$gene))
    # 同基因出现多次取最大分数
    df <- aggregate(score ~ gene, df, max)
    v <- df$score; names(v) <- df$gene
  } else {
    stop("`stats` 必须是命名数值向量或 data.frame/tibble。")
  }
  v <- v[is.finite(v)]
  if (to_upper) names(v) <- toupper(names(v))
  # 去重复：同名保留最大
  v <- tapply(v, names(v), max)
  v[order(v, decreasing = TRUE)]
}

# 画 enrichment + 正确标注 LE（红）与非LE（灰）
plot_enrichment_with_le <- function(pathway_name, pathways, stats,
                                    fg_res = NULL,
                                    gene_col = NULL, score_col = NULL,
                                    to_upper = TRUE,
                                    le_col = "red", other_col = "grey40") {
  
  # 1) 排行（命名向量，降序）
  geneList <- coerce_stats_vec(stats, gene_col, score_col, to_upper)
  
  # 2) 取集合（兼容名字可能是 "QKI" 或 "QKI_regulon"）
  set_genes <- pathways[[pathway_name]]
  if (is.null(set_genes)) {
    alt <- sub("_regulon$", "", pathway_name)
    if (!identical(alt, pathway_name)) set_genes <- pathways[[alt]]
  }
  if (is.null(set_genes)) stop("在 pathways 里找不到集合：", pathway_name)
  if (to_upper) set_genes <- toupper(set_genes)
  set_genes <- unique(set_genes)
  
  # 3) overlap 检查
  ov <- intersect(names(geneList), set_genes)
  message(sprintf("Set size = %d, universe = %d, overlap = %d",
                  length(set_genes), length(geneList), length(ov)))
  if (length(ov) == 0) stop("该集合与排行无交集（大小写或ID不一致）。")
  
  # 4) leading edge：优先从 fg_res 拿；拿不到就临时跑一次 fgsea
  le_genes <- NULL
  if (!is.null(fg_res) && "pathway" %in% names(fg_res)) {
    row <- fg_res[fg_res$pathway %in% c(pathway_name, sub("_regulon$", "", pathway_name)), ]
    if (nrow(row) > 0) le_genes <- unlist(row$leadingEdge)
  }
  if (is.null(le_genes)) {
    tmp <- fgsea(pathways = list(tmp = set_genes), stats = geneList, nperm = 10000)
    le_genes <- unlist(tmp$leadingEdge)
  }
  
  # 5) 位置（与 plotEnrichment 保持一致）
  all_idx <- which(names(geneList) %in% set_genes)
  le_idx  <- which(names(geneList) %in% le_genes)
  
  # 6) 绘图
  p <- fgsea::plotEnrichment(set_genes, geneList, ticksSize = 0) +
    ggtitle(paste0(pathway_name, " ~ OpenTargets ranks"))
  p +
    geom_segment(data = data.frame(x = setdiff(all_idx, le_idx)),
                 aes(x = x, xend = x, y = 0, yend = -0.035),
                 inherit.aes = FALSE, size = 0.15, color = other_col) +
    geom_segment(data = data.frame(x = le_idx),
                 aes(x = x, xend = x, y = 0, yend = -0.07),
                 inherit.aes = FALSE, size = 0.35, color = le_col) +
    coord_cartesian(clip = "off")
}

#---5UTR
plot_enrichment_with_le("QKI", regulons, ot_stats, fg_res)
plot_enrichment_with_le("DDX21", regulons, ot_stats, fg_res)
plot_enrichment_with_le("DDX3X", regulons, ot_stats, fg_res)
plot_enrichment_with_le("NIPBL", regulons, ot_stats, fg_res)
plot_enrichment_with_le("TRA2A", regulons, ot_stats, fg_res)


#---3UTR
plot_enrichment_with_le("CELF1", regulons, ot_stats, fg_res)
plot_enrichment_with_le("TAF15", regulons, ot_stats, fg_res)



#---Introns
plot_enrichment_with_le("RBM47", regulons, ot_stats, fg_res)
plot_enrichment_with_le("YBX3", regulons, ot_stats, fg_res)








# 把它们的 OpenTargets 分数与排名一起列出来
library(dplyr)
le_tab <- tibble(
  gene = le,
  score = ot_stats[le],
  rank  = match(le, names(ot_stats))  # ot_stats 已按降序
) |> arrange(rank)
le_tab


# stats 已按降序；pathways 是你的 regulon 列表
ord <- order(rank(ot_stats))                       # fgsea 用的顺序
pos_all <- which(names(ot_stats)[ord] %in% regulons[["CELF1"]])
length(pos_all)    # 应该等于 28
pos_all[1:20]      # 看前20个位置

plot_enrichment_with_le("CELF1", regulons, ot_stats, fg_res)

p <- fgsea::plotEnrichment(regulons[["CELF1"]], ot_stats, ticksSize = 0.15) +
  ggtitle("CELF1_regulon ~ OpenTargets ranks") +
  coord_cartesian(clip = "off")                # 防止被裁剪

# 额外在底部再画一层 rug（每个基因一个刻度；不容易被覆盖）
ticks_df <- data.frame(x = pos_all)
p <- p + ggplot2::geom_rug(data = ticks_df, aes(x = x), sides = "b",
                           inherit.aes = FALSE, alpha = 0.7)
p


# 拿到 leading edge
le <- unlist(fg_res[fg_res$pathway %in% c("CELF1","CELF1_regulon"), ]$leadingEdge)
pos_le <- which(names(ot_stats)[ord] %in% le)
pos_rest <- setdiff(pos_all, pos_le)

p +
  geom_segment(data = data.frame(x = pos_rest),
               aes(x = x, xend = x, y = 0, yend = -0.04),
               inherit.aes = FALSE, size = 0.15, color = "grey40") +
  geom_segment(data = data.frame(x = pos_le),
               aes(x = x, xend = x, y = 0, yend = -0.08),
               inherit.aes = FALSE, size = 0.3, color = "red")


# 导出更宽的图片以减少重叠
# ggsave("CELF1_gsea.png", p, width = 12, height = 5, dpi = 300)



#------Figure_2e

##_--------dotplot for regulon-opentargets enrichment
# ---- 数据 ----
txt <- "Regulon	NES	logp	padj
QKI_5UTR	1.307591657	5.000004343	0.000559994
DDX3X_5UTR	1.238209431	3.366535887	0.01203988
DDX21_5UTR	1.274739439	3.148745994	0.013253201
NIPBL_5UTR	1.255780212	2.567035052	0.025386413
TRA2A_5UTR	1.229844926	2.565435439	0.025386413
CELF1_3UTR	1.426940745	2.086416332	0.17544036
RBM47_Introns	1.238423435	2.037161662	0.183598164
TRA2A_3UTR	1.298532184	1.866146218	0.17544036
TAF15_3UTR	1.294268227	1.835949051	0.17544036
IGF2BP1_5UTR	1.424708484	1.478071846	0.155215445
DDX3X_3UTR	1.18704609	1.474182391	0.18910922
FAM120A_5UTR	1.256236336	1.296738167	0.200415774
YBX3_Introns	1.210306842	1.258218304	0.5518
PUM2_5UTR	1.127683633	1.234559325	0.200415774
WTAP_5UTR	1.250467089	1.202892334	0.200415774
DDX21_3UTR	1.141993404	1.078682099	0.32819231
RPRD2_5UTR	1.267089543	0.923869202	0.290128899
EIF4A1_5UTR	1.091410405	0.877886807	0.296729833
RYBP_3UTR	1.094149572	0.620517529	0.51800091
HNRNPH2_5UTR	1.130069427	0.560225305	0.468543335
YTHDF3_5UTR	1.089278448	0.558924323	0.468543335
ZNF148_Introns	1.047549317	0.515320216	0.614464295
PABPC1_Introns	1.033660783	0.453683819	0.670126632
HNRNPC_Introns	0.992828913	0.271212429	0.715745274
ILF3_3UTR	0.8350837	0.110167746	0.946277254
XRN2_Introns	0.914503258	0.07780326	0.926769111"

df <- read.table(text = txt, header = TRUE, sep = "\t", check.names = FALSE)

# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)

# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 4) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 3.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 12)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)





