
#2025-12-16

# Rheumatoid Arthritis (RA, N = 14,361 cases and 42,923 controls)

#based on all OpenTarget association genes

setwd("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/04_scRBP_trs_PBMC_eight_autoimmune_diseases/02_OpenTarget_results_based_on_all_genes")

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


#----RA
# =========================  
# 使用示例（替换成你的文件路径与列名）
# =========================
# 1) OpenTargets CSV：需要告诉哪一列是基因/分数
#   常见：gene 符号列类似 "gene" / "geneSymbol"
#        分数列类似 "overallScore" / "association_score" / "mlog10p"/"z"
ot_stats <- function_csv("01_opentargts_RA_all.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # 改成你的分数字段



# 2) Regulons GMT：包含多条 RBP_regulon
#3UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_3UTR.gmt")

# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) 导出结果
data.table::fwrite(fg_res, "RA_fgsea_RBP_regulons_vs_OpenTargets_3UTR_ALLgenes.csv")


###_-------GSEA富集图:
#load function: plot_enrichment_with_le.R
source("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/04_scRBP_trs_PBMC_eight_autoimmune_diseases/02_OpenTarget_results_based_on_all_genes/plot_enrichment_with_le.R")

regulon_rbp <- "DDX3X"
# 单张：
p <- plot_enrichment_with_le(regulon_rbp, regulons, ot_stats, fg_res = fg_res)
print(p)


regulon_rbp <- "RPS11"
# 单张：
p <- plot_enrichment_with_le(regulon_rbp, regulons, ot_stats, fg_res = fg_res)
print(p)





#5UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_5UTR.gmt")

# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "RA_fgsea_RBP_regulons_vs_OpenTargets_5UTR_ALLgenes.csv")

###_-------GSEA富集图:
#load function: plot_enrichment_with_le.R
source("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/04_scRBP_trs_PBMC_eight_autoimmune_diseases/02_OpenTarget_results_based_on_all_genes/plot_enrichment_with_le.R")

regulon_rbp <- "SF3B1"
# 单张：
p <- plot_enrichment_with_le(regulon_rbp, regulons, ot_stats, fg_res = fg_res)
print(p)




#CDS
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_CDS.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "RA_fgsea_RBP_regulons_vs_OpenTargets_CDS_ALLgenes.csv")


#Introns
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_Introns.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "RA_fgsea_RBP_regulons_vs_OpenTargets_Introns_ALLgenes.csv")




####-----MS

ot_stats <- function_csv("01_opentargts_MS_all.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # 改成你的分数字段

# 2) Regulons GMT：包含多条 RBP_regulon
#3UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_3UTR.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) 导出结果
data.table::fwrite(fg_res, "MS_fgsea_RBP_regulons_vs_OpenTargets_3UTR_ALLgenes.csv")

#5UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_5UTR.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "MS_fgsea_RBP_regulons_vs_OpenTargets_5UTR_ALLgenes.csv")

#CDS
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_CDS.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "MS_fgsea_RBP_regulons_vs_OpenTargets_CDS_ALLgenes.csv")


#Introns
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_Introns.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "MS_fgsea_RBP_regulons_vs_OpenTargets_Introns_ALLgenes.csv")




####-----IBD_ALLgenes

ot_stats <- function_csv("01_opentargts_IBD_all.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # 改成你的分数字段

# 2) Regulons GMT：包含多条 RBP_regulon
#3UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_3UTR.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) 导出结果
data.table::fwrite(fg_res, "IBD_fgsea_RBP_regulons_vs_OpenTargets_3UTR_ALLgenes.csv")

#5UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_5UTR.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "IBD_fgsea_RBP_regulons_vs_OpenTargets_5UTR_ALLgenes.csv")

#CDS
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_CDS.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "IBD_fgsea_RBP_regulons_vs_OpenTargets_CDS_ALLgenes.csv")


#Introns
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_Introns.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "IBD_fgsea_RBP_regulons_vs_OpenTargets_Introns_ALLgenes.csv")




####-----UC_ALLgenes

ot_stats <- function_csv("01_opentargts_UC_all.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # 改成你的分数字段

# 2) Regulons GMT：包含多条 RBP_regulon
#3UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_3UTR.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) 导出结果
data.table::fwrite(fg_res, "UC_fgsea_RBP_regulons_vs_OpenTargets_3UTR_ALLgenes.csv")


###_-------GSEA富集图:
#load function: plot_enrichment_with_le.R
source("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/04_scRBP_trs_PBMC_eight_autoimmune_diseases/02_OpenTarget_results_based_on_all_genes/plot_enrichment_with_le.R")

regulon_rbp <- "ZFP36"
# 单张：
p <- plot_enrichment_with_le(regulon_rbp, regulons, ot_stats, fg_res = fg_res)
print(p)


#5UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_5UTR.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "UC_fgsea_RBP_regulons_vs_OpenTargets_5UTR_ALLgenes.csv")

#CDS
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_CDS.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "UC_fgsea_RBP_regulons_vs_OpenTargets_CDS_ALLgenes.csv")






#Introns
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_Introns.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "UC_fgsea_RBP_regulons_vs_OpenTargets_Introns_ALLgenes.csv")




####-----SLE_ALLgenes

ot_stats <- function_csv("01_opentargts_SLE_all.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # 改成你的分数字段

# 2) Regulons GMT：包含多条 RBP_regulon
#3UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_3UTR.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) 导出结果
data.table::fwrite(fg_res, "SLE_fgsea_RBP_regulons_vs_OpenTargets_3UTR_ALLgenes.csv")

#5UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_5UTR.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "SLE_fgsea_RBP_regulons_vs_OpenTargets_5UTR_ALLgenes.csv")

#CDS
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_CDS.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "SLE_fgsea_RBP_regulons_vs_OpenTargets_CDS_ALLgenes.csv")


#Introns
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_Introns.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "SLE_fgsea_RBP_regulons_vs_OpenTargets_Introns_ALLgenes.csv")






####-----T1D_ALLgenes

ot_stats <- function_csv("01_opentargts_T1D_all.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # 改成你的分数字段

# 2) Regulons GMT：包含多条 RBP_regulon
#3UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_3UTR.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) 导出结果
data.table::fwrite(fg_res, "T1D_fgsea_RBP_regulons_vs_OpenTargets_3UTR_ALLgenes.csv")

#5UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_5UTR.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "T1D_fgsea_RBP_regulons_vs_OpenTargets_5UTR_ALLgenes.csv")

#CDS
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_CDS.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "T1D_fgsea_RBP_regulons_vs_OpenTargets_CDS_ALLgenes.csv")


#Introns
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_Introns.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "T1D_fgsea_RBP_regulons_vs_OpenTargets_Introns_ALLgenes.csv")





####-----PBC_ALLgenes

ot_stats <- function_csv("01_opentargts_PBC_all.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # 改成你的分数字段

# 2) Regulons GMT：包含多条 RBP_regulon
#3UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_3UTR.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) 导出结果
data.table::fwrite(fg_res, "PBC_fgsea_RBP_regulons_vs_OpenTargets_3UTR_ALLgenes.csv")

#5UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_5UTR.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "PBC_fgsea_RBP_regulons_vs_OpenTargets_5UTR_ALLgenes.csv")

#CDS
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_CDS.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "PBC_fgsea_RBP_regulons_vs_OpenTargets_CDS_ALLgenes.csv")


#Introns
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_Introns.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "PBC_fgsea_RBP_regulons_vs_OpenTargets_Introns_ALLgenes.csv")




####-----CD_ALLgenes

ot_stats <- function_csv("01_opentargts_CD_all.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # 改成你的分数字段

# 2) Regulons GMT：包含多条 RBP_regulon
#3UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_3UTR.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) 导出结果
data.table::fwrite(fg_res, "CD_fgsea_RBP_regulons_vs_OpenTargets_3UTR_ALLgenes.csv")

#5UTR
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_5UTR.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "CD_fgsea_RBP_regulons_vs_OpenTargets_5UTR_ALLgenes.csv")

#CDS
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_CDS.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "CD_fgsea_RBP_regulons_vs_OpenTargets_CDS_ALLgenes.csv")


#Introns
regulons <- read_gmt("../03_50runs_regulons_4RNA_Domains/scRBP_prune_matched_results_min10genes.symbol_default_Introns.gmt")
# 3) 运行 GSEA（建议 nperm 至少 10000；正式结果 100000）
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
data.table::fwrite(fg_res, "CD_fgsea_RBP_regulons_vs_OpenTargets_Introns_ALLgenes.csv")






###_-------GSEA富集图:
###_-------GSEA富集图:
###_-------GSEA富集图:
###_-------GSEA富集图:
#load function: plot_enrichment_with_le.R
source("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/04_scRBP_trs_PBMC_eight_autoimmune_diseases/02_OpenTarget_results_based_on_all_genes/plot_enrichment_with_le.R")

regulon_rbp <- "DDX3X"
# 单张：
p <- plot_enrichment_with_le(regulon_rbp, regulons, ot_stats, fg_res = fg_res)
print(p)

# 多张（每页一个图，保存 PDF）：
# plot_gsea_panels(
#   pathway_names = c("QKI", "DDX3X", "CELF1", "SF3B4"),
#   pathways      = regulons,
#   stats         = ot_stats,
#   fg_res        = fg_res,
#   out_pdf       = "GSEA_regulons_OpenTargets.pdf",
#   ncol          = 1
# )

# 多张拼版（需要 patchwork）：
# plot_gsea_panels(
#   pathway_names = c("QKI", "DDX3X", "CELF1", "SF3B4"),
#   pathways      = regulons,
#   stats         = ot_stats,
#   fg_res        = fg_res,
#   out_pdf       = "GSEA_regulons_OpenTargets_grid.pdf",
#   ncol          = 2
# )






##_--------dotplot for regulon-opentargets enrichment
##_--------dotplot for regulon-opentargets enrichment
##_--------dotplot for regulon-opentargets enrichment
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

#df <- read.csv(text = txt, header = TRUE, sep = "\t", check.names = FALSE)

 


########------------IBD
df <- read.csv("./02_opentargts_IBD_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon
keep_regs <- c(
  "PCBP1_3UTR",
  "DDX3X_Introns",
  "RPS11_5UTR",
  "RPS11_Introns",
  "EIF4G2_5UTR",
  "SF3B1_5UTR",
  "MAP1LC3B_CDS",
  "DDX3X_CDS",
  "ZFP36_3UTR",
  "HNRNPC_Introns",
  "RPS15A_5UTR",
  "LIN28B_5UTR",
  "HNRNPF_5UTR",
  "DDX3X_3UTR"
)

# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)







##-----------RA
df <- read.csv("./02_opentargts_RA_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon
keep_regs <- c(
  "RPS11_Introns",
  "SF3B1_5UTR",
  "EIF4G2_5UTR",
  "DDX3X_3UTR",
  "KHDRBS2_5UTR",
  "DDX3X_CDS",
  "RPS11_3UTR",
  "RPS11_5UTR"
)

# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)




########------------MS
df <- read.csv("./02_opentargts_MS_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)

# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon
keep_regs <- c(
  "HNRNPC_Introns",
  "KHDRBS2_5UTR",
  "PCBP1_3UTR",
  "PTBP1_CDS",
  "YBX1_3UTR",
  "RPS11_5UTR",
  "EIF4G2_5UTR",
  "DDX3X_3UTR",
  "SF3B1_5UTR",
  "RPS11_Introns"
)

# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)





########------------UC
df <- read.csv("./02_opentargts_UC_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon
keep_regs <- c(
  "DDX3X_Introns",
  "EIF4G2_5UTR",
  "SF3B1_5UTR",
  "RPS11_5UTR",
  "DDX3X_3UTR",
  "ZFP36L1_3UTR",
  "RPS11_Introns",
  "PCBP1_3UTR",
  "SBDS_3UTR",
  "ZFP36_3UTR",
  "HNRNPC_Introns",
  "MAP1LC3B_CDS",
  "NIPBL_Introns",
  "DDX3X_CDS",
  "RPS15A_5UTR"
)

# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)





########-----------SLE
df <- read.csv("./02_opentargts_SLE_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon
keep_regs <- c(
  "WTAP_CDS",
  "KHDRBS2_5UTR"
)

# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)





########-----------T1D
df <- read.csv("./02_opentargts_T1D_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon
keep_regs <- c(
  "PABPC1_Introns",
  "HNRNPM_5UTR",
  "RPS15A_5UTR",
  "HNRNPC_Introns",
  "HNRNPF_5UTR",
  "RPS11_5UTR",
  "APOBEC3C_5UTR",
  "HNRNPUL1_Introns",
  "KHDRBS2_5UTR"
)

# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)







########-----------PBC
df <- read.csv("./02_opentargts_PBC_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon
keep_regs <- c(
  "RPS11_Introns",
  "DDX3X_CDS",
  "MAP1LC3B_5UTR"
)

# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)







########----------CD
df <- read.csv("./02_opentargts_CD_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon
keep_regs <- c(
  "PCBP1_3UTR",
  "EIF4G2_5UTR",
  "RPS11_5UTR",
  "SF3B1_5UTR",
  "HNRNPC_Introns",
  "RPS11_Introns",
  "DDX3X_CDS",
  "HNRNPF_5UTR",
  "DDX3X_Introns",
  "ZFP36L1_3UTR",
  "MAP1LC3B_CDS",
  "TOP1_5UTR",
  "SBDS_3UTR",
  "APOBEC3C_5UTR",
  "DDX3X_3UTR"
)

# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)



