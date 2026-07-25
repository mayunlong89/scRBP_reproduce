#2025-12-29

# =========================  
# 9 blood cell traits — GSEA: OpenTargets CSV × RBP regulons (GMT)
# =========================

# Dependencies
# BiocManager::install("fgsea")

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
})

# ============================================================
# Helper functions
# ============================================================

# Read OpenTargets CSV: trim, deduplicate (keep max score), sort descending
read_opentargets_csv <- function(path, gene_col, score_col, to_upper = TRUE) {
  dt <- fread(path)
  stopifnot(gene_col %in% names(dt), score_col %in% names(dt))
  dt <- dt[, .(gene = get(gene_col), score = as.numeric(get(score_col)))]
  dt <- dt[!is.na(gene) & !is.na(score)]
  dt[, gene := str_trim(gene)]
  if (to_upper) dt[, gene := toupper(gene)]
  dt <- dt[, .(score = max(score)), by = gene][order(-score)]
  setNames(dt$score, dt$gene)
}

# Read GMT (compatible with NAME\tDESC\tgenes or NAME\tgenes)
read_gmt <- function(path, to_upper = TRUE) {
  lines <- readLines(path, warn = FALSE)
  sets <- lapply(lines, function(ln) {
    f <- strsplit(ln, "\t")[[1]]
    f <- f[nzchar(f)]
    if (length(f) < 2) return(NULL)
    genes <- f[-1]
    if (length(genes) > 1 && !grepl("^[A-Za-z0-9._-]+$", genes[1])) {
      genes <- genes[-1]
    }
    genes <- str_trim(genes)
    if (to_upper) genes <- toupper(genes)
    genes[genes != ""]
  })
  names(sets) <- sapply(strsplit(lines, "\t"), `[`, 1)
  sets <- Filter(Negate(is.null), sets)
  lapply(sets, function(v) unique(v[nzchar(v)]))
}

# Run fgsea
run_gsea <- function(stats_named_vector, pathways_list,
                     minSize = 10, maxSize = 5000, nperm = 10000) {
  universe <- names(stats_named_vector)
  pathways <- lapply(pathways_list, function(g) intersect(g, universe))
  sizes <- sapply(pathways, length)
  pathways <- pathways[sizes >= minSize & sizes <= maxSize]
  if (length(pathways) == 0) stop("No pathways left after filtering.")
  fg <- fgsea(pathways = pathways, stats = stats_named_vector,
              nperm = nperm, minSize = minSize, maxSize = maxSize)
  fg$overlap <- vapply(fg$pathway,
                       function(p) length(intersect(pathways[[p]], universe)), integer(1))
  fg$set_size <- vapply(fg$pathway, function(p) length(pathways[[p]]), integer(1))
  fg[order(padj, pval, -NES),
     c("pathway","NES","pval","padj","size","set_size","overlap","leadingEdge")]
}

# ============================================================
# Configuration: 9 traits (folder, OT csv, short name)
# ============================================================
traits <- list(
  list(dir = "00-OpenTargets_lymphocyte_count",   ot_csv = "01_opentargts_lymp_count.csv",   tag = "lymph_count"),
  list(dir = "00-OpenTargets_lymphocyte_percent", ot_csv = "01_opentargts_lymp_percent.csv",  tag = "lymph_percent"),
  list(dir = "01_MCV",                            ot_csv = "01_opentargts_MCV.csv",           tag = "MCV"),
  list(dir = "02_neutro",                         ot_csv = "01_opentargts_neutr_count.csv",   tag = "neutr_count"),
  list(dir = "03_WBC",                            ot_csv = "01_opentargts_WBC.csv",           tag = "WBC"),
  list(dir = "04_eosinophil",                     ot_csv = "01_opentargts_eosin_count.csv",   tag = "eosin_count"),
  list(dir = "05_basophil",                       ot_csv = "01_opentargts_baso_count.csv",    tag = "baso_count"),
  list(dir = "06_MCHC",                           ot_csv = "01_opentargts_MCHC.csv",          tag = "MCHC"),
  list(dir = "07_HL",                             ot_csv = "01_opentargts_HL.csv",            tag = "HL")
)

# 4 mRNA regions and their GMT files
gmt_dir <- "../00-6.5K_100runs_regulons-GO-KEGG"
regions <- c("3UTR", "5UTR", "CDS", "Introns")
gmt_files <- setNames(
  file.path(gmt_dir, paste0("scRBP_prune_matched_results_min10genes.symbol_default_6.5K_", regions, ".gmt")),
  regions
)

# GSEA parameters
NPERM    <- 100000
MIN_SIZE <- 5
MAX_SIZE <- 5000

# Regulons to highlight in dotplots
keep_regs <- c(
  "DDX21_5'UTR","DDX3X_5'UTR","QKI_5'UTR","YTHDF3_5'UTR","IGF2BP2_5'UTR",
  "AKAP1_5'UTR","TNRC6C_5'UTR","HNRNPC_5'UTR","SRSF5_5'UTR","DDX21_3'UTR",
  "YTHDF3_3'UTR","RYBP_3'UTR","CPEB4_CDS","RPS11_CDS","RPS15A_CDS",
  "VIM_CDS","NIPBL_CDS","FUBP3_CDS","CENPC_CDS","IFIT2_Intron",
  "PRPF8_Intron","YBX3_Intron","DDX6_Intron","SRSF1_Intron","ELAVL1_Intron",
  "DDX3X_3'UTR","IGF2BP2_3'UTR","CELF1_3'UTR","TRA2A_5'UTR","ZC3HAV1_5'UTR",
  "PABPC1_Intron","EIF3D_3'UTR","RPS19_5'UTR","ZRANB2_5'UTR","DICER1_CDS",
  "G3BP1_Intron","HNRNPA1_Intron","SRSF2_Intron","OAS1_Intron","FAM120A_3'UTR",
  "SRP14_3'UTR","NCBP3_3'UTR","RBM47_Intron","DEK_Intron","G3BP2_5'UTR",
  "HNRNPU_CDS","SAMD4A_Intron","MTDH_Intron","RPS15A_3'UTR","RPS19_3'UTR",
  "TRA2A_3'UTR","SERBP1_5'UTR","EIF4A1_5'UTR","KHSRP_5'UTR","RPRD2_5UTR",
  "FAM120A_5'UTR","ZRANB2_Intron","XRN2_Intron","ZRANB2_3'UTR","TAF15_3'UTR",
  "IGF2BP2_CDS","YBX3_CDS","RTCB_CDS","METAP2_3'UTR","RPS3A_Intron",
  "YWHAZ_CDS","NOP58_CDS","CHD7_CDS","TNRC6A_CDS","FMR1_3'UTR"
)

# Color palette for dotplots
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)

# ============================================================
# Main loop: iterate over traits × regions
# ============================================================
for (tr in traits) {
  cat("\n===== Processing:", tr$tag, "=====\n")

  # 1) Read OpenTargets scores
  ot_stats <- read_opentargets_csv(
    file.path(tr$dir, tr$ot_csv),
    gene_col = "symbol", score_col = "globalScore"
  )

  # 2) Run GSEA for each mRNA region
  for (reg in names(gmt_files)) {
    cat("  Region:", reg, "\n")
    regulons <- read_gmt(gmt_files[[reg]])
    fg_res   <- run_gsea(ot_stats, regulons,
                         minSize = MIN_SIZE, maxSize = MAX_SIZE, nperm = NPERM)
    fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
    fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

    outfile <- file.path(tr$dir,
                         paste0("fgsea_RBP_regulons_vs_OpenTargets_", reg, "_", tr$tag, ".csv"))
    fwrite(fg_res, outfile)
  }

  # 3) Dotplot (assumes merged dotplot CSV already exists)
  dotplot_csv <- file.path(tr$dir, "01_enrichment_GSEA_dotplot.csv")
  if (file.exists(dotplot_csv)) {
    df <- read.csv(dotplot_csv)
    df$logp    <- -log10(df$pval)
    df$NES     <- abs(df$NES)
    df$Regulon <- trimws(df$Regulon)
    df$label   <- ifelse(df$Regulon %in% keep_regs & df$pval < 0.05, df$Regulon, NA)

    p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
      geom_point(size = 2) +
      scale_color_gradientn(colors = pal, name = "NES") +
      geom_text_repel(aes(label = label),
                      max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                      min.segment.length = 0, seed = 42, size = 2.5) +
      labs(x = expression(-log[10](italic(P))), y = "NES",
           title = paste0("Regulon enrichment — ", tr$tag)) +
      theme_classic(base_size = 10)

    ggsave(file.path(tr$dir, paste0("GSEA_dotplot_", tr$tag, ".pdf")),
           p, width = 8, height = 6)
    print(p)
  }
}

cat("\nAll done!\n")



# =========================  

# Previously, we sparately run the code according to each blood cell trait.
                        
# =========================  
#2025-12-29

# =========================  
# 9 blood cell traits
# =========================


# =========================
# GSEA: OpenTargets CSV  ×  RBP regulons (GMT)
# =========================
# Dependencies
#BiocManager::install("fgsea")

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(stringr)
  library(ggplot2)
})

# ---- Read OpenTargets CSV: specify gene column and score field ----
# Will: trim whitespace, remove NA, take max score for duplicate genes, sort by score descending
read_opentargets_csv <- function_csv <- function(path, gene_col, score_col, to_upper = TRUE) {
  dt <- fread(path)
  stopifnot(gene_col %in% names(dt), score_col %in% names(dt))
  dt <- dt[, .(gene = get(gene_col), score = as.numeric(get(score_col)))]
  dt <- dt[!is.na(gene) & !is.na(score)]
  dt[, gene := str_trim(gene)]
  if (to_upper) dt[, gene := toupper(gene)]
  # If a gene appears multiple times: take max (can also use mean/median)
  dt <- dt[, .(score = max(score)), by = gene][order(-score)]
  setNames(dt$score, dt$gene)
}

# ---- Read GMT (compatible with NAME\tDESC\tg1\tg2... or NAME\tg1\tg2...) ----
read_gmt <- function(path, to_upper = TRUE) {
  lines <- readLines(path, warn = FALSE)
  sets <- lapply(lines, function(ln) {
    f <- strsplit(ln, "\t")[[1]]
    f <- f[nzchar(f)]
    if (length(f) < 2) return(NULL)
    nm <- f[1]
    genes <- f[-1]
    # If the second field is a description and more genes follow, skip it:
    # Heuristic: if the second field is not purely alphanumeric and >1 gene remains, treat it as a description column
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
  # Deduplicate
  sets <- lapply(sets, function(v) unique(v[nzchar(v)]))
  sets
}

# ---- Main function: run fgsea ----
run_gsea <- function(stats_named_vector, pathways_list,
                     minSize = 10, maxSize = 5000, nperm = 10000) {
  universe <- names(stats_named_vector)
  # Intersection filtering + size constraints
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
  # Add some useful columns
  fg$overlap <- vapply(fg$pathway,
                       function(p) length(intersect(pathways[[p]], universe)),
                       integer(1))
  fg$set_size <- vapply(fg$pathway, function(p) length(pathways[[p]]), integer(1))
  fg[order(padj, pval, -NES), c("pathway","NES","pval","padj","size","set_size","overlap","leadingEdge")]
}

# ---- Visualize enrichment curve for a single regulon ----
gsea_plot <- function(pathways_list, stats_named_vector, pathway_name, title = NULL) {
  stopifnot(pathway_name %in% names(pathways_list))
  p <- fgsea::plotEnrichment(pathways_list[[pathway_name]], stats_named_vector)
  if (!is.null(title)) p <- p + ggtitle(title)
  p + theme_minimal(base_size = 12)
}



#########################################################################
#########################################################################
#1) lymphocyte count
#########################################################################
#########################################################################

# 1) OpenTargets CSV: specify which column is gene/score
#   Common: gene symbol column like "gene" / "geneSymbol"
#           score column like "overallScore" / "association_score" / "mlog10p"/"z"
ot_stats <- function_csv("00-OpenTargets_lymphocyte_count/01_opentargts_lymp_count.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # Change to your score field


#--3'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_3UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 10, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./00-OpenTargets_lymphocyte_count/fgsea_RBP_regulons_vs_OpenTargets_3UTR_lymph_count.csv")




#--5'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_5UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 10, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./00-OpenTargets_lymphocyte_count/fgsea_RBP_regulons_vs_OpenTargets_5UTR_lymph_count.csv")



#--CDS
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_CDS.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 10, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./00-OpenTargets_lymphocyte_count/fgsea_RBP_regulons_vs_OpenTargets_CDS_lymph_count.csv")




#--Introns
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_Introns.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 10, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./00-OpenTargets_lymphocyte_count/fgsea_RBP_regulons_vs_OpenTargets_Introns_lymph_count.csv")



#---Plot---GSEA---Dotplot---for lymphocyte count
##_--------dotplot for regulon-opentargets enrichment


##_-----In this version, we merge results from all 4 mRNA regions, then display enrichment and label significant results
########------------monocyte count all L2G genes
df <- read.csv("./00-OpenTargets_lymphocyte_count/01_enrichment_GSEA_dotplot.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# Color palette: smaller padj = darker color
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# Labeling: only add labels for points with FDR < 0.05 (adjust threshold as needed)
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# Only label specified regulons
keep_regs <- c(
  "DDX21_5'UTR",
  "DDX3X_5'UTR",
  "QKI_5'UTR",
  "YTHDF3_5'UTR",
  "IGF2BP2_5'UTR",
  "AKAP1_5'UTR",
  "TNRC6C_5'UTR",
  "HNRNPC_5'UTR",
  "SRSF5_5'UTR",
  "DDX21_3'UTR",
  "YTHDF3_3'UTR",
  "RYBP_3'UTR",
  "CPEB4_CDS",
  "RPS11_CDS",
  "RPS15A_CDS",
  "VIM_CDS",
  "NIPBL_CDS",
  "FUBP3_CDS",
  "CENPC_CDS",
  "IFIT2_Intron",
  "PRPF8_Intron",
  "YBX3_Intron",
  "DDX6_Intron",
  "SRSF1_Intron",
  "ELAVL1_Intron",
  "DDX3X_3'UTR",
  "IGF2BP2_3'UTR",
  "CELF1_3'UTR",
  "TRA2A_5'UTR",
  "ZC3HAV1_5'UTR",
  "PABPC1_Intron",
  "EIF3D_3'UTR",
  "RPS19_5'UTR",
  "ZRANB2_5'UTR",
  "DICER1_CDS",
  "G3BP1_Intron",
  "HNRNPA1_Intron",
  "SRSF2_Intron",
  "OAS1_Intron",
  "FAM120A_3'UTR",
  "SRP14_3'UTR",
  "NCBP3_3'UTR",
  "RBM47_Intron",
  "DEK_Intron",
  "G3BP2_5'UTR",
  "HNRNPU_CDS",
  "SAMD4A_Intron",
  "MTDH_Intron",
  "RPS15A_3'UTR",
  "RPS19_3'UTR",
  "TRA2A_3'UTR",
  "SERBP1_5'UTR",
  "EIF4A1_5'UTR",
  "KHSRP_5'UTR",
  "RPRD2_5UTR",
  "FAM120A_5'UTR",
  "ZRANB2_Intron", 
  "XRN2_Intron",
  "ZRANB2_3'UTR",
  "TAF15_3'UTR",
  "IGF2BP2_CDS",
  "YBX3_CDS",
  "RTCB_CDS",
  "METAP2_3'UTR",
  "RPS3A_Intron",
  "YWHAZ_CDS",
  "NOP58_CDS",
  "CHD7_CDS",
  "TNRC6A_CDS",
  "FMR1_3'UTR"
)

# Optional: trim potential whitespace/invisible characters to avoid matching issues
df$Regulon <- trimws(df$Regulon)

#df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)
df$label <- ifelse(df$Regulon %in% keep_regs & df$pval < 0.05, df$Regulon, NA)

# ---- Plot ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # Reverse color scale: smaller padj = darker (use rev(pal))
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p













#########################################################################
#########################################################################
#2) lymphocyte percent
#########################################################################
#########################################################################

# 1) OpenTargets CSV: specify which column is gene/score
#   Common: gene symbol column like "gene" / "geneSymbol"
#           score column like "overallScore" / "association_score" / "mlog10p"/"z"
ot_stats <- function_csv("00-OpenTargets_lymphocyte_percent/01_opentargts_lymp_percent.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # Change to your score field


#--3'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_3UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./00-OpenTargets_lymphocyte_percent/fgsea_RBP_regulons_vs_OpenTargets_3UTR_lymph_percent.csv")




#--5'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_5UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./00-OpenTargets_lymphocyte_percent/fgsea_RBP_regulons_vs_OpenTargets_5UTR_lymph_percent.csv")



#--CDS
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_CDS.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./00-OpenTargets_lymphocyte_percent/fgsea_RBP_regulons_vs_OpenTargets_CDS_lymph_percent.csv")




#--Introns
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_Introns.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./00-OpenTargets_lymphocyte_percent/fgsea_RBP_regulons_vs_OpenTargets_Introns_lymph_percent.csv")



#---Plot---GSEA---Dotplot---for lymphocyte percent
##_--------dotplot for regulon-opentargets enrichment


##_-----In this version, we merge results from all 4 mRNA regions, then display enrichment and label significant results
########------------monocyte count all L2G genes
df <- read.csv("./00-OpenTargets_lymphocyte_percent/01_enrichment_GSEA_dotplot.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# Color palette: smaller padj = darker color
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# Labeling: only add labels for points with FDR < 0.05 (adjust threshold as needed)
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# Only label specified regulons
keep_regs <- c(
  "DDX21_5'UTR",
  "DDX3X_5'UTR",
  "QKI_5'UTR",
  "YTHDF3_5'UTR",
  "IGF2BP2_5'UTR",
  "AKAP1_5'UTR",
  "TNRC6C_5'UTR",
  "HNRNPC_5'UTR",
  "SRSF5_5'UTR",
  "DDX21_3'UTR",
  "YTHDF3_3'UTR",
  "RYBP_3'UTR",
  "CPEB4_CDS",
  "RPS11_CDS",
  "RPS15A_CDS",
  "VIM_CDS",
  "NIPBL_CDS",
  "FUBP3_CDS",
  "CENPC_CDS",
  "IFIT2_Intron",
  "PRPF8_Intron",
  "YBX3_Intron",
  "DDX6_Intron",
  "SRSF1_Intron",
  "ELAVL1_Intron",
  "DDX3X_3'UTR",
  "IGF2BP2_3'UTR",
  "CELF1_3'UTR",
  "TRA2A_5'UTR",
  "ZC3HAV1_5'UTR",
  "PABPC1_Intron",
  "EIF3D_3'UTR",
  "RPS19_5'UTR",
  "ZRANB2_5'UTR",
  "DICER1_CDS",
  "G3BP1_Intron",
  "HNRNPA1_Intron",
  "SRSF2_Intron",
  "OAS1_Intron",
  "FAM120A_3'UTR",
  "SRP14_3'UTR",
  "NCBP3_3'UTR",
  "RBM47_Intron",
  "DEK_Intron",
  "G3BP2_5'UTR",
  "HNRNPU_CDS",
  "SAMD4A_Intron",
  "MTDH_Intron",
  "RPS15A_3'UTR",
  "RPS19_3'UTR",
  "TRA2A_3'UTR",
  "SERBP1_5'UTR",
  "EIF4A1_5'UTR",
  "KHSRP_5'UTR",
  "RPRD2_5UTR",
  "FAM120A_5'UTR",
  "ZRANB2_Intron", 
  "XRN2_Intron",
  "ZRANB2_3'UTR",
  "TAF15_3'UTR",
  "IGF2BP2_CDS",
  "YBX3_CDS",
  "RTCB_CDS",
  "METAP2_3'UTR",
  "RPS3A_Intron",
  "YWHAZ_CDS",
  "NOP58_CDS",
  "CHD7_CDS",
  "TNRC6A_CDS",
  "FMR1_3'UTR"
)

# Optional: trim potential whitespace/invisible characters to avoid matching issues
df$Regulon <- trimws(df$Regulon)

#df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)
df$label <- ifelse(df$Regulon %in% keep_regs & df$pval < 0.05, df$Regulon, NA)

# ---- Plot ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # Reverse color scale: smaller padj = darker (use rev(pal))
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p









#########################################################################
#########################################################################
#3) MCV
#########################################################################
#########################################################################

# 1) OpenTargets CSV: specify which column is gene/score
#   Common: gene symbol column like "gene" / "geneSymbol"
#           score column like "overallScore" / "association_score" / "mlog10p"/"z"
ot_stats <- function_csv("01_MCV/01_opentargts_MCV.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # Change to your score field


#--3'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_3UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./01_MCV/fgsea_RBP_regulons_vs_OpenTargets_3UTR_MCV.csv")




#--5'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_5UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./01_MCV/fgsea_RBP_regulons_vs_OpenTargets_5UTR_MCV.csv")



#--CDS
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_CDS.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./01_MCV/fgsea_RBP_regulons_vs_OpenTargets_CDS_MCV.csv")




#--Introns
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_Introns.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./01_MCV/fgsea_RBP_regulons_vs_OpenTargets_Introns_MCV.csv")



#---Plot---GSEA---Dotplot---for lymphocyte percent
##_--------dotplot for regulon-opentargets enrichment


##_-----In this version, we merge results from all 4 mRNA regions, then display enrichment and label significant results
########------------monocyte count all L2G genes
df <- read.csv("./01_MCV/01_enrichment_GSEA_dotplot.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# Color palette: smaller padj = darker color
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# Labeling: only add labels for points with FDR < 0.05 (adjust threshold as needed)
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# Only label specified regulons
keep_regs <- c(
  "DDX21_5'UTR",
  "DDX3X_5'UTR",
  "QKI_5'UTR",
  "YTHDF3_5'UTR",
  "IGF2BP2_5'UTR",
  "AKAP1_5'UTR",
  "TNRC6C_5'UTR",
  "HNRNPC_5'UTR",
  "SRSF5_5'UTR",
  "DDX21_3'UTR",
  "YTHDF3_3'UTR",
  "RYBP_3'UTR",
  "CPEB4_CDS",
  "RPS11_CDS",
  "RPS15A_CDS",
  "VIM_CDS",
  "NIPBL_CDS",
  "FUBP3_CDS",
  "CENPC_CDS",
  "IFIT2_Intron",
  "PRPF8_Intron",
  "YBX3_Intron",
  "DDX6_Intron",
  "SRSF1_Intron",
  "ELAVL1_Intron",
  "DDX3X_3'UTR",
  "IGF2BP2_3'UTR",
  "CELF1_3'UTR",
  "TRA2A_5'UTR",
  "ZC3HAV1_5'UTR",
  "PABPC1_Intron",
  "EIF3D_3'UTR",
  "RPS19_5'UTR",
  "ZRANB2_5'UTR",
  "DICER1_CDS",
  "G3BP1_Intron",
  "HNRNPA1_Intron",
  "SRSF2_Intron",
  "OAS1_Intron",
  "FAM120A_3'UTR",
  "SRP14_3'UTR",
  "NCBP3_3'UTR",
  "RBM47_Intron",
  "DEK_Intron",
  "G3BP2_5'UTR",
  "HNRNPU_CDS",
  "SAMD4A_Intron",
  "MTDH_Intron",
  "RPS15A_3'UTR",
  "RPS19_3'UTR",
  "TRA2A_3'UTR",
  "SERBP1_5'UTR",
  "EIF4A1_5'UTR",
  "KHSRP_5'UTR",
  "RPRD2_5UTR",
  "FAM120A_5'UTR",
  "ZRANB2_Intron", 
  "XRN2_Intron",
  "ZRANB2_3'UTR",
  "TAF15_3'UTR",
  "IGF2BP2_CDS",
  "YBX3_CDS",
  "RTCB_CDS",
  "METAP2_3'UTR",
  "RPS3A_Intron",
  "YWHAZ_CDS",
  "NOP58_CDS",
  "CHD7_CDS",
  "TNRC6A_CDS",
  "FMR1_3'UTR"
)

# Optional: trim potential whitespace/invisible characters to avoid matching issues
df$Regulon <- trimws(df$Regulon)

#df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)
df$label <- ifelse(df$Regulon %in% keep_regs & df$pval < 0.05, df$Regulon, NA)

# ---- Plot ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # Reverse color scale: smaller padj = darker (use rev(pal))
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p






#########################################################################
#########################################################################
#4) Neutrophil count
#########################################################################
#########################################################################

# 1) OpenTargets CSV: specify which column is gene/score
#   Common: gene symbol column like "gene" / "geneSymbol"
#           score column like "overallScore" / "association_score" / "mlog10p"/"z"
ot_stats <- function_csv("02_neutro/01_opentargts_neutr_count.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # Change to your score field


#--3'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_3UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./02_neutro/fgsea_RBP_regulons_vs_OpenTargets_3UTR_neutr_count.csv")




#--5'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_5UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./02_neutro/fgsea_RBP_regulons_vs_OpenTargets_5UTR_neutr_count.csv")



#--CDS
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_CDS.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./02_neutro/fgsea_RBP_regulons_vs_OpenTargets_CDS_neutr_count.csv")




#--Introns
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_Introns.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./02_neutro/fgsea_RBP_regulons_vs_OpenTargets_Introns_neutr_count.csv")




#---Plot---GSEA---Dotplot---for neutrophil count
##_--------dotplot for regulon-opentargets enrichment


##_-----In this version, we merge results from all 4 mRNA regions, then display enrichment and label significant results
########------------monocyte count all L2G genes
df <- read.csv("./02_neutro/01_enrichment_GSEA_dotplot.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# Color palette: smaller padj = darker color
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# Labeling: only add labels for points with FDR < 0.05 (adjust threshold as needed)
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# Only label specified regulons
keep_regs <- c(
  "DDX21_5'UTR",
  "DDX3X_5'UTR",
  "QKI_5'UTR",
  "YTHDF3_5'UTR",
  "IGF2BP2_5'UTR",
  "AKAP1_5'UTR",
  "TNRC6C_5'UTR",
  "HNRNPC_5'UTR",
  "SRSF5_5'UTR",
  "DDX21_3'UTR",
  "YTHDF3_3'UTR",
  "RYBP_3'UTR",
  "CPEB4_CDS",
  "RPS11_CDS",
  "RPS15A_CDS",
  "VIM_CDS",
  "NIPBL_CDS",
  "FUBP3_CDS",
  "CENPC_CDS",
  "IFIT2_Intron",
  "PRPF8_Intron",
  "YBX3_Intron",
  "DDX6_Intron",
  "SRSF1_Intron",
  "ELAVL1_Intron",
  "DDX3X_3'UTR",
  "IGF2BP2_3'UTR",
  "CELF1_3'UTR",
  "TRA2A_5'UTR",
  "ZC3HAV1_5'UTR",
  "PABPC1_Intron",
  "EIF3D_3'UTR",
  "RPS19_5'UTR",
  "ZRANB2_5'UTR",
  "DICER1_CDS",
  "G3BP1_Intron",
  "HNRNPA1_Intron",
  "SRSF2_Intron",
  "OAS1_Intron",
  "FAM120A_3'UTR",
  "SRP14_3'UTR",
  "NCBP3_3'UTR",
  "RBM47_Intron",
  "DEK_Intron",
  "G3BP2_5'UTR",
  "HNRNPU_CDS",
  "SAMD4A_Intron",
  "MTDH_Intron",
  "RPS15A_3'UTR",
  "RPS19_3'UTR",
  "TRA2A_3'UTR",
  "SERBP1_5'UTR",
  "EIF4A1_5'UTR",
  "KHSRP_5'UTR",
  "RPRD2_5UTR",
  "FAM120A_5'UTR",
  "ZRANB2_Intron", 
  "XRN2_Intron",
  "ZRANB2_3'UTR",
  "TAF15_3'UTR",
  "IGF2BP2_CDS",
  "YBX3_CDS",
  "RTCB_CDS",
  "METAP2_3'UTR",
  "RPS3A_Intron",
  "YWHAZ_CDS",
  "NOP58_CDS",
  "CHD7_CDS",
  "TNRC6A_CDS",
  "FMR1_3'UTR"
)

# Optional: trim potential whitespace/invisible characters to avoid matching issues
df$Regulon <- trimws(df$Regulon)

#df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)
df$label <- ifelse(df$Regulon %in% keep_regs & df$pval < 0.05, df$Regulon, NA)

# ---- Plot ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # Reverse color scale: smaller padj = darker (use rev(pal))
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p








#########################################################################
#########################################################################
#5) WBC count
#########################################################################
#########################################################################

# 1) OpenTargets CSV: specify which column is gene/score
#   Common: gene symbol column like "gene" / "geneSymbol"
#           score column like "overallScore" / "association_score" / "mlog10p"/"z"
ot_stats <- function_csv("03_WBC/01_opentargts_WBC.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # Change to your score field


#--3'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_3UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./03_WBC/fgsea_RBP_regulons_vs_OpenTargets_3UTR_WBC.csv")




#--5'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_5UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./03_WBC/fgsea_RBP_regulons_vs_OpenTargets_5UTR_WBC.csv")



#--CDS
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_CDS.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./03_WBC/fgsea_RBP_regulons_vs_OpenTargets_CDS_WBC.csv")




#--Introns
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_Introns.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./03_WBC/fgsea_RBP_regulons_vs_OpenTargets_Introns_WBC.csv")






#---Plot---GSEA---Dotplot---for neutrophil count
##_--------dotplot for regulon-opentargets enrichment


##_-----In this version, we merge results from all 4 mRNA regions, then display enrichment and label significant results
########------------monocyte count all L2G genes
df <- read.csv("./03_WBC/01_enrichment_GSEA_dotplot.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# Color palette: smaller padj = darker color
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# Labeling: only add labels for points with FDR < 0.05 (adjust threshold as needed)
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# Only label specified regulons
keep_regs <- c(
  "DDX21_5'UTR",
  "DDX3X_5'UTR",
  "QKI_5'UTR",
  "YTHDF3_5'UTR",
  "IGF2BP2_5'UTR",
  "AKAP1_5'UTR",
  "TNRC6C_5'UTR",
  "HNRNPC_5'UTR",
  "SRSF5_5'UTR",
  "DDX21_3'UTR",
  "YTHDF3_3'UTR",
  "RYBP_3'UTR",
  "CPEB4_CDS",
  "RPS11_CDS",
  "RPS15A_CDS",
  "VIM_CDS",
  "NIPBL_CDS",
  "FUBP3_CDS",
  "CENPC_CDS",
  "IFIT2_Intron",
  "PRPF8_Intron",
  "YBX3_Intron",
  "DDX6_Intron",
  "SRSF1_Intron",
  "ELAVL1_Intron",
  "DDX3X_3'UTR",
  "IGF2BP2_3'UTR",
  "CELF1_3'UTR",
  "TRA2A_5'UTR",
  "ZC3HAV1_5'UTR",
  "PABPC1_Intron",
  "EIF3D_3'UTR",
  "RPS19_5'UTR",
  "ZRANB2_5'UTR",
  "DICER1_CDS",
  "G3BP1_Intron",
  "HNRNPA1_Intron",
  "SRSF2_Intron",
  "OAS1_Intron",
  "FAM120A_3'UTR",
  "SRP14_3'UTR",
  "NCBP3_3'UTR",
  "RBM47_Intron",
  "DEK_Intron",
  "G3BP2_5'UTR",
  "HNRNPU_CDS",
  "SAMD4A_Intron",
  "MTDH_Intron",
  "RPS15A_3'UTR",
  "RPS19_3'UTR",
  "TRA2A_3'UTR",
  "SERBP1_5'UTR",
  "EIF4A1_5'UTR",
  "KHSRP_5'UTR",
  "RPRD2_5UTR",
  "FAM120A_5'UTR",
  "ZRANB2_Intron", 
  "XRN2_Intron",
  "ZRANB2_3'UTR",
  "TAF15_3'UTR",
  "IGF2BP2_CDS",
  "YBX3_CDS",
  "RTCB_CDS",
  "METAP2_3'UTR",
  "RPS3A_Intron",
  "YWHAZ_CDS",
  "NOP58_CDS",
  "CHD7_CDS",
  "TNRC6A_CDS",
  "FMR1_3'UTR"
)

# Optional: trim potential whitespace/invisible characters to avoid matching issues
df$Regulon <- trimws(df$Regulon)

#df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)
df$label <- ifelse(df$Regulon %in% keep_regs & df$pval < 0.05, df$Regulon, NA)

# ---- Plot ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # Reverse color scale: smaller padj = darker (use rev(pal))
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p








#########################################################################
#########################################################################
#6) eosinophil count
#########################################################################
#########################################################################

# 1) OpenTargets CSV: specify which column is gene/score
#   Common: gene symbol column like "gene" / "geneSymbol"
#           score column like "overallScore" / "association_score" / "mlog10p"/"z"
ot_stats <- function_csv("04_eosinophil/01_opentargts_eosin_count.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # Change to your score field


#--3'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_3UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./04_eosinophil/fgsea_RBP_regulons_vs_OpenTargets_3UTR_eosin_count.csv")




#--5'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_5UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./04_eosinophil/fgsea_RBP_regulons_vs_OpenTargets_5UTR_eosin_count.csv")



#--CDS
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_CDS.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./04_eosinophil/fgsea_RBP_regulons_vs_OpenTargets_CDS_eosin_count.csv")




#--Introns
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_Introns.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./04_eosinophil/fgsea_RBP_regulons_vs_OpenTargets_Introns_eosin_count.csv")






#---Plot---GSEA---Dotplot---for eosinophil count
##_--------dotplot for regulon-opentargets enrichment


##_-----In this version, we merge results from all 4 mRNA regions, then display enrichment and label significant results
########------------monocyte count all L2G genes
df <- read.csv("./04_eosinophil/01_enrichment_GSEA_dotplot.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# Color palette: smaller padj = darker color
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# Labeling: only add labels for points with FDR < 0.05 (adjust threshold as needed)
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# Only label specified regulons
keep_regs <- c(
  "DDX21_5'UTR",
  "DDX3X_5'UTR",
  "QKI_5'UTR",
  "YTHDF3_5'UTR",
  "IGF2BP2_5'UTR",
  "AKAP1_5'UTR",
  "TNRC6C_5'UTR",
  "HNRNPC_5'UTR",
  "SRSF5_5'UTR",
  "DDX21_3'UTR",
  "YTHDF3_3'UTR",
  "RYBP_3'UTR",
  "CPEB4_CDS",
  "RPS11_CDS",
  "RPS15A_CDS",
  "VIM_CDS",
  "NIPBL_CDS",
  "FUBP3_CDS",
  "CENPC_CDS",
  "IFIT2_Intron",
  "PRPF8_Intron",
  "YBX3_Intron",
  "DDX6_Intron",
  "SRSF1_Intron",
  "ELAVL1_Intron",
  "DDX3X_3'UTR",
  "IGF2BP2_3'UTR",
  "CELF1_3'UTR",
  "TRA2A_5'UTR",
  "ZC3HAV1_5'UTR",
  "PABPC1_Intron",
  "EIF3D_3'UTR",
  "RPS19_5'UTR",
  "ZRANB2_5'UTR",
  "DICER1_CDS",
  "G3BP1_Intron",
  "HNRNPA1_Intron",
  "SRSF2_Intron",
  "OAS1_Intron",
  "FAM120A_3'UTR",
  "SRP14_3'UTR",
  "NCBP3_3'UTR",
  "RBM47_Intron",
  "DEK_Intron",
  "G3BP2_5'UTR",
  "HNRNPU_CDS",
  "SAMD4A_Intron",
  "MTDH_Intron",
  "RPS15A_3'UTR",
  "RPS19_3'UTR",
  "TRA2A_3'UTR",
  "SERBP1_5'UTR",
  "EIF4A1_5'UTR",
  "KHSRP_5'UTR",
  "RPRD2_5UTR",
  "FAM120A_5'UTR",
  "ZRANB2_Intron", 
  "XRN2_Intron",
  "ZRANB2_3'UTR",
  "TAF15_3'UTR",
  "IGF2BP2_CDS",
  "YBX3_CDS",
  "RTCB_CDS",
  "METAP2_3'UTR",
  "RPS3A_Intron",
  "YWHAZ_CDS",
  "NOP58_CDS",
  "CHD7_CDS",
  "TNRC6A_CDS",
  "FMR1_3'UTR"
)

# Optional: trim potential whitespace/invisible characters to avoid matching issues
df$Regulon <- trimws(df$Regulon)

#df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)
df$label <- ifelse(df$Regulon %in% keep_regs & df$pval < 0.05, df$Regulon, NA)

# ---- Plot ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # Reverse color scale: smaller padj = darker (use rev(pal))
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p






#########################################################################
#########################################################################
#7) basophil count
#########################################################################
#########################################################################

# 1) OpenTargets CSV: specify which column is gene/score
#   Common: gene symbol column like "gene" / "geneSymbol"
#           score column like "overallScore" / "association_score" / "mlog10p"/"z"
ot_stats <- function_csv("./05_basophil/01_opentargts_baso_count.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # Change to your score field


#--3'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_3UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./05_basophil/fgsea_RBP_regulons_vs_OpenTargets_3UTR_baso_count.csv")




#--5'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_5UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./05_basophil/fgsea_RBP_regulons_vs_OpenTargets_5UTR_baso_count.csv")



#--CDS
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_CDS.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./05_basophil/fgsea_RBP_regulons_vs_OpenTargets_CDS_baso_count.csv")




#--Introns
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_Introns.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./05_basophil/fgsea_RBP_regulons_vs_OpenTargets_Introns_baso_count.csv")





#---Plot---GSEA---Dotplot---for basophil count
##_--------dotplot for regulon-opentargets enrichment


##_-----In this version, we merge results from all 4 mRNA regions, then display enrichment and label significant results
########------------monocyte count all L2G genes
df <- read.csv("./05_basophil/01_enrichment_GSEA_dotplot.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# Color palette: smaller padj = darker color
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# Labeling: only add labels for points with FDR < 0.05 (adjust threshold as needed)
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# Only label specified regulons
keep_regs <- c(
  "DDX21_5'UTR",
  "DDX3X_5'UTR",
  "QKI_5'UTR",
  "YTHDF3_5'UTR",
  "IGF2BP2_5'UTR",
  "AKAP1_5'UTR",
  "TNRC6C_5'UTR",
  "HNRNPC_5'UTR",
  "SRSF5_5'UTR",
  "DDX21_3'UTR",
  "YTHDF3_3'UTR",
  "RYBP_3'UTR",
  "CPEB4_CDS",
  "RPS11_CDS",
  "RPS15A_CDS",
  "VIM_CDS",
  "NIPBL_CDS",
  "FUBP3_CDS",
  "CENPC_CDS",
  "IFIT2_Intron",
  "PRPF8_Intron",
  "YBX3_Intron",
  "DDX6_Intron",
  "SRSF1_Intron",
  "ELAVL1_Intron",
  "DDX3X_3'UTR",
  "IGF2BP2_3'UTR",
  "CELF1_3'UTR",
  "TRA2A_5'UTR",
  "ZC3HAV1_5'UTR",
  "PABPC1_Intron",
  "EIF3D_3'UTR",
  "RPS19_5'UTR",
  "ZRANB2_5'UTR",
  "DICER1_CDS",
  "G3BP1_Intron",
  "HNRNPA1_Intron",
  "SRSF2_Intron",
  "OAS1_Intron",
  "FAM120A_3'UTR",
  "SRP14_3'UTR",
  "NCBP3_3'UTR",
  "RBM47_Intron",
  "DEK_Intron",
  "G3BP2_5'UTR",
  "HNRNPU_CDS",
  "SAMD4A_Intron",
  "MTDH_Intron",
  "RPS15A_3'UTR",
  "RPS19_3'UTR",
  "TRA2A_3'UTR",
  "SERBP1_5'UTR",
  "EIF4A1_5'UTR",
  "KHSRP_5'UTR",
  "RPRD2_5UTR",
  "FAM120A_5'UTR",
  "ZRANB2_Intron", 
  "XRN2_Intron",
  "ZRANB2_3'UTR",
  "TAF15_3'UTR",
  "IGF2BP2_CDS",
  "YBX3_CDS",
  "RTCB_CDS",
  "METAP2_3'UTR",
  "RPS3A_Intron",
  "YWHAZ_CDS",
  "NOP58_CDS",
  "CHD7_CDS",
  "TNRC6A_CDS",
  "FMR1_3'UTR"
)

# Optional: trim potential whitespace/invisible characters to avoid matching issues
df$Regulon <- trimws(df$Regulon)

#df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)
df$label <- ifelse(df$Regulon %in% keep_regs & df$pval < 0.05, df$Regulon, NA)

# ---- Plot ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # Reverse color scale: smaller padj = darker (use rev(pal))
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p



#########################################################################
#########################################################################
#8) MCHC  
#########################################################################
#########################################################################

# 1) OpenTargets CSV: specify which column is gene/score
#   Common: gene symbol column like "gene" / "geneSymbol"
#           score column like "overallScore" / "association_score" / "mlog10p"/"z"
ot_stats <- function_csv("./06_MCHC/01_opentargts_MCHC.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # Change to your score field


#--3'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_3UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./06_MCHC/fgsea_RBP_regulons_vs_OpenTargets_3UTR_MCHC.csv")




#--5'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_5UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./06_MCHC/fgsea_RBP_regulons_vs_OpenTargets_5UTR_MCHC.csv")




#--CDS
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_CDS.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./06_MCHC/fgsea_RBP_regulons_vs_OpenTargets_CDS_MCHC.csv")





#--Introns
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_Introns.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./06_MCHC/fgsea_RBP_regulons_vs_OpenTargets_Introns_MCHC.csv")






#---Plot---GSEA---Dotplot---for basophil count
##_--------dotplot for regulon-opentargets enrichment


##_-----In this version, we merge results from all 4 mRNA regions, then display enrichment and label significant results
########------------monocyte count all L2G genes
df <- read.csv("./06_MCHC/01_enrichment_GSEA_dotplot.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# Color palette: smaller padj = darker color
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# Labeling: only add labels for points with FDR < 0.05 (adjust threshold as needed)
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# Only label specified regulons
keep_regs <- c(
  "DDX21_5'UTR",
  "DDX3X_5'UTR",
  "QKI_5'UTR",
  "YTHDF3_5'UTR",
  "IGF2BP2_5'UTR",
  "AKAP1_5'UTR",
  "TNRC6C_5'UTR",
  "HNRNPC_5'UTR",
  "SRSF5_5'UTR",
  "DDX21_3'UTR",
  "YTHDF3_3'UTR",
  "RYBP_3'UTR",
  "CPEB4_CDS",
  "RPS11_CDS",
  "RPS15A_CDS",
  "VIM_CDS",
  "NIPBL_CDS",
  "FUBP3_CDS",
  "CENPC_CDS",
  "IFIT2_Intron",
  "PRPF8_Intron",
  "YBX3_Intron",
  "DDX6_Intron",
  "SRSF1_Intron",
  "ELAVL1_Intron",
  "DDX3X_3'UTR",
  "IGF2BP2_3'UTR",
  "CELF1_3'UTR",
  "TRA2A_5'UTR",
  "ZC3HAV1_5'UTR",
  "PABPC1_Intron",
  "EIF3D_3'UTR",
  "RPS19_5'UTR",
  "ZRANB2_5'UTR",
  "DICER1_CDS",
  "G3BP1_Intron",
  "HNRNPA1_Intron",
  "SRSF2_Intron",
  "OAS1_Intron",
  "FAM120A_3'UTR",
  "SRP14_3'UTR",
  "NCBP3_3'UTR",
  "RBM47_Intron",
  "DEK_Intron",
  "G3BP2_5'UTR",
  "HNRNPU_CDS",
  "SAMD4A_Intron",
  "MTDH_Intron",
  "RPS15A_3'UTR",
  "RPS19_3'UTR",
  "TRA2A_3'UTR",
  "SERBP1_5'UTR",
  "EIF4A1_5'UTR",
  "KHSRP_5'UTR",
  "RPRD2_5UTR",
  "FAM120A_5'UTR",
  "ZRANB2_Intron", 
  "XRN2_Intron",
  "ZRANB2_3'UTR",
  "TAF15_3'UTR",
  "IGF2BP2_CDS",
  "YBX3_CDS",
  "RTCB_CDS",
  "METAP2_3'UTR",
  "RPS3A_Intron",
  "YWHAZ_CDS",
  "NOP58_CDS",
  "CHD7_CDS",
  "TNRC6A_CDS",
  "FMR1_3'UTR"
)

# Optional: trim potential whitespace/invisible characters to avoid matching issues
df$Regulon <- trimws(df$Regulon)

#df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)
df$label <- ifelse(df$Regulon %in% keep_regs & df$pval < 0.05, df$Regulon, NA)

# ---- Plot ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # Reverse color scale: smaller padj = darker (use rev(pal))
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p






#########################################################################
#########################################################################
#9) HL
#########################################################################
#########################################################################

# 1) OpenTargets CSV: specify which column is gene/score
#   Common: gene symbol column like "gene" / "geneSymbol"
#           score column like "overallScore" / "association_score" / "mlog10p"/"z"
ot_stats <- function_csv("./07_HL/01_opentargts_HL.csv",
                         gene_col  = "symbol",
                         score_col = "globalScore")         # Change to your score field


#--3'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_3UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./07_HL/fgsea_RBP_regulons_vs_OpenTargets_3UTR_HL.csv")




#--5'UTR
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_5UTR.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)


# 5) Export results
data.table::fwrite(fg_res, "./07_HL/fgsea_RBP_regulons_vs_OpenTargets_5UTR_HL.csv")




#--CDS
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_CDS.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./07_HL/fgsea_RBP_regulons_vs_OpenTargets_CDS_HL.csv")





#--Introns
# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_Introns.gmt")


# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 5, maxSize = 5000, nperm = 100000)

fg_res$padj_BY <- p.adjust(fg_res$pval, method = "BY")
fg_res$padj_BH <- p.adjust(fg_res$pval, method = "BH")

head(fg_res)
# 5) Export results
data.table::fwrite(fg_res, "./07_HL/fgsea_RBP_regulons_vs_OpenTargets_Introns_HL.csv")






#---Plot---GSEA---Dotplot---for basophil count
##_--------dotplot for regulon-opentargets enrichment


##_-----In this version, we merge results from all 4 mRNA regions, then display enrichment and label significant results
########------------monocyte count all L2G genes
df <- read.csv("./07_HL/01_enrichment_GSEA_dotplot.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# Color palette: smaller padj = darker color
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# Labeling: only add labels for points with FDR < 0.05 (adjust threshold as needed)
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# Only label specified regulons
keep_regs <- c(
  "DDX21_5'UTR",
  "DDX3X_5'UTR",
  "QKI_5'UTR",
  "YTHDF3_5'UTR",
  "IGF2BP2_5'UTR",
  "AKAP1_5'UTR",
  "TNRC6C_5'UTR",
  "HNRNPC_5'UTR",
  "SRSF5_5'UTR",
  "DDX21_3'UTR",
  "YTHDF3_3'UTR",
  "RYBP_3'UTR",
  "CPEB4_CDS",
  "RPS11_CDS",
  "RPS15A_CDS",
  "VIM_CDS",
  "NIPBL_CDS",
  "FUBP3_CDS",
  "CENPC_CDS",
  "IFIT2_Intron",
  "PRPF8_Intron",
  "YBX3_Intron",
  "DDX6_Intron",
  "SRSF1_Intron",
  "ELAVL1_Intron",
  "DDX3X_3'UTR",
  "IGF2BP2_3'UTR",
  "CELF1_3'UTR",
  "TRA2A_5'UTR",
  "ZC3HAV1_5'UTR",
  "PABPC1_Intron",
  "EIF3D_3'UTR",
  "RPS19_5'UTR",
  "ZRANB2_5'UTR",
  "DICER1_CDS",
  "G3BP1_Intron",
  "HNRNPA1_Intron",
  "SRSF2_Intron",
  "OAS1_Intron",
  "FAM120A_3'UTR",
  "SRP14_3'UTR",
  "NCBP3_3'UTR",
  "RBM47_Intron",
  "DEK_Intron",
  "G3BP2_5'UTR",
  "HNRNPU_CDS",
  "SAMD4A_Intron",
  "MTDH_Intron",
  "RPS15A_3'UTR",
  "RPS19_3'UTR",
  "TRA2A_3'UTR",
  "SERBP1_5'UTR",
  "EIF4A1_5'UTR",
  "KHSRP_5'UTR",
  "RPRD2_5UTR",
  "FAM120A_5'UTR",
  "ZRANB2_Intron", 
  "XRN2_Intron",
  "ZRANB2_3'UTR",
  "TAF15_3'UTR",
  "IGF2BP2_CDS",
  "YBX3_CDS",
  "RTCB_CDS",
  "METAP2_3'UTR",
  "RPS3A_Intron",
  "YWHAZ_CDS",
  "NOP58_CDS",
  "CHD7_CDS",
  "TNRC6A_CDS",
  "FMR1_3'UTR"
)

# Optional: trim potential whitespace/invisible characters to avoid matching issues
df$Regulon <- trimws(df$Regulon)

#df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)
df$label <- ifelse(df$Regulon %in% keep_regs & df$pval < 0.05, df$Regulon, NA)

# ---- Plot ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # Reverse color scale: smaller padj = darker (use rev(pal))
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
