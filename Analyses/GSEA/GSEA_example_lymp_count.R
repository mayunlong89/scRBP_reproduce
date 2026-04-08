#2025-09-08

#lymphocyte count

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
                     minSize = 10, maxSize = 5000, nperm = 1000) {
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

# =========================
# Usage example (replace with your file paths and column names)
# =========================
# 1) OpenTargets CSV: specify which column is gene/score
#   Common: gene symbol column like "gene" / "geneSymbol"
#           score column like "overallScore" / "association_score" / "mlog10p"/"z"
ot_stats <- function_csv("01_opentargts_lymp_count.csv",
                         gene_col  = "symbol",
                         score_col = "gwasCredibleSets")         # Change to your score field

# 2) Regulons GMT: contains multiple RBP_regulon entries
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_3UTR.gmt")
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_5UTR.gmt")
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_CDS.gmt")
regulons <- read_gmt("../00-6.5K_100runs_regulons-GO-KEGG/scRBP_prune_matched_results_min10genes.symbol_default_6.5K_Introns.gmt")

# 3) Run GSEA (recommend nperm at least 10000; use 100000 for final results)
fg_res <- run_gsea(ot_stats, regulons, minSize = 10, maxSize = 5000, nperm = 100000)
head(fg_res)

  
# 5) Export results
data.table::fwrite(fg_res, "fgsea_RBP_regulons_vs_OpenTargets_3UTR_lymphocyte_count.csv")
data.table::fwrite(fg_res, "fgsea_RBP_regulons_vs_OpenTargets_5UTR_lymphocyte_count.csv")
data.table::fwrite(fg_res, "fgsea_RBP_regulons_vs_OpenTargets_CDS_lymphocyte_count.csv")
data.table::fwrite(fg_res, "fgsea_RBP_regulons_vs_OpenTargets_Introns_lymphocyte_count.csv")


# Find the CELF1 regulon row (name may be "CELF1" or "CELF1_regulon")
row <- fg_res[fg_res$pathway %in% c("CELF1", "CELF1_regulon"), ][1, ]

# Extract leading edge gene vector
le <- unlist(row$leadingEdge)
length(le); le[1:10]  


# Assuming fg_res is the result from run_gsea(...); ot_stats is a named vector (sorted descending)
row <- fg_res[fg_res$pathway %in% c("QKI","QKI_regulon"), ][1, ]
row[, c("pathway","NES","size","pval","padj")]

# Extract leading edge and list by OpenTargets rank (front to back)
le <- unlist(row$leadingEdge)
ord <- order(rank(ot_stats))                               # order used by fgsea
le_tab <- data.frame(
  gene = le,
  score = ot_stats[le],
  rank  = match(le, names(ot_stats))
)[order(le_tab$rank), ]
nrow(le_tab)                         # number of leading edge genes
min(le_tab$rank) / length(ot_stats)  # coverage start (proportion of total ranking)
max(le_tab$rank) / length(ot_stats)  # coverage end
head(le_tab, 15)                     # top 15 leading edge genes


# 4) Plot enrichment curve for a specific regulon (e.g., YBX3_regulon)
# 4) Plot enrichment curve for a specific regulon (e.g., YBX3_regulon)
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

# Utility to coerce data.frame / vector into a named numeric vector
coerce_stats_vec <- function(stats, gene_col = NULL, score_col = NULL, to_upper = TRUE) {
  if (is.numeric(stats) && !is.null(names(stats))) {
    v <- stats
  } else if (is.data.frame(stats)) {
    stopifnot(!is.null(gene_col), !is.null(score_col))
    df <- stats[, c(gene_col, score_col)]
    colnames(df) <- c("gene","score")
    df <- df[!is.na(df$gene) & !is.na(df$score), ]
    if (to_upper) df$gene <- toupper(trimws(df$gene))
    # Take max score for duplicate genes
    df <- aggregate(score ~ gene, df, max)
    v <- df$score; names(v) <- df$gene
  } else {
    stop("`stats` must be a named numeric vector or a data.frame/tibble.")
  }
  v <- v[is.finite(v)]
  if (to_upper) names(v) <- toupper(names(v))
  # Deduplicate: keep max for same name
  v <- tapply(v, names(v), max)
  v[order(v, decreasing = TRUE)]
}

# Plot enrichment with correctly annotated LE (red) and non-LE (grey)
plot_enrichment_with_le <- function(pathway_name, pathways, stats,
                                    fg_res = NULL,
                                    gene_col = NULL, score_col = NULL,
                                    to_upper = TRUE,
                                    le_col = "red", other_col = "grey40") {
  
  # 1) Ranked list (named vector, descending)
  geneList <- coerce_stats_vec(stats, gene_col, score_col, to_upper)
  
  # 2) Get gene set (compatible with names like "QKI" or "QKI_regulon")
  set_genes <- pathways[[pathway_name]]
  if (is.null(set_genes)) {
    alt <- sub("_regulon$", "", pathway_name)
    if (!identical(alt, pathway_name)) set_genes <- pathways[[alt]]
  }
  if (is.null(set_genes)) stop("Cannot find gene set in pathways: ", pathway_name)
  if (to_upper) set_genes <- toupper(set_genes)
  set_genes <- unique(set_genes)
  
  # 3) Overlap check
  ov <- intersect(names(geneList), set_genes)
  message(sprintf("Set size = %d, universe = %d, overlap = %d",
                  length(set_genes), length(geneList), length(ov)))
  if (length(ov) == 0) stop("No overlap between gene set and ranked list (check case or ID mismatch).")
  
  # 4) Leading edge: prefer from fg_res; if unavailable, run fgsea temporarily
  le_genes <- NULL
  if (!is.null(fg_res) && "pathway" %in% names(fg_res)) {
    row <- fg_res[fg_res$pathway %in% c(pathway_name, sub("_regulon$", "", pathway_name)), ]
    if (nrow(row) > 0) le_genes <- unlist(row$leadingEdge)
  }
  if (is.null(le_genes)) {
    tmp <- fgsea(pathways = list(tmp = set_genes), stats = geneList, nperm = 10000)
    le_genes <- unlist(tmp$leadingEdge)
  }
  
  # 5) Positions (consistent with plotEnrichment)
  all_idx <- which(names(geneList) %in% set_genes)
  le_idx  <- which(names(geneList) %in% le_genes)
  
  # 6) Plot
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
plot_enrichment_with_le("PUM1", regulons, ot_stats, fg_res)
plot_enrichment_with_le("HNRNPC", regulons, ot_stats, fg_res)
plot_enrichment_with_le("QKI", regulons, ot_stats, fg_res)

#---Introns
plot_enrichment_with_le("RBM47", regulons, ot_stats, fg_res)








# List their OpenTargets scores and ranks together
library(dplyr)
le_tab <- tibble(
  gene = le,
  score = ot_stats[le],
  rank  = match(le, names(ot_stats))  # ot_stats is already sorted descending
) |> arrange(rank)
le_tab


# stats is sorted descending; pathways is your regulon list
ord <- order(rank(ot_stats))                       # order used by fgsea
pos_all <- which(names(ot_stats)[ord] %in% regulons[["CELF1"]])
length(pos_all)    # should equal 28
pos_all[1:20]      # inspect the first 20 positions


p <- fgsea::plotEnrichment(regulons[["CELF1"]], ot_stats, ticksSize = 0.15) +
  ggtitle("CELF1_regulon ~ OpenTargets ranks") +
  coord_cartesian(clip = "off")                # prevent clipping

# Add an extra rug layer at the bottom (one tick per gene; less likely to be obscured)
ticks_df <- data.frame(x = pos_all)
p <- p + ggplot2::geom_rug(data = ticks_df, aes(x = x), sides = "b",
                           inherit.aes = FALSE, alpha = 0.7)
p


# Get leading edge
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


# Export a wider figure to reduce overlap
# ggsave("CELF1_gsea.png", p, width = 12, height = 5, dpi = 300)




##_--------dotplot for regulon-opentargets enrichment
## -------- dotplot for regulon-OpenTargets enrichment (new) --------
txt <- "Regulon\tNES\tlogp\tpadj
DDX3X_3UTR\t1.306921251\t2.616188977\t0.079501292
TAF15_3UTR\t1.174488275\t1.039295791\t0.284197158
DDX3X_5UTR\t1.346593683\t5.000004343\t0.000359996
QKI_5UTR\t1.30623445\t5.000004343\t0.000359996
PUM1_5UTR\t1.402208323\t4.522883088\t0.000539995
HNRNPC_5UTR\t1.246188233\t3.769555422\t0.00203998
TNRC6C_5UTR\t1.280640031\t2.931818481\t0.012034165
FMR1_5UTR\t1.306860281\t1.292317846\t0.244863127
NIPBL_5UTR\t1.143172285\t1.237100814\t0.245348135
TAF15_5UTR\t1.212696977\t0.816234503\t0.37905305
EEF2_5UTR\t0.719425736\t0.048221226\t0.986888245
NIPBL_CDS\t1.237409638\t3.050614336\t0.037222283
RBM47_Introns\t1.268714439\t2.790489328\t0.089099109
SRSF1_Introns\t1.015484889\t0.344242474\t0.842251337
LARP7_Introns\t0.969313125\t0.265447445\t0.863661961
HNRNPC_Introns\t0.953586735\t0.133250462\t0.94704208
DEK_Introns\t0.833076888\t0.032418929\t0.981613261
"

df <- read.table(text = txt, header = TRUE, sep = "\t", check.names = FALSE,
                 stringsAsFactors = FALSE)

# Color scale: larger NES = darker color
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)

# Only add labels for significant points (adjust threshold as needed)
df$label <- ifelse(df$logp > -log10(0.05), df$Regulon, NA)

library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color =  NES)) +
  geom_point(size = 4) +
  scale_color_gradientn(colors = pal, trans = "log10",
                        breaks = c(0.6, 0.8, 1.0, 1.2, 1.4),
                        #labels = scales::label_scientific(digits = 2),
                        name = " NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 3.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 12)

p
# ggsave("regulon_scatter_new.png", p, width = 6, height = 6, dpi = 300)
