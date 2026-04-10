suppressPackageStartupMessages({
  library(fgsea)
  library(ggplot2)
})

# ---------- 1) robust: coerce stats to named numeric vector ----------
coerce_stats_vec <- function(stats, gene_col = NULL, score_col = NULL, to_upper = TRUE) {
  
  if (is.data.frame(stats)) {
    stopifnot(!is.null(gene_col), !is.null(score_col))
    df <- stats[, c(gene_col, score_col)]
    colnames(df) <- c("gene", "score")
    df <- df[!is.na(df$gene) & !is.na(df$score), ]
    df$gene <- as.character(df$gene)
    if (to_upper) df$gene <- toupper(trimws(df$gene))
    df$gene <- trimws(df$gene)
    
    # remove empty gene
    df <- df[df$gene != "" & !is.na(df$gene), ]
    
    # collapse duplicates by max score
    df <- aggregate(score ~ gene, df, max)
    
    v <- df$score
    names(v) <- df$gene
    
  } else {
    # try to coerce vector / matrix / array -> numeric vector
    v <- as.numeric(stats)
    nm <- names(stats)
    
    if (is.null(nm)) {
      # if stats is matrix with rownames
      if (!is.null(rownames(stats)) && length(rownames(stats)) == length(v)) {
        nm <- rownames(stats)
      } else {
        stop("`stats` 必须是带 names 的数值向量，或 data.frame(指定 gene_col/score_col)。")
      }
    }
    
    nm <- as.character(nm)
    if (length(nm) != length(v)) {
      stop(sprintf("stats 长度(%d)与names长度(%d)不一致，请检查 ot_stats 是否被转成了 matrix/df。",
                   length(v), length(nm)))
    }
    
    if (to_upper) nm <- toupper(trimws(nm))
    nm <- trimws(nm)
    
    keep <- is.finite(v) & !is.na(nm) & nm != ""
    v <- v[keep]
    nm <- nm[keep]
    names(v) <- nm
    
    # collapse duplicates by max (robust, no tapply-length pitfalls)
    sp <- split(v, names(v))
    v <- vapply(sp, function(x) max(x, na.rm = TRUE), numeric(1))
  }
  
  v <- sort(v, decreasing = TRUE)
  return(v)
}

# ---------- 2) main plotting function ----------
plot_enrichment_with_le <- function(pathway_name, pathways, stats,
                                    fg_res = NULL,
                                    gene_col = NULL, score_col = NULL,
                                    to_upper = TRUE,
                                    nperm_if_need = 20000,
                                    le_col = "red3",
                                    other_col = "grey50",
                                    tick_other_h = 0.035,
                                    tick_le_h = 0.07,
                                    add_label = TRUE) {
  
  # (A) gene ranks
  geneList <- coerce_stats_vec(stats, gene_col, score_col, to_upper)
  
  # (B) get pathway genes (allow "DDX3X" or "DDX3X_regulon")
  set_genes <- pathways[[pathway_name]]
  alt1 <- sub("_regulon$", "", pathway_name)
  alt2 <- paste0(pathway_name, "_regulon")
  if (is.null(set_genes) && !identical(alt1, pathway_name)) set_genes <- pathways[[alt1]]
  if (is.null(set_genes)) set_genes <- pathways[[alt2]]
  if (is.null(set_genes)) stop("在 pathways 里找不到集合：", pathway_name)
  
  set_genes <- unique(as.character(set_genes))
  if (to_upper) set_genes <- toupper(trimws(set_genes))
  
  # overlap check
  ov <- intersect(names(geneList), set_genes)
  message(sprintf("[overlap] set=%d, universe=%d, overlap=%d",
                  length(set_genes), length(geneList), length(ov)))
  if (length(ov) == 0) stop("集合与排行无交集：请检查 gene ID/大小写。")
  
  # (C) leading edge + NES/pval (prefer fg_res, otherwise run fgsea once)
  le_genes <- NULL
  NES <- NA_real_
  pval <- NA_real_
  
  if (!is.null(fg_res) && "pathway" %in% colnames(fg_res)) {
    hit <- fg_res[fg_res$pathway %in% c(pathway_name, alt1, alt2), ]
    if (nrow(hit) > 0) {
      hit <- hit[order(hit$pval), ][1, ]
      le_genes <- unlist(hit$leadingEdge)
      NES <- hit$NES
      pval <- hit$pval
    }
  }
  
  if (is.null(le_genes)) {
    tmp <- fgsea(pathways = list(tmp = set_genes),
                 stats = geneList,
                 nperm = nperm_if_need)
    le_genes <- unlist(tmp$leadingEdge)
    NES <- tmp$NES[1]
    pval <- tmp$pval[1]
  }
  
  if (to_upper) le_genes <- toupper(le_genes)
  
  # (D) positions in ranked list
  all_idx <- which(names(geneList) %in% set_genes)
  le_idx  <- which(names(geneList) %in% le_genes)
  other_idx <- setdiff(all_idx, le_idx)
  
  # (E) base enrichment plot
  p <- fgsea::plotEnrichment(set_genes, geneList, ticksSize = 0) +
    ggtitle(paste0(pathway_name, " ~ OpenTargets ranks")) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(5.5, 5.5, 18, 5.5)  # bottom extra for ticks
    ) +
    coord_cartesian(clip = "off")
  
  # (F) add custom ticks (grey=non-LE, red=LE)
  p <- p +
    geom_segment(
      data = data.frame(x = other_idx),
      aes(x = x, xend = x, y = 0, yend = -tick_other_h),
      inherit.aes = FALSE,
      linewidth = 0.25,
      color = other_col
    ) +
    geom_segment(
      data = data.frame(x = le_idx),
      aes(x = x, xend = x, y = 0, yend = -tick_le_h),
      inherit.aes = FALSE,
      linewidth = 0.4,
      color = le_col
    )
  
  # (G) optional label
  if (add_label && is.finite(pval) && is.finite(NES)) {
    lab <- sprintf("P = %.2g, NES = %.2f", pval, NES)
    p <- p + annotate("text",
                      x = floor(length(geneList) * 0.08),
                      y = max(p$data$enrichmentScore, na.rm = TRUE) * 0.85,
                      label = lab,
                      hjust = 0,
                      size = 4)
  }
  
  return(p)
}
