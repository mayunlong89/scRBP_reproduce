# --- 2026-01-14
# GO enrichment analysis for developmental regulon trend clusters
# Two gene set strategies per cluster:
#   (1) Union: all target genes across regulons in the cluster
#   (2) Core: genes appearing in >= threshold fraction of regulons (frequency-filtered)

# ============================================================
# Load required packages
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(forcats)
})

# ============================================================
# 0) Input / Output / Parameters
# ============================================================
clusters_csv <- "RegulonTrends_global_v3.clusters.csv"
gmt_file     <- "../03_fetal_brain_regulons_4_mRNA_regions/summary_all_4regions_regulon_min10genes.gmt"
out_dir      <- "GO_by_trend_clusters3_freqCore"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Number of representative GO terms per cluster
top_k <- 5

# GO enrichment thresholds
p_cut           <- 0.05
q_cut           <- 0.20
min_gs_go       <- 10     # Minimum gene set size for enrichGO (recommend >= 20 for stability)
simplify_cutoff <- 0.7

# Which gene set strategies to run
run_union_GO <- TRUE
run_core_GO  <- TRUE

# Core gene selection rule
# "prop": genes in >= core_min_prop fraction of regulons; "k": genes in >= core_min_k regulons
core_rule_mode <- "prop"
core_min_prop  <- 0.05    # Recommended: 0.15-0.30; use 0.10-0.20 for large clusters
core_min_k     <- 3       # Recommended: >= 3 or >= 5
core_max_genes <- Inf     # Optional cap on core set size (set Inf to disable)

# ============================================================
# 1) Read cluster assignments
# ============================================================
clu <- fread(clusters_csv, data.table = FALSE)
stopifnot(all(c("regulon", "cluster") %in% colnames(clu)))
clu$cluster <- as.integer(clu$cluster)

# ============================================================
# 2) Read GMT file: regulon -> target gene list
#    GMT format: setName <TAB> description <TAB> gene1 <TAB> gene2 ...
# ============================================================
read_gmt_local <- function(gmt_path) {
  lines <- readLines(gmt_path)
  lst <- lapply(lines, function(x) {
    parts <- strsplit(x, "\t", fixed = TRUE)[[1]]
    set_name <- parts[1]
    genes <- unique(parts[-c(1, 2)])
    genes <- genes[genes != "" & !is.na(genes)]
    list(term = set_name, genes = genes)
  })
  names(lst) <- vapply(lst, `[[`, character(1), "term")
  lapply(lst, `[[`, "genes")
}
gmt <- read_gmt_local(gmt_file)

# Check overlap between cluster regulons and GMT entries
common_regs <- intersect(clu$regulon, names(gmt))
if (length(common_regs) < 10) {
  stop("Too few regulons overlap between clusters_csv and gmt. Please check naming consistency.")
}
message("[INFO] regulons in clusters: ", nrow(clu))
message("[INFO] regulons in gmt:      ", length(gmt))
message("[INFO] overlap regulons:     ", length(common_regs))

# Keep only regulons present in GMT
clu2 <- clu %>% dplyr::filter(regulon %in% names(gmt))

# ============================================================
# 3) Build cluster-level union gene sets + regulon size statistics
# ============================================================
reg_size_df <- data.frame(
  regulon = names(gmt),
  n_genes_regulon = vapply(gmt, length, integer(1)),
  stringsAsFactors = FALSE
)

cluster_union_tbl <- clu2 %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarize(
    regulons        = list(regulon),
    genes_union     = list(unique(unlist(gmt[regulon]))),
    n_regulons      = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    n_targets_union     = length(genes_union),
    regulon_sizes       = list(reg_size_df$n_genes_regulon[match(regulons, reg_size_df$regulon)]),
    mean_regulon_size   = mean(unlist(regulon_sizes), na.rm = TRUE),
    median_regulon_size = stats::median(unlist(regulon_sizes), na.rm = TRUE),
    union_over_mean     = n_targets_union / mean_regulon_size,
    union_over_median   = n_targets_union / median_regulon_size
  ) %>%
  dplyr::ungroup()

# Print union size summary (sanity check for overly large unions)
cluster_summary <- cluster_union_tbl %>%
  dplyr::select(cluster, n_regulons, n_targets_union,
                mean_regulon_size, median_regulon_size,
                union_over_mean, union_over_median) %>%
  dplyr::arrange(dplyr::desc(n_targets_union))

message("\n[SUMMARY] Cluster UNION target-gene sizes:")
print(as.data.frame(cluster_summary), row.names = FALSE)

fwrite(cluster_summary, file = file.path(out_dir, "trendCluster_geneCounts.summary.csv"))

# Save union gene sets
cluster_union_out <- cluster_union_tbl %>%
  dplyr::mutate(
    regulons    = sapply(regulons, function(x) paste(x, collapse = ";")),
    genes_union = sapply(genes_union, function(x) paste(x, collapse = ";"))
  )
fwrite(cluster_union_out, file = file.path(out_dir, "trendCluster_geneSets.union.csv"))

# ============================================================
# 4) Compute within-cluster gene frequency and build CORE gene sets
#    Frequency = number of regulons (within cluster) containing each gene
#    Proportion = frequency / n_regulons
# ============================================================
calc_gene_freq_one_cluster <- function(cluster_id, regulons_vec, gmt_list) {
  dt_list <- lapply(regulons_vec, function(rg) {
    genes <- gmt_list[[rg]]
    if (is.null(genes) || length(genes) == 0) return(NULL)
    data.table(regulon = rg, gene = genes)
  })
  dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  if (nrow(dt) == 0) return(NULL)

  n_reg <- length(regulons_vec)
  freq_dt <- dt[, .(freq = uniqueN(regulon)), by = gene][, prop := freq / n_reg][]
  freq_dt[, cluster := cluster_id]
  freq_dt[, n_regulons := n_reg]
  setcolorder(freq_dt, c("cluster", "gene", "freq", "prop", "n_regulons"))
  freq_dt[order(-prop, -freq, gene)]
}

freq_all  <- list()
core_list <- list()

for (i in seq_len(nrow(cluster_union_tbl))) {
  cid  <- cluster_union_tbl$cluster[i]
  regs <- cluster_union_tbl$regulons[[i]]

  freq_dt <- calc_gene_freq_one_cluster(cid, regs, gmt)
  if (is.null(freq_dt)) next

  # Save per-cluster gene frequency table
  fwrite(freq_dt, file = file.path(out_dir, paste0("trendCluster_geneFreq.cluster", cid, ".csv")))

  # Apply core gene selection rule
  if (core_rule_mode == "prop") {
    core_dt   <- freq_dt[prop >= core_min_prop]
    freq_rule <- paste0("prop>=", core_min_prop)
  } else if (core_rule_mode == "k") {
    core_dt   <- freq_dt[freq >= core_min_k]
    freq_rule <- paste0("freq>=", core_min_k)
  } else {
    stop("core_rule_mode must be 'prop' or 'k'")
  }

  # Optional cap on core set size
  if (is.finite(core_max_genes) && nrow(core_dt) > core_max_genes) {
    core_dt   <- core_dt[order(-prop, -freq)][seq_len(core_max_genes)]
    freq_rule <- paste0(freq_rule, ",cap=", core_max_genes)
  }

  freq_all[[as.character(cid)]] <- as.data.frame(freq_dt)

  core_list[[as.character(cid)]] <- list(
    cluster = cid, regulons = regs,
    genes_core = core_dt$gene, n_targets_core = nrow(core_dt),
    n_regulons = length(regs), freq_rule = freq_rule
  )
}

# Assemble core gene set table (keep regulons/genes as list-columns)
core_tbl <- dplyr::bind_rows(lapply(core_list, function(x) {
  tibble::tibble(
    cluster        = x$cluster,
    regulons       = list(x$regulons),
    genes_core     = list(x$genes_core),
    n_targets_core = x$n_targets_core,
    n_regulons     = x$n_regulons,
    freq_rule      = x$freq_rule
  )
}))

if (nrow(core_tbl) == 0) {
  stop("No core gene sets generated. Please relax core thresholds (core_min_prop / core_min_k).")
}

# Join union size for comparison
core_tbl <- core_tbl %>%
  dplyr::left_join(
    cluster_union_tbl %>% dplyr::select(cluster, n_targets_union),
    by = "cluster"
  ) %>%
  dplyr::mutate(
    core_over_union = ifelse(n_targets_union > 0, n_targets_core / n_targets_union, NA_real_)
  )

# Save core gene sets
core_out <- core_tbl %>%
  dplyr::mutate(
    regulons   = sapply(regulons, function(x) paste(x, collapse = ";")),
    genes_core = sapply(genes_core, function(x) paste(x, collapse = ";"))
  )
fwrite(core_out, file = file.path(out_dir, "trendCluster_geneSets.core.csv"))

message("\n[SUMMARY] CORE gene set sizes:")
print(core_tbl %>%
        dplyr::select(cluster, n_regulons, n_targets_union, n_targets_core, core_over_union, freq_rule) %>%
        dplyr::arrange(dplyr::desc(n_targets_core)) %>%
        as.data.frame(), row.names = FALSE)

# ============================================================
# 5) Prepare GO background and SYMBOL -> ENTREZ mapping
#    Background = union of all target genes across all analyzed regulons
# ============================================================
bg_symbols <- unique(unlist(gmt[unique(clu2$regulon)], use.names = FALSE))

all_symbols <- unique(c(
  bg_symbols,
  unlist(cluster_union_tbl$genes_union, use.names = FALSE),
  unlist(core_tbl$genes_core, use.names = FALSE)
))

map <- suppressMessages(
  bitr(all_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
)
sym2ent <- split(map$ENTREZID, map$SYMBOL)

to_entrez <- function(symbols) {
  unique(unlist(sym2ent[symbols], use.names = FALSE))
}
bg_entrez <- to_entrez(bg_symbols)

# ============================================================
# 6) GO enrichment per cluster (BP / MF / CC)
# ============================================================
run_one_cluster <- function(genes_symbol, cluster_id, ont, set_type,
                            n_targets, n_regulons, freq_rule) {
  gene_entrez <- to_entrez(genes_symbol)
  if (length(gene_entrez) < min_gs_go) return(NULL)

  er <- enrichGO(
    gene = gene_entrez, universe = bg_entrez,
    OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
    ont = ont, pvalueCutoff = p_cut, qvalueCutoff = q_cut,
    pAdjustMethod = "BH", readable = TRUE
  )
  if (is.null(er) || nrow(as.data.frame(er)) == 0) return(NULL)

  # Simplify redundant terms (especially for BP)
  er_s <- tryCatch(
    simplify(er, cutoff = simplify_cutoff, by = "p.adjust", select_fun = min),
    error = function(e) er
  )

  df <- as.data.frame(er_s) %>%
    dplyr::mutate(
      cluster = cluster_id, ontology = ont, set_type = set_type,
      n_targets = n_targets, n_regulons = n_regulons, freq_rule = freq_rule
    ) %>%
    dplyr::relocate(cluster, ontology, set_type, n_targets, n_regulons, freq_rule)

  top_df <- df %>%
    dplyr::arrange(p.adjust, dplyr::desc(Count)) %>%
    dplyr::slice_head(n = top_k)

  list(all = df, top = top_df)
}

onts <- c("BP", "MF", "CC")

run_settype_GO <- function(set_type = c("union", "core")) {
  set_type <- match.arg(set_type)

  if (set_type == "union") {
    if (!run_union_GO) return(NULL)
    set_tbl <- cluster_union_tbl %>%
      dplyr::transmute(cluster, genes = genes_union, n_targets = n_targets_union,
                       n_regulons, freq_rule = "union")
  } else {
    if (!run_core_GO) return(NULL)
    set_tbl <- core_tbl %>%
      dplyr::transmute(cluster, genes = genes_core, n_targets = n_targets_core,
                       n_regulons, freq_rule)
  }

  set_tbl <- set_tbl %>% dplyr::filter(n_targets >= min_gs_go)
  message("\n[INFO] GO enrichment for ", set_type, " (clusters passing min_gs=", min_gs_go, "): ", nrow(set_tbl))

  res_all <- list()
  res_top <- list()

  for (ont in onts) {
    message("[INFO] Running GO ", ont, " for ", set_type, " ...")
    all_list <- list()
    top_list <- list()

    for (i in seq_len(nrow(set_tbl))) {
      cid          <- set_tbl$cluster[i]
      genes_symbol <- set_tbl$genes[[i]]

      rr <- run_one_cluster(genes_symbol, cid, ont, set_type,
                            set_tbl$n_targets[i], set_tbl$n_regulons[i], set_tbl$freq_rule[i])
      if (is.null(rr)) next

      all_list[[as.character(cid)]] <- rr$all
      top_list[[as.character(cid)]] <- rr$top

      fwrite(rr$all, file = file.path(out_dir, paste0("GO_", ont, "_", set_type, "_cluster", cid, ".csv")))
    }

    res_all[[ont]] <- bind_rows(all_list)
    res_top[[ont]] <- bind_rows(top_list)

    fwrite(res_all[[ont]], file = file.path(out_dir, paste0("GO_", ont, "_", set_type, "_ALLclusters.full.csv")))
    fwrite(res_top[[ont]], file = file.path(out_dir, paste0("GO_", ont, "_", set_type, "_ALLclusters.top", top_k, ".csv")))
  }

  list(all = res_all, top = res_top)
}

res_union <- run_settype_GO("union")
res_core  <- run_settype_GO("core")

# ============================================================
# 7) Dot plot for top BP terms per cluster
# ============================================================
make_bp_dotplot <- function(top_df, set_type_tag) {
  if (is.null(top_df) || nrow(top_df) == 0) return(NULL)

  bp_top <- top_df %>%
    dplyr::mutate(
      neglog10FDR = -log10(p.adjust),
      GeneRatio_num = sapply(GeneRatio, function(x) {
        parts <- strsplit(x, "/", fixed = TRUE)[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
      }),
      term_label = paste0("C", cluster, ": ", Description)
    )

  p <- ggplot(bp_top, aes(x = factor(cluster),
                           y = forcats::fct_reorder(term_label, neglog10FDR))) +
    geom_point(aes(size = GeneRatio_num, color = neglog10FDR)) +
    scale_color_viridis_c() +
    labs(
      x     = "Trend cluster",
      y     = NULL,
      color = expression(-log[10]("FDR")),
      size  = "GeneRatio",
      title = paste0("GO-BP per trend cluster (", set_type_tag, ", top ", top_k, ", simplified)")
    ) +
    theme_classic(base_size = 12)

  out_pdf <- file.path(out_dir, paste0("GO_BP_", set_type_tag, ".top", top_k, ".dotplot.pdf"))
  ggsave(out_pdf, p, width = 9.5, height = 6.8)
  out_pdf
}

if (!is.null(res_union)) make_bp_dotplot(res_union$top$BP, "union")
if (!is.null(res_core))  make_bp_dotplot(res_core$top$BP, "core")

cat("\nDONE. Outputs in:", out_dir, "\n")
cat("Key files:\n",
    " - trendCluster_geneCounts.summary.csv\n",
    " - trendCluster_geneSets.union.csv / .core.csv\n",
    " - trendCluster_geneFreq.cluster*.csv\n",
    " - GO_{BP,MF,CC}_{union/core}_ALLclusters.full.csv\n",
    " - GO_{BP,MF,CC}_{union/core}_ALLclusters.top", top_k, ".csv\n",
    " - GO_BP_{union/core}.top", top_k, ".dotplot.pdf\n", sep = "")
