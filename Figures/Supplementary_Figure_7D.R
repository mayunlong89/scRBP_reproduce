# ============================================================
# 2026-03-02
# TF–RBP Jaccard comparison
# min-rank aggregation across motifs
# fixed top 1500 targets per protein
# ============================================================

#hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather

 
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(pheatmap)
  library(arrow)
})

# =========================
# 0) Inputs
# =========================

tf_map_file  <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins/66TFs_motifs_links_v9_hgnc.csv"
rbp_map_file <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins/66RBPs_motifs_links.csv"

tf_rank_file <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/SCENIC_python/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"

rbp_rank_files <- list(
  UTR5   = "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/01_gene_level_final_25044motifs/Cluster_Buster_matrix_5UTR_gene_final/motif_gene_rank_formatted.feather",
  UTR3   = "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/01_gene_level_final_25044motifs/Cluster_Buster_matrix_3UTR_gene_final/motif_gene_rank_formatted.feather",
  CDS    = "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/01_gene_level_final_25044motifs/Cluster_Buster_matrix_CDS_gene_final/motif_gene_rank_formatted.feather",
  Intron = "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/01_gene_level_final_25044motifs/Cluster_Buster_matrix_Introns_gene_final/motif_gene_rank_formatted.feather"
)

regions <- names(rbp_rank_files)
top_n   <- 1500

out_csv  <- paste0("66_sharedProteins_Jaccard_RBPregions_vs_TF__10kb_up_and_down_tss.minRankAgg_top", top_n, ".csv")
out_size <- paste0("66_sharedProteins_targetSizes__10kb_up_and_down_tss.minRankAgg_top", top_n, ".csv")
out_pdf  <- paste0("66_sharedProteins_Jaccard_RBPregions_vs_TF__10kb_up_and_down_tss.minRankAgg_top", top_n, ".rowScaled.pdf")

# =========================
# 1) Helper functions
# =========================

msg_stop <- function(...) { message(...); quit(save="no", status=1) }

top_genes_from_rankvec_n <- function(rank_vec, gene_names, top_n = 1500) {
  ok <- !is.na(rank_vec)
  rank_vec   <- rank_vec[ok]
  gene_names <- gene_names[ok]
  if (length(rank_vec) == 0) return(character(0))

  k <- min(top_n, length(rank_vec))
  ord <- order(rank_vec, decreasing = FALSE)
  unique(gene_names[ord][seq_len(k)])
}

aggregate_min_rank <- function(df, motifs, id_col) {
  if (!id_col %in% colnames(df)) msg_stop("Missing id column: ", id_col)

  motifs <- unique(motifs)
  idx <- which(df[[id_col]] %in% motifs)
  if (length(idx) == 0) return(NULL)

  gene_cols <- setdiff(colnames(df), id_col)
  submat <- as.matrix(df[idx, gene_cols, drop = FALSE])
  suppressWarnings(storage.mode(submat) <- "numeric")

  agg <- apply(submat, 2, function(v) suppressWarnings(min(v, na.rm = TRUE)))
  names(agg) <- gene_cols
  agg
}

get_protein_targets_minRank_n <- function(df, motifs, id_col, top_n = 1500) {
  agg <- aggregate_min_rank(df, motifs, id_col)
  if (is.null(agg)) return(character(0))
  top_genes_from_rankvec_n(agg, names(agg), top_n)
}

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

# =========================
# 2) Load motif maps
# =========================

tf_map <- read_csv(tf_map_file, show_col_types = FALSE) %>%
  rename_with(str_trim) %>%
  mutate(TF = str_trim(TF), Motif = str_trim(Motif)) %>%
  filter(TF != "", Motif != "")

rbp_map <- read_csv(rbp_map_file, show_col_types = FALSE) %>%
  rename_with(str_trim) %>%
  mutate(RBP = str_trim(RBP), Motif = str_trim(Motif)) %>%
  filter(RBP != "", Motif != "")

shared_proteins <- intersect(unique(tf_map$TF), unique(rbp_map$RBP))
cat("Shared proteins:", length(shared_proteins), "\n")

if (length(shared_proteins) == 0) msg_stop("No shared proteins found.")

# =========================
# 3) Load ranking matrices
# =========================

cat("Loading TF ranking...\n")
tf_rank <- read_feather(tf_rank_file)
if (!"motifs" %in% colnames(tf_rank)) msg_stop("TF ranking missing 'motifs' column")

rbp_rank <- list()
for (r in regions) {
  cat("Loading RBP ranking:", r, "\n")
  rbp_rank[[r]] <- read_feather(rbp_rank_files[[r]])
  if (!"features" %in% colnames(rbp_rank[[r]]))
    msg_stop("RBP ranking missing 'features' column for region: ", r)
}

# =========================
# 4) Compute Jaccard matrix
# =========================

res <- matrix(
  NA_real_,
  nrow = length(shared_proteins),
  ncol = length(regions),
  dimnames = list(shared_proteins, paste0("RBP_", regions, "_vs_TF"))
)

size_tbl <- tibble(
  Protein = shared_proteins,
  TF_motifs = NA_integer_,
  RBP_motifs = NA_integer_,
  TF_targets = NA_integer_
)

for (r in regions)
  size_tbl[[paste0("RBP_targets_", r)]] <- NA_integer_

for (i in seq_along(shared_proteins)) {

  p <- shared_proteins[i]

  tf_motifs  <- tf_map  %>% filter(TF == p)  %>% pull(Motif)
  rbp_motifs <- rbp_map %>% filter(RBP == p) %>% pull(Motif)

  tf_targets <- get_protein_targets_minRank_n(tf_rank, tf_motifs, "motifs", top_n)

  for (r in regions) {
    rbp_targets_r <- get_protein_targets_minRank_n(rbp_rank[[r]], rbp_motifs, "features", top_n)
    res[p, paste0("RBP_", r, "_vs_TF")] <- jaccard(rbp_targets_r, tf_targets)
    size_tbl[[paste0("RBP_targets_", r)]][i] <- length(rbp_targets_r)
  }

  size_tbl$TF_motifs[i]  <- length(tf_motifs)
  size_tbl$RBP_motifs[i] <- length(rbp_motifs)
  size_tbl$TF_targets[i] <- length(tf_targets)

  cat("[", p, "] TF targets=", length(tf_targets), "\n")
}

# Save outputs
write_csv(as.data.frame(res) %>% rownames_to_column("Protein"), out_csv)
write_csv(size_tbl, out_size)

# =========================
# 5) Heatmap (row scaled)
# =========================

pdf(out_pdf, width = 8, height = 10)

pheatmap(
  res,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 6,
  main = paste0(
    "Jaccard similarity (min-rank aggregation; top ",
    top_n,
    " targets per protein)"
  )
)

dev.off()

cat("Done.\n")
cat("Outputs:\n", out_csv, "\n", out_size, "\n", out_pdf, "\n")



