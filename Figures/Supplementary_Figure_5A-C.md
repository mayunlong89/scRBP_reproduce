# Assessment of Cell Population Size on Regulon Stability

## Overview

This pipeline evaluates how cell population size affects the stability and reproducibility of scRBP-inferred regulons across three tissue types. For each tissue, we subsample cells at seven sizes (2K, 5K, 10K, 25K, 50K, 100K, and a full-size reference) using geometric sketching with five independent random seeds. We then run the full scRBP pipeline (infer → modules → prune → GMT conversion) at each size and compare the resulting regulons against the full-size reference using both RBP-centric and target-centric similarity metrics (F1, Precision, Recall, Jaccard).

**Tissues analyzed (run in parallel):**

| # | Tissue | Total cells | Reference size | Working directory |
|---|--------|-------------|----------------|-------------------|
| 1 | Adult brain | 209,801 | 200K | `.../Palmer_hardwick_short_read/z_04_cell_sizes_regulons` |
| 2 | PBMC | 172,601 | 170K | `.../10X_covid19_PBMC` |
| 3 | Fetal brain | 599,221 | 200K | `.../replication_2_developing_brain_cells` |

**Experimental design:**
- Seeds: `42, 123, 124, 1234, 12345`
- Subsample sizes: `2K, 5K, 10K, 25K, 50K, 100K` + full reference
- mRNA regions: `3UTR, 5UTR, CDS, Introns`
- Total prune runs per tissue: 5 seeds × 7 sizes × 4 regions = **140**
- Grand total across 3 tissues: **420 prune runs**

---

## Tissue-specific configuration

```bash
# ============================================================
# Tissue-specific variables (set ONE block per tissue)
# ============================================================

# --- (1) Adult brain ---
TISSUE="adult_brain"
BASE=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/Palmer_hardwick_short_read
WORK=${BASE}/z_04_cell_sizes_regulons
INPUT_H5AD=""  # Steps 1-2 already completed for adult brain
REF_SIZE="200K"
REF_NCELLS=200000
LABEL="palmer_adult_brain"

# --- (2) PBMC ---
TISSUE="pbmc"
BASE=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_covid19_PBMC
WORK=${BASE}
INPUT_H5AD=${BASE}/PBMC_healthy_subset.h5ad
REF_SIZE="170K"
REF_NCELLS=173684
LABEL="PBMC_healthy"

# --- (3) Fetal brain ---
TISSUE="fetal_brain"
BASE=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells
WORK=${BASE}
INPUT_H5AD=${BASE}/fetal_brain_all_cells.h5ad
REF_SIZE="200K"
REF_NCELLS=200000
LABEL="fetal_brain"
```

---

## Part 1: scRBP Pipeline (Steps 1–5, per seed)

All five steps below are repeated for **each seed** (42, 123, 124, 1234, 12345). The code uses loop variables `${SEED}` and `${SIZE}` throughout.

### Step 1 — Geometric sketching (cell subsampling)

Use `geosketch_general_v2.py` to subsample cells from the full dataset at each target size, preserving the geometric structure of the original population.

```bash
# Subsample cells at each size for a given seed
# Run for each SEED in {42, 123, 124, 1234, 12345}
SEED=42
SIZES=(2000 5000 10000 25000 50000 100000 ${REF_NCELLS})
TOOL_GEOSKETCH=${BASE}  # location of geosketch_general_v2.py

for N in "${SIZES[@]}"; do
  SIZE_LABEL=$(echo "${N}" | awk '{if($1>=1000) printf "%dK", $1/1000; else print $1}')
  OUTPUT=${WORK}/${LABEL}_subset${SIZE_LABEL}cells_seed${SEED}_v2.feather

  srun --pty --mem=200G --cpus-per-task=20 --time=10-00:00:00 \
    python ${TOOL_GEOSKETCH}/geosketch_general_v2.py \
      --input  ${INPUT_H5AD} \
      --output ${OUTPUT} \
      --n_cells ${N} \
      --seed ${SEED}
done
```

> **Note for Adult brain:** Steps 1–2 were previously completed; subsampled GRN files are stored under `z_01_cell_population_subset_analysis/z_geometric/z_01_version/`.


### Step 2 — scRBP infer (GRN inference via GRNBoost2)

Run `scRBP_infer_v6_cor.py` (or v5 for Adult brain) to infer RBP-target regulatory networks from the subsampled expression matrices.

```bash
# GRN inference for each cell size
TOOL_INFER=${BASE}/../Palmer_hardwick_short_read  # location of scRBP_infer scripts
RBP_LIST=${TOOL_INFER}/rbp_list_616RBPs.tsv

for N in "${SIZES[@]}"; do
  SIZE_LABEL=$(echo "${N}" | awk '{if($1>=1000) printf "%dK", $1/1000; else print $1}')
  INPUT_MATRIX=${WORK}/${LABEL}_subset${SIZE_LABEL}cells_seed${SEED}_v2.feather
  OUTPUT_GRN=${WORK}/z_GRNBoost2_${LABEL}_subset${SIZE_LABEL}cells_seed${SEED}_v2_grn_modules.tsv

  srun --mem=150G --time=10-00:00:00 --cpus-per-task=20 --pty \
    python ${TOOL_INFER}/scRBP_infer_v6_cor.py \
      --matrix   ${INPUT_MATRIX} \
      --rbp_list ${RBP_LIST} \
      --output   ${OUTPUT_GRN} \
      --n_workers 10 \
      --batch_size 5 \
      --threshold 0.03 \
      --method grnboost2 \
      --log ${OUTPUT_GRN%.tsv}.log
done
```


### Step 3 — Extract modules (6 strategies)

Apply `scRBP_merge_grn_modules.py` to extract RBP regulatory modules from each GRN network using 6 merging strategies (combinations of top-N filtering and percentile thresholds).

```bash
TOOL_MODULE=${BASE}/../Palmer_hardwick_short_read/z_03_regulon_construction
SEED_DIR=${WORK}/seed${SEED}_results  # output directory per seed

for N in "${SIZES[@]}"; do
  SIZE_LABEL=$(echo "${N}" | awk '{if($1>=1000) printf "%dK", $1/1000; else print $1}')

  # Input GRN file path (tissue-specific naming)
  # Adult brain: stored in z_01_cell_population_subset_analysis (for <200K) or main dir (200K)
  # PBMC / Fetal brain: stored in WORK directory
  INPUT_GRN=${WORK}/z_GRNBoost2_${LABEL}_subset${SIZE_LABEL}cells_seed${SEED}*_modules.tsv
  OUTPUT_MODULE=${SEED_DIR}/z_GRNBoost2_${LABEL}_subset${SIZE_LABEL}cells_seed${SEED}_modules.tsv

  python ${TOOL_MODULE}/scRBP_merge_grn_modules.py \
    --input ${INPUT_GRN} \
    --importance_threshold 0.005 \
    --top_n_list 5,10,50 \
    --target_top_n 50 \
    --percentile 0.75,0.90 \
    --output_merged ${OUTPUT_MODULE} \
    --verbose
done
```


### Step 4 — scRBP prune (motif-based regulon refinement)

Run `scRBP_prune_ctxscore_mt_v4.py` across all 4 mRNA regions for each cell size. This is the most computationally intensive step — each run requires ~50GB memory and 8+ cores.

```bash
# Shared tool and annotation paths
TOOL_PRUNE=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/01_gene_level_final_25044motifs/Cluster_Buster_matrix_3UTR_gene_final
ANNO=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/01_gene_level_final_25044motifs
MOTIF_LINKS=${ANNO}/11_616RBPs_20746motifs_links_annotation.csv

REGIONS=(3UTR 5UTR CDS Introns)

for N in "${SIZES[@]}"; do
  SIZE_LABEL=$(echo "${N}" | awk '{if($1>=1000) printf "%dK", $1/1000; else print $1}')
  MODULE_TSV=${SEED_DIR}/z_GRNBoost2_${LABEL}_subset${SIZE_LABEL}cells_seed${SEED}_modules.tsv
  SIZE_RESULT_DIR=${SEED_DIR}/${SIZE_LABEL}_results

  for REGION in "${REGIONS[@]}"; do
    RANK_DATA=${ANNO}/Cluster_Buster_matrix_${REGION}_gene_final
    SAVE_DIR=${SIZE_RESULT_DIR}/Results_final_multiprocessing_all6strategies_${SIZE_LABEL}Cells_seed${SEED}_top1500_${REGION}

    srun --mem=50G --ntasks=1 --cpus-per-task=32 --time=10-00:00:00 \
      python ${TOOL_PRUNE}/scRBP_prune_ctxscore_mt_v4.py \
        --rbp_targets      ${MODULE_TSV} \
        --motif_rbp_links  ${MOTIF_LINKS} \
        --motif_target_ranks ${RANK_DATA}/motif_gene_rank_formatted.feather \
        --save_dir         ${SAVE_DIR}/ \
        --rank_threshold   1500 \
        --auc_threshold    0.05 \
        --min_genes        20 \
        --nes_threshold    3.0 \
        --n_jobs           8 \
        --chunksize        1
  done
done
```

> **Per-tissue total runs for Step 4:**
> 5 seeds × 7 sizes × 4 regions = **140 jobs** (each ~50GB, 8 CPUs)


### Step 5 — Convert regulons to GMT format

Convert the pruned regulon CSV output to GMT format (both gene-symbol and Entrez-ID versions) for downstream pathway analysis.

```bash
# Run in each prune output directory containing scRBP_prune_matched_results.csv
TOOL_CONVERT=${ANNO}/Cluster_Buster_matrix_3UTR_gene_final
MAGMA_ANNO=/mnt/isilon/gandal_lab/mayl/01_GWAS_tools/MAGMA

# Apply to each seed × size × region output directory
cd ${SAVE_DIR}  # the prune output directory

python ${TOOL_CONVERT}/scRBP_convert_format_v2.py \
  --input      scRBP_prune_matched_results.csv \
  --out-symbol scRBP_prune_matched_results_min5genes.symbol_default.gmt \
  --out-entrez scRBP_prune_matched_results_min5genes.entrez_default.gmt \
  --map-custom ${MAGMA_ANNO}/NCBI38.gene.loc \
  --min_genes 5 \
  --sanitize-rbp \
  --drop-unmapped-genes \
  --drop-empty-sets
```

> **Gene mapping rate:** ~96.5% (7,301/7,563 genes mapped via NCBI38.gene.loc)


---

## Part 2: Regulon Stability Evaluation

After running all seeds × sizes × regions, we evaluate how well smaller cell populations recover the regulons inferred from the full-size reference. The evaluation code below is **tissue-agnostic** — only the `root` path and `REF_SIZE` variable differ across tissues.

### Configuration (tissue-specific)

```r
# Set ONE of the following per tissue:

# (1) Adult brain
root <- "/mnt/isilon/.../Palmer_hardwick_short_read/z_04_cell_sizes_regulons"
REF_SIZE <- "200K"
OUTPUT_PREFIX <- "adult_brain"

# (2) PBMC
root <- "/mnt/isilon/.../10X_covid19_PBMC"
REF_SIZE <- "170K"
OUTPUT_PREFIX <- "PBMC"

# (3) Fetal brain
root <- "/mnt/isilon/.../replication_2_developing_brain_cells"
REF_SIZE <- "200K"
OUTPUT_PREFIX <- "fetal_brain"
```


### 2.1 RBP-centric evaluation (seed-to-seed comparison)

For each cell size, compare the set of RBPs recovered in each seed against the corresponding reference-size result from the same seed. Metrics: F1, Precision, Recall, Jaccard over the RBP identity sets.

```r
rm(list = ls()); library(tidyverse)

sizes_all <- c("2K", "5K", "10K", "25K", "50K", "100K", REF_SIZE)
regions   <- c("5UTR", "3UTR", "CDS", "Introns")
gmt_name  <- "scRBP_prune_matched_results_min5genes.symbol_default.gmt"

# Region pattern matching (case-insensitive)
region_pat <- list(
  "5UTR"    = "(?i)(^|/)(5UTR|UTR5|5[_-]?prime[_-]?UTR)(/|$)",
  "3UTR"    = "(?i)(^|/)(3UTR|UTR3|3[_-]?prime[_-]?UTR)(/|$)",
  "CDS"     = "(?i)(^|/)(CDS|coding)(/|$)",
  "Introns" = "(?i)(^|/)(Introns?|intron)(/|$)"
)

# --- Utility functions ---
read_gmt_rbps <- function(gmt_file) {
  if (!file.exists(gmt_file)) return(character())
  lines <- readr::read_lines(gmt_file)
  if (!length(lines)) return(character())
  lines |> strsplit("\t", fixed = TRUE) |> purrr::map_chr(~ .x[[1]]) |> unique()
}

size_lab <- function(p) sub("_results$", "", basename(p))

find_region_gmt <- function(size_dir, region) {
  files <- list.files(size_dir, pattern = paste0("^", gmt_name, "$"),
                      recursive = TRUE, full.names = TRUE)
  if (!length(files)) return(NA_character_)
  pat  <- region_pat[[region]]
  cand <- files[grepl(pat, dirname(files), perl = TRUE)]
  if (!length(cand)) cand <- files[grepl(region, files, ignore.case = TRUE)]
  if (!length(cand)) return(NA_character_)
  cand[1]
}

pair_metrics <- function(A, B) {
  A <- unique(A); B <- unique(B)
  inter <- length(intersect(A, B))
  union <- length(union(A, B))
  p  <- if (!length(A)) NA_real_ else inter / length(A)
  r  <- if (!length(B)) NA_real_ else inter / length(B)
  f1 <- if (is.na(p) || is.na(r) || (p + r) == 0) NA_real_ else 2 * p * r / (p + r)
  j  <- if (!union) NA_real_ else inter / union
  tibble(precision = p, recall = r, f1 = f1, jaccard = j,
         n_test = length(A), n_ref = length(B), n_inter = inter)
}

# --- Collect metrics across all seeds ---
seed_to_long_df <- function(seed_dir) {
  seed_id   <- sub("_results$", "", basename(seed_dir))
  size_dirs <- list.dirs(seed_dir, recursive = FALSE, full.names = TRUE)
  size_dirs <- size_dirs[dir.exists(size_dirs) & str_detect(basename(size_dirs), "K_results$")]
  size_map  <- setNames(size_dirs, size_lab(size_dirs))
  if (!REF_SIZE %in% names(size_map)) return(tibble())

  rows <- list()
  for (rg in regions) {
    ref_gmt <- find_region_gmt(size_map[[REF_SIZE]], rg)
    if (is.na(ref_gmt)) next
    ref_set <- read_gmt_rbps(ref_gmt)
    for (sz in setdiff(intersect(names(size_map), sizes_all), REF_SIZE)) {
      tst_gmt <- find_region_gmt(size_map[[sz]], rg)
      if (is.na(tst_gmt)) next
      tst_set <- read_gmt_rbps(tst_gmt)
      rows[[length(rows) + 1]] <-
        pair_metrics(tst_set, ref_set) |>
        mutate(seed = seed_id, region = rg, size = sz)
    }
  }
  bind_rows(rows) |>
    pivot_longer(c(precision, recall, f1, jaccard),
                 names_to = "metric", values_to = "score") |>
    mutate(
      size   = factor(size, levels = c("2K","5K","10K","25K","50K","100K")),
      region = factor(region, levels = regions),
      metric = factor(metric,
                      levels = c("f1","jaccard","precision","recall"),
                      labels = c("F1 score","Jaccard","Precision","Recall")),
      seed   = factor(seed_id)
    )
}

seed_dirs <- dir(root, pattern = "^seed\\d+_results$", full.names = TRUE, recursive = FALSE)
long_all  <- bind_rows(lapply(seed_dirs, seed_to_long_df))

# --- Plot: boxplot + jitter ---
set.seed(2024)
p <- ggplot(long_all, aes(x = size, y = score)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, color = "black", fill = "white") +
  geom_point(aes(color = region, shape = seed),
             position = position_jitter(width = 0.32), size = 2, alpha = 0.85) +
  facet_wrap(~ metric, ncol = 2) +
  scale_y_continuous(limits = c(0.25, 0.75)) +
  scale_color_brewer(palette = "Set2", name = "mRNA region") +
  labs(x = "Cell count", y = "Score") +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

ggsave(paste0("rbp_similarity_seed_to_seed_", OUTPUT_PREFIX, ".pdf"), p, width = 8.5, height = 6.5)
```


### 2.2 RBP-centric evaluation (consensus reference)

Construct a **consensus reference** by majority vote (≥3/5 seeds) at the full reference size, then evaluate each seed × size against this consensus.

```r
# --- Build consensus reference from all seeds at reference size ---
get_seed_ref_sets <- function(seed_dir) {
  seed_id   <- sub("_results$", "", basename(seed_dir))
  size_dirs <- list.dirs(seed_dir, recursive = FALSE, full.names = TRUE)
  size_dirs <- size_dirs[dir.exists(size_dirs) & str_detect(basename(size_dirs), "K_results$")]
  size_map  <- setNames(size_dirs, size_lab(size_dirs))
  if (!REF_SIZE %in% names(size_map)) return(NULL)
  tibble(seed = seed_id, region = regions) |>
    mutate(gmt  = map_chr(region, ~ find_region_gmt(size_map[[REF_SIZE]], .x)),
           rbps = map(gmt, read_gmt_rbps)) |>
    select(seed, region, rbps)
}

seed_ref <- bind_rows(compact(map(seed_dirs, get_seed_ref_sets)))

# Majority vote: RBPs present in ≥3 of 5 seeds
consensus_ref <- seed_ref %>%
  group_by(region) %>%
  summarise(rbps_ref = list({ tb <- table(unlist(rbps)); names(tb)[tb >= 3] }),
            .groups = "drop")
consensus_list <- setNames(consensus_ref$rbps_ref, consensus_ref$region)

# --- Evaluate 2K–100K against consensus ---
sizes_eval <- c("2K", "5K", "10K", "25K", "50K", "100K")
rows <- list()
for (sd in seed_dirs) {
  seed_id   <- sub("_results$", "", basename(sd))
  size_dirs <- list.dirs(sd, recursive = FALSE, full.names = TRUE)
  size_dirs <- size_dirs[dir.exists(size_dirs) & str_detect(basename(size_dirs), "K_results$")]
  size_map  <- setNames(size_dirs, size_lab(size_dirs))
  for (sz in intersect(names(size_map), sizes_eval)) {
    for (rg in regions) {
      gmt <- find_region_gmt(size_map[[sz]], rg)
      if (is.na(gmt)) next
      tst_set <- read_gmt_rbps(gmt)
      ref_set <- consensus_list[[rg]]
      rows[[length(rows) + 1]] <-
        pair_metrics(tst_set, ref_set) |>
        mutate(seed = seed_id, region = rg, size = sz)
    }
  }
}
metrics_df <- bind_rows(rows)

long_all <- metrics_df |>
  pivot_longer(c(precision, recall, f1, jaccard),
               names_to = "metric", values_to = "score") |>
  mutate(
    size   = factor(size, levels = sizes_eval),
    region = factor(region, levels = regions),
    metric = factor(metric,
                    levels = c("f1","jaccard","precision","recall"),
                    labels = c("F1 score","Jaccard","Precision","Recall")),
    seed   = factor(seed)
  )

# --- Compute split-half P5 threshold from reference-size pairwise comparisons ---
pairs  <- combn(unique(seed_ref$seed), 2, simplify = FALSE)
rows2  <- list()
for (pp in pairs) {
  for (rg in regions) {
    A <- seed_ref %>% filter(seed == pp[1], region == rg) %>% pull(rbps) %>% purrr::pluck(1)
    B <- seed_ref %>% filter(seed == pp[2], region == rg) %>% pull(rbps) %>% purrr::pluck(1)
    if (is.null(A) || is.null(B)) next
    rows2[[length(rows2) + 1]] <- pair_metrics(A, B) %>% mutate(region = rg)
  }
}
long_ref <- bind_rows(rows2) %>%
  pivot_longer(c(precision, recall, f1, jaccard),
               names_to = "metric", values_to = "score") %>%
  mutate(metric = factor(metric,
                         levels = c("f1","jaccard","precision","recall"),
                         labels = c("F1 score","Jaccard","Precision","Recall")))

# Median P5 threshold per metric and threshold band per region
theta_median <- long_ref %>%
  group_by(region, metric) %>%
  summarise(p5 = quantile(score, 0.05, na.rm = TRUE), .groups = "drop") %>%
  group_by(metric) %>%
  summarise(p5 = median(p5), .groups = "drop")

band_df <- long_ref %>%
  group_by(region, metric) %>%
  summarise(p5 = quantile(score, 0.05, na.rm = TRUE), .groups = "drop") %>%
  group_by(metric) %>%
  summarise(ymin = min(p5), ymax = max(p5), .groups = "drop")

# --- Plot with threshold band + reference line ---
set.seed(2024)
p <- ggplot(long_all, aes(x = size, y = score)) +
  geom_rect(data = band_df,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#CFE8FF", alpha = 0.35) +
  geom_boxplot(width = 0.65, outlier.shape = NA, color = "black", fill = "white") +
  geom_point(aes(color = region, shape = seed),
             position = position_jitter(width = 0.32), size = 2, alpha = 0.85) +
  geom_hline(data = theta_median, aes(yintercept = p5),
             linetype = "dashed", linewidth = 0.6, color = "#f68989") +
  facet_wrap(~ metric, ncol = 2) +
  scale_y_continuous(limits = c(0, 0.8)) +
  scale_color_brewer(palette = "Set2", name = "mRNA region") +
  labs(x = "Cell count", y = "Score",
       title = paste0("Against ", REF_SIZE, " consensus reference (RBP-centric)")) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

ggsave(paste0("rbp_similarity_vs", REF_SIZE, "_consensus_", OUTPUT_PREFIX, ".pdf"),
       p, width = 8.5, height = 6.5)
```


### 2.3 Ceiling normalization (optional)

Normalize raw scores by the split-half reproducibility ceiling at the reference size, so that scores approaching 1.0 indicate performance as good as the intrinsic seed-to-seed variability.

```r
CEIL_QUANT <- 0.50  # median ceiling

ceil_pooled <- long_ref %>%
  group_by(metric) %>%
  summarise(ceil = quantile(score, CEIL_QUANT, na.rm = TRUE), .groups = "drop")

long_rel <- long_all %>%
  left_join(ceil_pooled, by = "metric") %>%
  mutate(score_rel = ifelse(ceil > 0, pmin(score / ceil, 1), NA_real_))

p_rel <- ggplot(long_rel, aes(x = size, y = score_rel)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, color = "black", fill = "white") +
  geom_point(aes(color = region, shape = seed),
             position = position_jitter(width = 0.32), size = 2, alpha = 0.85) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
  facet_wrap(~ metric, ncol = 2) +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_color_brewer(palette = "Set2", name = "mRNA region") +
  labs(x = "Cell count", y = "Relative score",
       title = "Ceiling-normalized vs split-half ceiling") +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

ggsave(paste0("rbp_relative_vs", REF_SIZE, "_ceiling_", OUTPUT_PREFIX, ".pdf"),
       p_rel, width = 8.8, height = 6.6)
```


### 2.4 Target-centric evaluation (per-RBP regulon overlap)

Instead of comparing RBP identity sets, compare the **target gene sets** for each shared RBP between the test size and the reference. This evaluates whether regulon content (not just existence) is preserved.

```r
# Read GMT as RBP -> gene list mapping
read_gmt_map <- function(gmt_file) {
  if (!file.exists(gmt_file)) return(list())
  lines <- readr::read_lines(gmt_file)
  if (!length(lines)) return(list())
  toks  <- strsplit(lines, "\t", fixed = TRUE)
  rbps  <- vapply(toks, function(x) x[[1]], character(1))
  genes <- lapply(toks, function(x) {
    if (length(x) <= 2) character() else unique(x[3:length(x)])
  })
  setNames(genes[vapply(genes, length, integer(1)) > 0], rbps[vapply(genes, length, integer(1)) > 0])
}

# Compute per-RBP target metrics across all seeds
seed_to_target_df <- function(seed_dir) {
  seed_id   <- sub("_results$", "", basename(seed_dir))
  size_dirs <- list.dirs(seed_dir, recursive = FALSE, full.names = TRUE)
  size_dirs <- size_dirs[dir.exists(size_dirs) & str_detect(basename(size_dirs), "K_results$")]
  size_map  <- setNames(size_dirs, size_lab(size_dirs))
  if (!REF_SIZE %in% names(size_map)) return(tibble())

  rows <- list()
  for (rg in regions) {
    ref_gmt <- find_region_gmt(size_map[[REF_SIZE]], rg)
    if (is.na(ref_gmt)) next
    ref_map <- read_gmt_map(ref_gmt)
    for (sz in setdiff(intersect(names(size_map), sizes_all), REF_SIZE)) {
      tst_gmt <- find_region_gmt(size_map[[sz]], rg)
      if (is.na(tst_gmt)) next
      tst_map <- read_gmt_map(tst_gmt)
      shared  <- intersect(names(tst_map), names(ref_map))
      for (r in shared) {
        rows[[length(rows) + 1]] <-
          pair_metrics(unique(tst_map[[r]]), unique(ref_map[[r]])) |>
          mutate(RBP = r, seed = seed_id, region = rg, size = sz)
      }
    }
  }
  bind_rows(rows) |>
    pivot_longer(c(precision, recall, f1, jaccard),
                 names_to = "metric", values_to = "score") |>
    mutate(
      size   = factor(size, levels = c("2K","5K","10K","25K","50K","100K")),
      region = factor(region, levels = regions),
      metric = factor(metric,
                      levels = c("f1","jaccard","precision","recall"),
                      labels = c("F1 score","Jaccard","Precision","Recall")),
      seed   = factor(seed_id)
    )
}

long_target <- bind_rows(lapply(seed_dirs, seed_to_target_df))
```


### 2.5 Target-centric: average across seeds

Average the per-RBP target-centric scores across the 5 seeds to reduce noise, then plot.

```r
long_avg <- long_target %>%
  filter(!is.na(score)) %>%
  group_by(RBP, region, size, metric) %>%
  summarise(score   = mean(score, na.rm = TRUE),
            n_seeds = n_distinct(seed), .groups = "drop") %>%
  filter(n_seeds >= 1) %>%
  mutate(
    size   = factor(size, levels = c("2K","5K","10K","25K","50K","100K")),
    region = factor(region, levels = regions),
    metric = factor(metric, levels = c("F1 score","Jaccard","Precision","Recall"))
  )

p_avg <- ggplot(long_avg, aes(x = size, y = score)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, color = "black", fill = "white") +
  geom_point(aes(color = region),
             position = position_jitter(width = 0.22), size = 1.2, alpha = 0.65) +
  facet_wrap(~ metric, ncol = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = "Set2", name = "mRNA region") +
  labs(x = "Cell count", y = "Score",
       subtitle = "Per-RBP target-centric (averaged across seeds)") +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

ggsave(paste0("target_centric_avg_", OUTPUT_PREFIX, ".pdf"), p_avg, width = 9, height = 6.5)
```


### 2.6 Target-centric: top-100 RBPs by mean F1 (main result)

Select the top-100 RBPs per cell size (ranked by mean F1 across 4 mRNA regions) and plot their target-level similarity metrics.

```r
top_k <- 100
rank_metric <- "F1 score"

# Rank RBPs by mean F1 at each cell size
rank_tbl <- long_avg %>%
  filter(metric == rank_metric) %>%
  group_by(size, RBP) %>%
  summarise(rank_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  group_by(size) %>%
  slice_max(order_by = rank_score, n = top_k, with_ties = FALSE) %>%
  ungroup()

long_top <- long_avg %>%
  semi_join(rank_tbl %>% select(size, RBP), by = c("size", "RBP"))

p_top <- ggplot(long_top, aes(x = size, y = score)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, color = "black", fill = "white") +
  geom_point(aes(color = region),
             position = position_jitter(width = 0.22), size = 1.2, alpha = 0.65) +
  facet_wrap(~ metric, ncol = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = "Set2", name = "mRNA region") +
  labs(x = "Cell count", y = "Score",
       subtitle = paste0("Top ", top_k, " RBPs by mean ", rank_metric, " (seed-averaged)")) +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

ggsave(paste0("target_centric_top100_", OUTPUT_PREFIX, ".pdf"), p_top, width = 9, height = 6.5)
```
