# F1 / Precision / Recall Benchmarking Pipeline (Isoform-Centric)

> **Date:** 2025-06-21 / 2025-06-24  
> **Author:** Yunlong Ma

---

## Table of Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Step 1.0 — Top-k threshold sensitivity test (Clustered motifs only)](#step-10--top-k-threshold-sensitivity-test)
4. [Step 1.1–1.5 — F1 evaluation across 5 methods × 4 regions](#step-1115--f1-evaluation-across-5-methods--4-regions)
5. [Step 2 — Visualization: per-region comparison across methods](#step-2--visualization-per-region-comparison)
6. [Step 3 — Visualization: best-region performance comparison](#step-3--visualization-best-region-performance-comparison)

---

## Overview

### Approach: Target-Centric Evaluation

Each RBP has ~25–50 motifs, and each motif produces a transcriptome-wide target gene ranking. For evaluation:

1. Extract eCLIP gold-standard target genes for each RBP as **positives**.
2. Randomly sample an equal number of background genes (1:1 ratio) as **negatives**.
3. Using the existing motif-based ranking, predict the **top-k** ranked genes as predicted positives.
4. Compute **F1**, **Precision**, and **Recall**.
5. Repeat N times with bootstrap resampling; report the mean.
6. For each RBP, retain only the **best-performing motif** (highest F1).

### Design Matrix

| Dimension | Values |
|-----------|--------|
| **Methods** | Clustered, Singleton, FIMO, HOMER2, Archetype |
| **Regions** | 3UTR, 5UTR, CDS, Introns |
| **Cell lines** | K562 (120 RBPs), HepG2 (103 RBPs) |
| **Top-k cutoffs** | 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000 (Step 1.0 only) |
| **Final top-k** | 2000 (Steps 1.1–1.5) |

**Total runs (Steps 1.1–1.5):** 5 methods × 4 regions × 2 cell lines = **40 scripts**

---

## Prerequisites

### Python dependencies

```bash
pip install pandas numpy scikit-learn matplotlib seaborn tqdm pyarrow
```

### Input files

| File type | Description |
|-----------|-------------|
| `motif_sorted_rank_*.feather` | Pre-computed motif-by-gene rank matrix (rows = motifs, cols = genes) |
| `encode_rbp_motif_link_clustered_motifs_only.csv` | Ground-truth RBP–motif mapping (for Clustered/Singleton/FIMO/HOMER) |
| `encode_rbp_archetype_motif_links_only.csv` | Ground-truth RBP–motif mapping (for Archetype) |
| `ENCODE_K562_120RBPs_targets_isoforms.csv` | eCLIP gold-standard targets for K562 |
| `ENCODE_HepG2_103RBPs_targets_isoforms.csv` | eCLIP gold-standard targets for HepG2 |

### Rank matrix naming convention

| Method | Rank matrix pattern |
|--------|-------------------|
| Clustered | `motif_sorted_rank_8578clustered_motifs_{REGION}_subset_3866motifs_isoforms.feather` |
| Singleton | `motif_sorted_rank_8578singleton_motifs_{REGION}_subset_3866motifs_isoforms.feather` |
| Archetype | `motif_sorted_rank_8578clustered_motifs_{REGION}_subset_1050archetypemotifs_isoform.feather` |
| FIMO | `motif_sorted_rank_3866clustered_motifs_FIMO_{REGION}.feather` |
| HOMER | `motif_sorted_rank_3866clustered_motifs_HOMER_{REGION}.feather` |

---

## Step 1.0 — Top-k threshold sensitivity test

Test different top-k cutoffs (500–5000) using **Clustered motifs on 3UTR only** to determine the optimal threshold.

### Consolidated Python script: `run_F1_topk_sensitivity.py`

```python
#!/usr/bin/env python3
"""
Top-k sensitivity test for F1/Precision/Recall evaluation.
Tests multiple top-k cutoffs on Clustered motifs, 3UTR region.
"""
import pandas as pd
import numpy as np
from sklearn.metrics import precision_score, recall_score, f1_score
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import sys

# ============================================================
# Configuration
# ============================================================
BASE = "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE"
ISO_BASE = f"{BASE}_isoforms"

RANK_MATRIX_PATH = f"{ISO_BASE}/motif_sorted_rank_8578clustered_motifs_3UTR_subset_3866motifs_isoforms.feather"
MOTIF_RBP_MAP_PATH = f"{BASE}/encode_rbp_motif_link_clustered_motifs_only.csv"

ECLIP_FILES = {
    "K562":  f"{ISO_BASE}/ENCODE_K562_120RBPs_targets_isoforms.csv",
    "HepG2": f"{ISO_BASE}/ENCODE_HepG2_103RBPs_targets_isoforms.csv",
}

TOP_K_VALUES = [500, 1000, 1500, 2000, 2500, 3000, 4000, 5000]
RATIO = 1       # Positive-to-negative ratio
N_REPEATS = 1   # Bootstrap repeats (1 for quick sensitivity scan)

# ============================================================
# Core evaluation function
# ============================================================
def evaluate_motif_predictions(rank_matrix, motif_rbp_map, eclip_df,
                                max_k, ratio=1, N=1):
    """Compute F1/Precision/Recall for each motif at a given top-k cutoff."""
    rbp_to_targets = eclip_df.groupby("RBP")["Gene"].apply(set).to_dict()
    motif_to_rbp = motif_rbp_map.set_index("Motif")["RBP"].to_dict()

    # Build motif -> gene rank dict (lower rank = higher priority)
    motif_rank_dict = {
        motif: {gene: rank for rank, gene in enumerate(row.dropna().values)}
        for motif, row in rank_matrix.iterrows()
    }

    results = []
    for motif, gene_rank in tqdm(motif_rank_dict.items(), desc=f"top-{max_k}"):
        if motif not in motif_to_rbp:
            continue
        rbp = motif_to_rbp[motif]
        if rbp not in rbp_to_targets:
            continue

        gold_targets = rbp_to_targets[rbp]
        motif_genes = set(gene_rank.keys())
        positives = list(motif_genes & gold_targets)

        if len(positives) < 10:
            continue

        negatives_pool = list(motif_genes - set(positives))
        if len(negatives_pool) < ratio * len(positives):
            continue

        # Bootstrap evaluation
        f1_list, precision_list, recall_list = [], [], []
        for _ in range(N):
            negatives = np.random.choice(negatives_pool, size=ratio * len(positives), replace=False)
            test_genes = positives + list(negatives)
            test_ranks = [(gene, gene_rank[gene]) for gene in test_genes]
            test_ranks.sort(key=lambda x: x[1])

            top_k = max_k
            predicted_positives = set([g for g, _ in test_ranks[:top_k]])

            y_true = [1 if gene in positives else 0 for gene in test_genes]
            y_pred = [1 if gene in predicted_positives else 0 for gene in test_genes]

            precision_list.append(precision_score(y_true, y_pred, zero_division=0))
            recall_list.append(recall_score(y_true, y_pred, zero_division=0))
            f1_list.append(f1_score(y_true, y_pred, zero_division=0))

        results.append({
            "Motif": motif, "RBP": rbp,
            "Positive_Genes": len(positives),
            "F1": np.mean(f1_list),
            "Precision": np.mean(precision_list),
            "Recall": np.mean(recall_list),
        })

    return pd.DataFrame(results)


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    print("Loading rank matrix...")
    rank_matrix = pd.read_feather(RANK_MATRIX_PATH)
    motif_rbp_map = pd.read_csv(MOTIF_RBP_MAP_PATH)

    for cell_line, eclip_path in ECLIP_FILES.items():
        eclip_df = pd.read_csv(eclip_path)

        for max_k in TOP_K_VALUES:
            print(f"\n=== {cell_line} | top-{max_k} ===")

            result_df = evaluate_motif_predictions(
                rank_matrix, motif_rbp_map, eclip_df,
                max_k=max_k, ratio=RATIO, N=N_REPEATS
            )

            out_csv = f"motif_eval_boxplot_{cell_line}_top{max_k}_ratio{RATIO}_repeat{N_REPEATS}times_clustered_3UTR_isoforms.csv"
            result_df.to_csv(out_csv, index=False)

            # Keep best motif per RBP and plot
            best = result_df.sort_values("F1", ascending=False).groupby("RBP").head(1)
            plot_df = best.melt(
                id_vars=["Motif", "RBP"],
                value_vars=["F1", "Precision", "Recall"],
                var_name="Metric", value_name="Score"
            )

            plt.figure(figsize=(8, 6))
            sns.boxplot(x="Metric", y="Score", data=plot_df, palette="Set2")
            sns.stripplot(x="Metric", y="Score", data=plot_df, color="black", size=2, jitter=True, alpha=0.3)
            plt.title(f"Motif Prediction Evaluation (Top-{max_k}, {cell_line})")
            plt.tight_layout()
            plt.savefig(out_csv.replace(".csv", ".pdf"))
            plt.close()

            print(f"  Saved: {out_csv}")

    print("\nDone!")
```

### Run

```bash
srun --mem=100G --time=20-00:00:00 --pty python run_F1_topk_sensitivity.py
```

### Top-k sensitivity visualization (R)

```r
# Visualize how F1/Precision/Recall change across different top-k thresholds
library(dplyr); library(ggplot2); library(readr); library(tidyr); library(RColorBrewer)

files <- list.files(path = "./", pattern = "*.csv", full.names = TRUE)

# Batch read and annotate source
all_data <- lapply(files, function(file) {
  df <- read_csv(file)
  parts <- strsplit(basename(file), "_")[[1]]
  df$Method <- gsub(".csv", "", parts[5])
  df$CellLine <- parts[4]
  return(df)
}) %>% bind_rows()

# Retain best motif per RBP
best_motifs <- all_data %>%
  group_by(RBP, CellLine, Method) %>%
  slice_max(order_by = F1, n = 1, with_ties = FALSE) %>%
  ungroup()

# Convert to long format
plot_data <- best_motifs %>%
  select(Method, CellLine, F1, Precision, Recall) %>%
  pivot_longer(cols = c("F1", "Precision", "Recall"),
               names_to = "Metric", values_to = "Score")

plot_data$Method <- factor(plot_data$Method,
                           levels = c("top500","top1000","top1500","top2000",
                                      "top2500","top3000","top4000","top5000"))

# Boxplot with jitter
p <- ggplot(plot_data, aes(x = Method, y = Score, fill = Method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.4, color = "black") +
  facet_wrap(~Metric, scales = "free_y", nrow = 3) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14, face = "bold")) +
  labs(title = "Top-k Sensitivity Analysis (Clustered, 3UTR)",
       x = "Top-k Cutoff", y = "Score")

pdf("topK_sensitivity_boxplot_3UTR_isoforms.pdf", width = 10, height = 8)
print(p)
dev.off()
```

---

## Step 1.1–1.5 — F1 evaluation across 5 methods × 4 regions

Using the selected **top-k = 2000** cutoff, evaluate all 5 methods across all 4 mRNA regions and both cell lines.

### Consolidated Python script: `run_F1_all_methods.py`

```python
#!/usr/bin/env python3
"""
F1/Precision/Recall evaluation for 5 motif methods × 4 mRNA regions × 2 cell lines.
Consolidates Sections 1.1–1.5 into a single parameterized script.
"""
import pandas as pd
import numpy as np
from sklearn.metrics import precision_score, recall_score, f1_score
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import os

# ============================================================
# Configuration
# ============================================================
BASE = "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE"
ISO_BASE = f"{BASE}_isoforms"

MAX_K = 2000      # Top-k cutoff (selected from Step 1.0 sensitivity analysis)
RATIO = 1         # Positive-to-negative ratio
N_REPEATS = 30    # Bootstrap repeats

REGIONS = ["3UTR", "5UTR", "CDS", "Introns"]

CELL_LINES = {
    "K562":  f"{ISO_BASE}/ENCODE_K562_120RBPs_targets_isoforms.csv",
    "HepG2": f"{ISO_BASE}/ENCODE_HepG2_103RBPs_targets_isoforms.csv",
}

# Method configurations: method_tag -> (rank_matrix_pattern, motif_rbp_map)
METHODS = {
    "clustered": {
        "rank_pattern": f"{ISO_BASE}/motif_sorted_rank_8578clustered_motifs_{{REGION}}_subset_3866motifs_isoforms.feather",
        "motif_map": f"{BASE}/encode_rbp_motif_link_clustered_motifs_only.csv",
    },
    "singleton": {
        "rank_pattern": f"{ISO_BASE}/motif_sorted_rank_8578singleton_motifs_{{REGION}}_subset_3866motifs_isoforms.feather",
        "motif_map": f"{BASE}/encode_rbp_motif_link_clustered_motifs_only.csv",
    },
    "archetype": {
        "rank_pattern": f"{ISO_BASE}/motif_sorted_rank_8578clustered_motifs_{{REGION}}_subset_1050archetypemotifs_isoform.feather",
        "motif_map": f"{BASE}/encode_rbp_archetype_motif_links_only.csv",
    },
    "FIMO": {
        "rank_pattern": f"{ISO_BASE}/motif_sorted_rank_3866clustered_motifs_FIMO_{{REGION}}.feather",
        "motif_map": f"{BASE}/encode_rbp_motif_link_clustered_motifs_only.csv",
    },
    "HOMER2": {
        "rank_pattern": f"{ISO_BASE}/motif_sorted_rank_3866clustered_motifs_HOMER_{{REGION}}.feather",
        "motif_map": f"{BASE}/encode_rbp_motif_link_clustered_motifs_only.csv",
    },
}


def evaluate_motif_predictions(rank_matrix, motif_rbp_map, eclip_df,
                                max_k, ratio=1, N=1):
    """Compute F1/Precision/Recall for each motif at a given top-k cutoff."""
    rbp_to_targets = eclip_df.groupby("RBP")["Gene"].apply(set).to_dict()
    motif_to_rbp = motif_rbp_map.set_index("Motif")["RBP"].to_dict()

    motif_rank_dict = {
        motif: {gene: rank for rank, gene in enumerate(row.dropna().values)}
        for motif, row in rank_matrix.iterrows()
    }

    results = []
    for motif, gene_rank in tqdm(motif_rank_dict.items()):
        if motif not in motif_to_rbp:
            continue
        rbp = motif_to_rbp[motif]
        if rbp not in rbp_to_targets:
            continue

        gold_targets = rbp_to_targets[rbp]
        motif_genes = set(gene_rank.keys())
        positives = list(motif_genes & gold_targets)

        if len(positives) < 10:
            continue
        negatives_pool = list(motif_genes - set(positives))
        if len(negatives_pool) < ratio * len(positives):
            continue

        f1_list, prec_list, rec_list = [], [], []
        for _ in range(N):
            negatives = np.random.choice(negatives_pool, size=ratio * len(positives), replace=False)
            test_genes = positives + list(negatives)
            test_ranks = sorted([(g, gene_rank[g]) for g in test_genes], key=lambda x: x[1])

            predicted_pos = set([g for g, _ in test_ranks[:max_k]])
            y_true = [1 if g in positives else 0 for g in test_genes]
            y_pred = [1 if g in predicted_pos else 0 for g in test_genes]

            prec_list.append(precision_score(y_true, y_pred, zero_division=0))
            rec_list.append(recall_score(y_true, y_pred, zero_division=0))
            f1_list.append(f1_score(y_true, y_pred, zero_division=0))

        results.append({
            "Motif": motif, "RBP": rbp,
            "Positive_Genes": len(positives),
            "F1": np.mean(f1_list),
            "Precision": np.mean(prec_list),
            "Recall": np.mean(rec_list),
        })

    return pd.DataFrame(results)


# ============================================================
# Main loop: 5 methods × 4 regions × 2 cell lines = 40 runs
# ============================================================
if __name__ == "__main__":
    for method, config in METHODS.items():
        motif_rbp_map = pd.read_csv(config["motif_map"])

        for region in REGIONS:
            rank_path = config["rank_pattern"].replace("{REGION}", region)

            if not os.path.exists(rank_path):
                print(f"[SKIP] {rank_path} not found")
                continue

            print(f"\nLoading rank matrix: {method} / {region}")
            rank_matrix = pd.read_feather(rank_path)

            for cell_line, eclip_path in CELL_LINES.items():
                print(f"  Evaluating: {method} | {region} | {cell_line}")
                eclip_df = pd.read_csv(eclip_path)

                result_df = evaluate_motif_predictions(
                    rank_matrix, motif_rbp_map, eclip_df,
                    max_k=MAX_K, ratio=RATIO, N=N_REPEATS
                )

                # Save results
                out_csv = f"motif_eval_boxplot_{cell_line}_top{MAX_K}_ratio{RATIO}_repeat{N_REPEATS}times_{method}_{region}_isoforms.csv"
                result_df.to_csv(out_csv, index=False)

                # Keep best motif per RBP and plot
                best = result_df.sort_values("F1", ascending=False).groupby("RBP").head(1)
                plot_df = best.melt(
                    id_vars=["Motif", "RBP"],
                    value_vars=["F1", "Precision", "Recall"],
                    var_name="Metric", value_name="Score"
                )

                plt.figure(figsize=(8, 6))
                sns.boxplot(x="Metric", y="Score", data=plot_df, palette="Set2")
                sns.stripplot(x="Metric", y="Score", data=plot_df, color="black", size=2, jitter=True, alpha=0.3)
                plt.title(f"Motif Prediction (Top-{MAX_K}, {method}, {region}, {cell_line})")
                plt.tight_layout()
                plt.savefig(out_csv.replace(".csv", ".pdf"))
                plt.close()

                print(f"    Saved: {out_csv}")

    print("\nAll evaluations complete!")
```

### Run

```bash
srun --mem=100G --time=20-00:00:00 --pty python run_F1_all_methods.py
```

> **Note:** The rank matrix is loaded once per method × region pair and reused for both cell lines, which saves significant I/O time compared to the original 40 separate scripts.

---

## Step 2 — Visualization: per-region comparison across methods

Compare F1/Precision/Recall across 5 methods within each mRNA region (merging K562 + HepG2 results).

### R script: `plot_F1_per_region.R`

```r
library(dplyr); library(ggplot2); library(readr); library(tidyr); library(RColorBrewer)

REGIONS <- c("3UTR", "5UTR", "CDS", "Introns")
method_colors <- brewer.pal(5, "Set2")

for (region in REGIONS) {
  files <- list.files(path = paste0("./", region, "/"), pattern = "*.csv", full.names = TRUE)

  # Batch read and annotate
  all_data <- lapply(files, function(file) {
    df <- read_csv(file)
    parts <- strsplit(basename(file), "_")[[1]]
    df$Method <- gsub(".csv", "", parts[8])
    df$CellLine <- parts[4]
    return(df)
  }) %>% bind_rows()

  # Retain best motif per RBP
  best_motifs <- all_data %>%
    group_by(RBP, CellLine, Method) %>%
    slice_max(order_by = F1, n = 1, with_ties = FALSE) %>%
    ungroup()

  # Convert to long format
  plot_data <- best_motifs %>%
    select(Method, CellLine, F1, Precision, Recall) %>%
    pivot_longer(cols = c("F1", "Precision", "Recall"),
                 names_to = "Metric", values_to = "Score")

  plot_data$Method <- factor(plot_data$Method,
                             levels = c("clustered", "singleton", "archetype", "FIMO", "HOMER2"))

  # Boxplot with jitter
  p <- ggplot(plot_data, aes(x = Method, y = Score, fill = Method)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.4, color = "black") +
    facet_wrap(~Metric, scales = "free_y", nrow = 3) +
    scale_fill_manual(values = method_colors) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 14, face = "bold")) +
    labs(title = paste("Method Comparison -", region),
         x = "Motif Scanning Method", y = "Score")

  pdf(paste0("motif_method_comparison_", region, "_isoforms.pdf"), width = 10, height = 8)
  print(p)
  dev.off()
}
```

---

## Step 3 — Visualization: best-region performance comparison

For each RBP, select the mRNA region where **Clustered motifs** achieve the highest F1, then compare all 5 methods using that same region. This shows the best-case performance of the Clustered method and whether other methods match it.

### R script: `plot_F1_best_region.R`

```r
library(dplyr); library(ggplot2); library(readr); library(tidyr)
library(stringr); library(RColorBrewer)

# Step 1: Batch read all CSV files across all regions
files <- list.files(path = "./z_3UTR_5UTR_CDS_Introns_results",
                    pattern = "*.csv", full.names = TRUE)

extract_info <- function(filename) {
  parts <- strsplit(basename(filename), "_")[[1]]
  list(method = gsub(".csv", "", parts[8]),
       cell_line = parts[4],
       region = parts[9])
}

all_data <- lapply(files, function(file) {
  df <- read_csv(file)
  info <- extract_info(file)
  df$Method <- info$method
  df$CellLine <- info$cell_line
  df$Region <- info$region
  return(df)
}) %>% bind_rows()

# Step 1.1: Per method x region x RBP, keep the best motif
best_per_method_region <- all_data %>%
  group_by(RBP, CellLine, Method, Region) %>%
  slice_max(order_by = F1, n = 1, with_ties = FALSE) %>%
  ungroup()

# Step 2: For each RBP + CellLine, find the region where Clustered achieves the best F1
best_clustered_regions <- best_per_method_region %>%
  filter(Method == "clustered") %>%
  group_by(RBP, CellLine) %>%
  slice_max(order_by = F1, n = 1, with_ties = FALSE) %>%
  select(RBP, CellLine, Region) %>%
  rename(BestRegion = Region)

# Step 3: Filter all methods to use Clustered's best region per RBP
best_data <- best_per_method_region %>%
  inner_join(best_clustered_regions, by = c("RBP", "CellLine")) %>%
  filter(Region == BestRegion)

# Step 4: Convert to long format
plot_data <- best_data %>%
  select(Method, CellLine, F1, Precision, Recall) %>%
  pivot_longer(cols = c("F1", "Precision", "Recall"),
               names_to = "Metric", values_to = "Score")

plot_data$Method <- factor(plot_data$Method,
                           levels = c("clustered", "singleton", "archetype", "FIMO", "HOMER2"))
method_colors <- brewer.pal(5, "Set2")

# Step 5: Plot
p <- ggplot(plot_data, aes(x = Method, y = Score, fill = Method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.4, color = "black") +
  facet_wrap(~Metric, scales = "free_y", nrow = 3) +
  scale_fill_manual(values = method_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14, face = "bold")) +
  labs(title = "Method Comparison (Each RBP Uses Clustered's Best Region)",
       x = "Motif Scanning Method", y = "Score")

pdf("motif_method_comparison_best_region_isoforms.pdf", width = 10, height = 8)
print(p)
dev.off()
```

---

## Pipeline Flowchart

```
Step 1.0: Top-k sensitivity test
          (Clustered × 3UTR × [500..5000])
          → Select optimal top-k = 2000
                    ↓
Step 1.1–1.5: F1 evaluation
          (5 methods × 4 regions × 2 cell lines = 40 runs)
          → Per-run CSV + boxplot PDF
                    ↓
Step 2: Per-region visualization
          (Merge K562 + HepG2 per region, compare 5 methods)
          → 4 comparison PDFs (one per region)
                    ↓
Step 3: Best-region visualization
          (Select Clustered's best region per RBP,
           compare all 5 methods at that region)
          → Final comparison PDF
```

## Key Improvement

The original pipeline contained **40+ nearly identical Python scripts** (~100 lines each, ~4000+ total lines) that differed only in 3 variables: `method`, `region`, and `cell_line`. This refactored version consolidates them into **a single parameterized script** with a triple nested loop, reducing the codebase from ~4000 lines to ~120 lines while producing identical outputs.
