# Isoform-Centric Stability Assessment: Jaccard Similarity of 73 Common RBPs Between K562 and HepG2

> **Date:** 2025-06-20 / 2025-06-21  
> **Author:** Yunlong Ma

---

## Table of Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Step 1 — Extract real and background motif-RBP links](#step-1--extract-real-and-background-motif-rbp-links)
4. [Step 2 — Calculate scaled Jaccard similarity](#step-2--calculate-scaled-jaccard-similarity)
5. [Step 3 — Bootstrap significance testing](#step-3--bootstrap-significance-testing)

---

## Overview

### Motivation

Previous analyses (gene-level benchmarking, Steps 7–13) established that **clustered motifs** achieve the best **accuracy** among the 5 ranking strategies. This pipeline extends the evaluation to **stability** — i.e., how reproducible are the results across independent biological contexts?

### Approach

The 223 gold-standard eCLIP datasets (150 RBPs) come from two cell lines: **K562** (120 RBPs) and **HepG2** (103 RBPs), with **73 common RBPs** shared between them. For each motif scanning method, we:

1. **Extract** the real (eCLIP-validated) and background motif-RBP links from each cell line's scRBP results.
2. **Compute** the Jaccard similarity of target gene sets between K562 and HepG2 for the 73 common RBPs, using z-score transformation to correct for baseline differences between the two cell lines.
3. **Test significance** via bootstrap, comparing real Jaccard scores against background null distributions.

By comparing Jaccard similarity across the 5 methods, we assess **stability** (cross-cell-line reproducibility).

### Design matrix

| Dimension | Values |
|-----------|--------|
| **Methods** | FIMO, HOMER, Clustered, Singleton (individual), Archetype |
| **Regions** | 3UTR, 5UTR, CDS, Introns |
| **Cell lines** | K562, HepG2 |
| **Link types** | Real (eCLIP-validated), Background |

Total runs per step: **5 methods × 4 regions × 2 cell lines × 2 link types = 80** (Step 1) or **5 × 4 × 2 = 40** (Steps 2–3).

---

## Prerequisites

### Required scripts

| Script | Step | Language |
|--------|------|----------|
| `filter_real_background_links.R` | 1 | R |
| `calc_scaled_jaccard.py` | 2 | Python |
| `bootstrap_cal.R` | 3 | R |

### Directory variables

```bash
# Tool scripts
TOOL=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2

# eCLIP ground-truth references
Real_link_DIR=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE
Reference_DIR=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE_isoforms

# Base for scRBP isoform results
base=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE_isoforms/00_scRBP_isoforms

# Output base
jaccard_base=$base/00_Jaccard_similarity_isoform_between_twoCellLines
```

### Ground-truth files

| Method | Real-link file |
|--------|---------------|
| FIMO, HOMER, Clustered, Singleton | `encode_rbp_motif_link_clustered_motifs_only.csv` |
| Archetype | `encode_rbp_archetype_motif_links_only.csv` |

### eCLIP reference files (for background scaling)

| Cell line | File |
|-----------|------|
| K562 | `ENCODE_K562_120RBPs_targets_isoforms.csv` |
| HepG2 | `ENCODE_HepG2_103RBPs_targets_isoforms.csv` |

---

## Step 1 — Extract real and background motif-RBP links

For each method × region × cell line, extract both eCLIP-validated (real) and background links from the scRBP enrichment results.

### Usage

```bash
Rscript $TOOL/filter_real_background_links.R \
  <scRBP_results>/global_nes_leading_edge_results.csv \
  <ground_truth_links>.csv \
  <output_real_links>.csv \
  <output_background_links>.csv
```

### Automated loop

```bash
TOOL=/mnt/isilon/gandal_lab/mayl/.../z_jaccard_similarity_nes_K562_HepG2
Real_link_DIR=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE
base=/mnt/isilon/gandal_lab/mayl/.../00_scRBP_isoforms
jaccard_base=$base/00_Jaccard_similarity_isoform_between_twoCellLines

# ---- Method configurations ----
# Format: method_tag | K562_dir_pattern | HepG2_dir_pattern | truth_file
# Note: K562 uses "_<REGION>_" pattern; HepG2 uses "<REGION>_" (no leading underscore for some)

declare -A truth_files=(
  ["FIMO"]="encode_rbp_motif_link_clustered_motifs_only.csv"
  ["HOMER"]="encode_rbp_motif_link_clustered_motifs_only.csv"
  ["Clustered"]="encode_rbp_motif_link_clustered_motifs_only.csv"
  ["individual"]="encode_rbp_motif_link_clustered_motifs_only.csv"
  ["archetype"]="encode_rbp_archetype_motif_links_only.csv"
)

# K562 directory name patterns (method portion)
declare -A k562_patterns=(
  ["FIMO"]="FIMO"
  ["HOMER"]="HOMER"
  ["Clustered"]="clustered_motifs"
  ["individual"]="singleton_motifs"
  ["archetype"]="archetype_motifs"
)

# HepG2 directory name patterns (method portion)
declare -A hepg2_patterns=(
  ["FIMO"]="FIMO"
  ["HOMER"]="HOMER"
  ["Clustered"]="clustered_motifs"
  ["individual"]="singleton_motifs"
  ["archetype"]="archetype_motifs"
)

REGIONS=("3UTR" "5UTR" "CDS" "Introns")

for method in FIMO HOMER Clustered individual archetype; do
  truth=$Real_link_DIR/${truth_files[$method]}

  for region in "${REGIONS[@]}"; do
    output_DIR=$jaccard_base/$region
    mkdir -p $output_DIR

    # K562 scRBP results directory
    k562_DIR=$base/z_ENCODE_K562_120RBPs_top0.05_5mingenes_adjustedMethods_none_${k562_patterns[$method]}_${region}_isoforms

    # HepG2 scRBP results directory
    hepg2_DIR=$base/z_ENCODE_HepG2_103RBPs_top0.05_5mingenes_adjustedMethods_none_${hepg2_patterns[$method]}${region}_isoforms

    # K562
    Rscript $TOOL/filter_real_background_links.R \
      $k562_DIR/global_nes_leading_edge_results.csv \
      $truth \
      $output_DIR/${method}_K562_filtered_real_links.csv \
      $output_DIR/${method}_K562_filtered_background_links.csv

    # HepG2
    Rscript $TOOL/filter_real_background_links.R \
      $hepg2_DIR/global_nes_leading_edge_results.csv \
      $truth \
      $output_DIR/${method}_HepG2_filtered_real_links.csv \
      $output_DIR/${method}_HepG2_filtered_background_links.csv
  done
done
```

> **Note on directory naming:** The K562 directories use `_<REGION>_isoforms` (with leading underscore before region), while some HepG2 directories use `<REGION>_isoforms` (no leading underscore). Verify paths match your filesystem before running.

---

## Step 2 — Calculate scaled Jaccard similarity

For each method × region, compute the scaled Jaccard similarity between K562 and HepG2 for both real and background links. The scaling uses eCLIP reference data to correct for baseline differences.

### Usage

```bash
python $TOOL/calc_scaled_jaccard.py \
  --k562_real  <K562_filtered_links>.csv \
  --hepg2_real <HepG2_filtered_links>.csv \
  --k562_bg    <K562_eCLIP_reference>.csv \
  --hepg2_bg   <HepG2_eCLIP_reference>.csv \
  --output     <scaled_jaccard_output>.csv
```

### Automated loop

```bash
TOOL=/mnt/isilon/gandal_lab/mayl/.../z_jaccard_similarity_nes_K562_HepG2
Reference_DIR=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE_isoforms
jaccard_base=$base/00_Jaccard_similarity_isoform_between_twoCellLines

METHODS=("FIMO" "HOMER" "Clustered" "individual" "archetype")
REGIONS=("3UTR" "5UTR" "CDS" "Introns")

for method in "${METHODS[@]}"; do
  for region in "${REGIONS[@]}"; do
    output_DIR=$jaccard_base/$region

    # Real links
    python $TOOL/calc_scaled_jaccard.py \
      --k562_real  $output_DIR/${method}_K562_filtered_real_links.csv \
      --hepg2_real $output_DIR/${method}_HepG2_filtered_real_links.csv \
      --k562_bg    $Reference_DIR/ENCODE_K562_120RBPs_targets_isoforms.csv \
      --hepg2_bg   $Reference_DIR/ENCODE_HepG2_103RBPs_targets_isoforms.csv \
      --output     $output_DIR/${method}_K562_HepG2_scaled_jaccard_similarity.csv

    # Background links
    python $TOOL/calc_scaled_jaccard.py \
      --k562_real  $output_DIR/${method}_K562_filtered_background_links.csv \
      --hepg2_real $output_DIR/${method}_HepG2_filtered_background_links.csv \
      --k562_bg    $Reference_DIR/ENCODE_K562_120RBPs_targets_isoforms.csv \
      --hepg2_bg   $Reference_DIR/ENCODE_HepG2_103RBPs_targets_isoforms.csv \
      --output     $output_DIR/${method}_K562_HepG2_scaled_jaccard_similarity_background.csv
  done
done
```

---

## Step 3 — Bootstrap significance testing

For each method × region, test whether the real scaled Jaccard similarity is significantly higher than the background null distribution using bootstrap resampling.

### Usage

```bash
Rscript $TOOL/bootstrap_cal.R \
  <real_jaccard>.csv \
  <background_jaccard>.csv \
  Z_Scaled_Jaccard \
  <output_plot>.pdf \
  <output_pvalue>.txt
```

### Automated loop

```bash
TOOL=/mnt/isilon/gandal_lab/mayl/.../z_jaccard_similarity_nes_K562_HepG2
jaccard_base=$base/00_Jaccard_similarity_isoform_between_twoCellLines

METHODS=("FIMO" "HOMER" "Clustered" "individual" "archetype")
REGIONS=("3UTR" "5UTR" "CDS" "Introns")

for method in "${METHODS[@]}"; do
  for region in "${REGIONS[@]}"; do
    output_DIR=$jaccard_base/$region

    Rscript $TOOL/bootstrap_cal.R \
      $output_DIR/${method}_K562_HepG2_scaled_jaccard_similarity.csv \
      $output_DIR/${method}_K562_HepG2_scaled_jaccard_similarity_background.csv \
      Z_Scaled_Jaccard \
      $output_DIR/${method}_K562_HepG2_bootstrap_null_distribution_scaled_isoform.pdf \
      $output_DIR/${method}_K562_HepG2_pval_scaled.txt
  done
done
```

---

## Output Structure

After running all steps, each region directory will contain:

```
jaccard_base/
├── 3UTR/
│   ├── {method}_K562_filtered_real_links.csv          # Step 1
│   ├── {method}_K562_filtered_background_links.csv    # Step 1
│   ├── {method}_HepG2_filtered_real_links.csv         # Step 1
│   ├── {method}_HepG2_filtered_background_links.csv   # Step 1
│   ├── {method}_K562_HepG2_scaled_jaccard_similarity.csv            # Step 2
│   ├── {method}_K562_HepG2_scaled_jaccard_similarity_background.csv # Step 2
│   ├── {method}_K562_HepG2_bootstrap_null_distribution_scaled_isoform.pdf  # Step 3
│   └── {method}_K562_HepG2_pval_scaled.txt            # Step 3
├── 5UTR/
│   └── ... (same structure)
├── CDS/
│   └── ... (same structure)
└── Introns/
    └── ... (same structure)
```

Where `{method}` ∈ {`FIMO`, `HOMER`, `Clustered`, `individual`, `archetype`}.

---

## Pipeline Flowchart

```
Step 1: Extract real & background links
        (per method × region × cell line)
                    ↓
Step 2: Compute scaled Jaccard similarity
        (K562 vs HepG2, corrected by eCLIP background)
        → Real Jaccard + Background Jaccard
                    ↓
Step 3: Bootstrap significance test
        (Real vs Background null distribution)
        → p-value + null distribution plot
```

**Total runs:**
- Step 1: 5 methods × 4 regions × 2 cell lines = **40 calls**
- Step 2: 5 methods × 4 regions × 2 (real + background) = **40 calls**
- Step 3: 5 methods × 4 regions = **20 calls**
