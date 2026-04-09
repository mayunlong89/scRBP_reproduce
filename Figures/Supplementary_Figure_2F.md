

# ENCODE eCLIP Benchmarking Pipeline (Steps 7–12)

> **Date:** 2025-06-19 / 2025-06-20  
> **Purpose:** Benchmark scRBP motif enrichment results against ENCODE eCLIP data across 5 motif scanning methods and 4 mRNA regions.

---

## Table of Contents

1. [Overview](#overview)  
2. [Prerequisites](#prerequisites)  
3. [Step 7 — Convert pruned results to motif-by-RBP matrix](#step-7--convert-pruned-results-to-motif-by-rbp-matrix)  
4. [Step 8 — Add RBP names to the matrix (per cell line)](#step-8--add-rbp-names-to-the-matrix-per-cell-line)  
5. [Step 9 — Merge two cell lines into a single matrix](#step-9--merge-two-cell-lines-into-a-single-matrix)  
6. [Step 10 — Z-score normalization (row-wise, by motif)](#step-10--z-score-normalization-row-wise-by-motif)  
7. [Step 11 — Compute AUC and plot recovery curves (per method × region)](#step-11--compute-auc-and-plot-recovery-curves-per-method--region)  
8. [Step 12 — Compare recovery curves across 5 methods](#step-12--compare-recovery-curves-across-5-methods)  
9. [Step 13 — Merge all 4 mRNA regions and run final comparison](#step-13--merge-all-4-mrna-regions-and-run-final-comparison)

---

## Overview

This pipeline evaluates **5 motif scanning methods** across **4 mRNA regions** using ENCODE eCLIP data from **2 cell lines** (K562 and HepG2).

**5 Motif methods:**

| # | Method | Description |
|---|--------|-------------|
| 1 | Clustered | Clustered motifs |
| 2 | Singleton | Singleton motifs |
| 3 | Archetype | Archetype motifs |
| 4 | FIMO | FIMO-scanned motifs |
| 5 | HOMER2 | HOMER2-scanned motifs |

**4 mRNA regions:** `3UTR`, `5UTR`, `CDS`, `Introns`

**2 Cell lines:** `K562` (120 RBPs), `HepG2` (103 RBPs)

---

## Prerequisites

### Directory structure

```
# Tool scripts directory
tool=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE

# Base directory for isoform-level analysis
base=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE_isoforms/00_scRBP_isoforms
```

### Required Python scripts

| Script | Used in |
|--------|---------|
| `convert_prune_to_matrix_final.py` | [Step 7](https://github.com/mayunlong89/scRBP_reproduce/blob/main/Analyses/Isoform_motif_ranking/convert_prune_to_matrix_final.py) |
| `adding_RBP_names_with_HepG2.py` | [Step 8](https://github.com/mayunlong89/scRBP_reproduce/blob/main/Analyses/Isoform_motif_ranking/adding_RBP_names_with_HepG2.py) |
| `adding_RBP_names_with_K562.py` | [Step 8](https://github.com/mayunlong89/scRBP_reproduce/blob/main/Analyses/Isoform_motif_ranking/adding_RBP_names_with_K562.py) |
| `merge_motif_rbp_matrix_twocell_lines.py` | [Step 9](https://github.com/mayunlong89/scRBP_reproduce/blob/main/Analyses/Isoform_motif_ranking/merge_motif_rbp_matrix_twocell_lines.py) |
| `zscore_motif_by_rbp_matrix_Motifrow_scale.py` | [Step 10](https://github.com/mayunlong89/scRBP_reproduce/blob/main/Analyses/Isoform_motif_ranking/zscore_motif_by_rbp_matrix_Motifrow_scale.py) |
| `python_recovery_RBP_plot_final_main.py` | [Step 11](https://github.com/mayunlong89/scRBP_reproduce/blob/main/Analyses/Isoform_motif_ranking/python_recovery_RBP_plot_final_main.py) |
| `compare_recovery_curves.py` | [Step 12](https://github.com/mayunlong89/scRBP_reproduce/blob/main/Analyses/Isoform_motif_ranking/compare_recovery_curves.py) |
| `merge_rbp_best_rank.py` | [Step 13](https://github.com/mayunlong89/scRBP_reproduce/blob/main/Analyses/Isoform_motif_ranking/merge_rbp_best_rank.py) |
| `compare_recovery_curves_plotpdf.py` | [Step 13](https://github.com/mayunlong89/scRBP_reproduce/blob/main/Analyses/Isoform_motif_ranking/compare_recovery_curves_plotpdf.py) |

### Ground-truth files

| File | Used for |
|------|----------|
| `encode_rbp_motif_link_clustered_motifs_only_with_celltypes.csv` | Clustered / Singleton / FIMO / HOMER |
| `encode_rbp_archetype_motif_links_only_with_celltypes.csv` | Archetype |

---

## Step 7 — Convert pruned results to motif-by-RBP matrix

Convert the GSEA leading-edge results into a motif-by-RBP enrichment matrix.

```bash
tool=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE

srun --pty --mem=50G --time=10-00:00:00 \
  python $tool/convert_prune_to_matrix_final.py \
    -i global_nes_leading_edge_results.csv \
    -o motif_by_rbp_matrix_5mingenes.csv
```

> **Note:** Run this inside each method × region working directory. The input `global_nes_leading_edge_results.csv` is generated from earlier steps (Steps 1–6, not shown here).

---

## Step 8 — Add RBP names to the matrix (per cell line)

Annotate each motif-by-RBP matrix with proper RBP names for each cell line.

```bash
tool=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE

# HepG2
python $tool/adding_RBP_names_with_HepG2.py

# K562
python $tool/adding_RBP_names_with_K562.py
```

> **Note:** Run these inside each method × region working directory. Output files are named `motif_by_rbp_matrix_5mingenes_HepG2.csv` and `motif_by_rbp_matrix_5mingenes_K562.csv`.

---

## Step 9 — Merge two cell lines into a single matrix

Merge HepG2 and K562 matrices for each method × region combination.

### Usage pattern

```bash
tool=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE

python $tool/merge_motif_rbp_matrix_twocell_lines.py \
  -i1 <K562_DIR>/motif_by_rbp_matrix_5mingenes_K562.csv \
  -i2 <HepG2_DIR>/motif_by_rbp_matrix_5mingenes_HepG2.csv \
  -o  <output_DIR>/merged_motif_by_rbp_matrix_5mingenes_twocell_lines_<REGION>.csv
```

### Full commands (5 methods × 4 regions = 20 runs)

<details>
<summary><b>9.1 Clustered motifs (click to expand)</b></summary>

```bash
tool=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE
base=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE_isoforms/00_scRBP_isoforms
output_dir=$base/01_clustered

for region in 3UTR 5UTR CDS Introns; do
  k562_DIR=$base/z_ENCODE_K562_120RBPs_top0.05_5mingenes_adjustedMethods_none_clustered_motifs_${region}_isoforms
  hepg2_DIR=$base/z_ENCODE_HepG2_103RBPs_top0.05_5mingenes_adjustedMethods_none_clustered_motifs${region}_isoforms

  python $tool/merge_motif_rbp_matrix_twocell_lines.py \
    -i1 $k562_DIR/motif_by_rbp_matrix_5mingenes_K562.csv \
    -i2 $hepg2_DIR/motif_by_rbp_matrix_5mingenes_HepG2.csv \
    -o  $output_dir/merged_motif_by_rbp_matrix_5mingenes_twocell_lines_${region}.csv
done
```

</details>

<details>
<summary><b>9.2 Singleton motifs (click to expand)</b></summary>

```bash
output_dir=$base/02_singleton

for region in 3UTR 5UTR CDS Introns; do
  k562_DIR=$base/z_ENCODE_K562_120RBPs_top0.05_5mingenes_adjustedMethods_none_singleton_motifs_${region}_isoforms
  hepg2_DIR=$base/z_ENCODE_HepG2_103RBPs_top0.05_5mingenes_adjustedMethods_none_singleton_motifs${region}_isoforms

  python $tool/merge_motif_rbp_matrix_twocell_lines.py \
    -i1 $k562_DIR/motif_by_rbp_matrix_5mingenes_K562.csv \
    -i2 $hepg2_DIR/motif_by_rbp_matrix_5mingenes_HepG2.csv \
    -o  $output_dir/merged_motif_by_rbp_matrix_5mingenes_twocell_lines_${region}.csv
done
```

</details>

<details>
<summary><b>9.3 Archetype motifs (click to expand)</b></summary>

```bash
output_dir=$base/03_archetype

for region in 3UTR 5UTR CDS Introns; do
  k562_DIR=$base/z_ENCODE_K562_120RBPs_top0.05_5mingenes_adjustedMethods_none_archetype_motifs_${region}_isoforms
  hepg2_DIR=$base/z_ENCODE_HepG2_103RBPs_top0.05_5mingenes_adjustedMethods_none_archetype_motifs${region}_isoforms

  python $tool/merge_motif_rbp_matrix_twocell_lines.py \
    -i1 $k562_DIR/motif_by_rbp_matrix_5mingenes_K562.csv \
    -i2 $hepg2_DIR/motif_by_rbp_matrix_5mingenes_HepG2.csv \
    -o  $output_dir/merged_motif_by_rbp_matrix_5mingenes_twocell_lines_${region}.csv
done
```

</details>

<details>
<summary><b>9.4 FIMO motifs (click to expand)</b></summary>

```bash
output_dir=$base/04_FIMO

for region in 3UTR 5UTR CDS Introns; do
  k562_DIR=$base/z_ENCODE_K562_120RBPs_top0.05_5mingenes_adjustedMethods_none_FIMO_${region}_isoforms
  hepg2_DIR=$base/z_ENCODE_HepG2_103RBPs_top0.05_5mingenes_adjustedMethods_none_FIMO_${region}_isoforms

  python $tool/merge_motif_rbp_matrix_twocell_lines.py \
    -i1 $k562_DIR/motif_by_rbp_matrix_5mingenes_K562.csv \
    -i2 $hepg2_DIR/motif_by_rbp_matrix_5mingenes_HepG2.csv \
    -o  $output_dir/merged_motif_by_rbp_matrix_5mingenes_twocell_lines_${region}.csv
done
```

</details>

<details>
<summary><b>9.5 HOMER2 motifs (click to expand)</b></summary>

```bash
output_dir=$base/05_HOMER2

for region in 3UTR 5UTR CDS Introns; do
  k562_DIR=$base/z_ENCODE_K562_120RBPs_top0.05_5mingenes_adjustedMethods_none_HOMER_${region}_isoforms
  hepg2_DIR=$base/z_ENCODE_HepG2_103RBPs_top0.05_5mingenes_adjustedMethods_none_HOMER_${region}_isoforms

  python $tool/merge_motif_rbp_matrix_twocell_lines.py \
    -i1 $k562_DIR/motif_by_rbp_matrix_5mingenes_K562.csv \
    -i2 $hepg2_DIR/motif_by_rbp_matrix_5mingenes_HepG2.csv \
    -o  $output_dir/merged_motif_by_rbp_matrix_5mingenes_twocell_lines_${region}.csv
done
```

</details>

---

## Step 10 — Z-score normalization (row-wise, by motif)

Scale and normalize the merged motif-by-RBP matrix. Z-scores are computed row-wise (per motif across all RBPs).

### Usage pattern

```bash
python $tool/zscore_motif_by_rbp_matrix_Motifrow_scale.py \
  -i <input_merged_matrix>.csv \
  -o <output_normalized_matrix>_nes.csv
```

### Full commands (5 methods × 4 regions = 20 runs)

```bash
tool=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE
base=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE_isoforms/00_scRBP_isoforms

# Method directories
declare -A method_dirs=(
  ["clustered"]="01_clustered"
  ["singleton"]="02_singleton"
  ["archetype"]="03_archetype"
  ["FIMO"]="04_FIMO"
  ["HOMER2"]="05_HOMER2"
)

for method in clustered singleton archetype FIMO HOMER2; do
  output_DIR=$base/${method_dirs[$method]}
  for region in 3UTR 5UTR CDS Introns; do
    python $tool/zscore_motif_by_rbp_matrix_Motifrow_scale.py \
      -i $output_DIR/merged_motif_by_rbp_matrix_5mingenes_twocell_lines_${region}.csv \
      -o $output_DIR/merged_motif_by_rbp_matrix_5mingenes_twocell_lines_${region}_nes.csv
  done
done
```

---

## Step 11 — Compute AUC and plot recovery curves (per method × region)

For each method and mRNA region, compute the RBP recovery curve and AUC score.

### Usage pattern

```bash
python $tool/python_recovery_RBP_plot_final_main.py \
  -m <normalized_matrix>_nes.csv \
  -t <ground_truth>.csv \
  -r <output_rank>.csv \
  -f <output_figure>.png \
  -a <output_auc>.txt \
  --max_rank 100
```

### Ground-truth file mapping

| Method | Ground-truth file |
|--------|------------------|
| Clustered | `encode_rbp_motif_link_clustered_motifs_only_with_celltypes.csv` |
| Singleton | `encode_rbp_motif_link_clustered_motifs_only_with_celltypes.csv` |
| Archetype | `encode_rbp_archetype_motif_links_only_with_celltypes.csv` |
| FIMO | `encode_rbp_motif_link_clustered_motifs_only_with_celltypes.csv` |
| HOMER2 | `encode_rbp_motif_link_clustered_motifs_only_with_celltypes.csv` |

### Full commands (5 methods × 4 regions = 20 runs)

```bash
tool=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE
base=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE_isoforms/00_scRBP_isoforms

# Define method configs: work_dir, output_dir, truth_file
declare -A work_dirs=( ["clustered"]="01_clustered" ["singleton"]="02_singleton" ["archetype"]="03_archetype" ["FIMO"]="04_FIMO" ["HOMER2"]="05_HOMER2" )
declare -A out_dirs=(  ["clustered"]="1_clustered_output" ["singleton"]="2_singleton_output" ["archetype"]="3_archetype_output" ["FIMO"]="4_FIMO_output" ["HOMER2"]="5_HOMER2_output" )
declare -A truth_files=(
  ["clustered"]="encode_rbp_motif_link_clustered_motifs_only_with_celltypes.csv"
  ["singleton"]="encode_rbp_motif_link_clustered_motifs_only_with_celltypes.csv"
  ["archetype"]="encode_rbp_archetype_motif_links_only_with_celltypes.csv"
  ["FIMO"]="encode_rbp_motif_link_clustered_motifs_only_with_celltypes.csv"
  ["HOMER2"]="encode_rbp_motif_link_clustered_motifs_only_with_celltypes.csv"
)

for method in clustered singleton archetype FIMO HOMER2; do
  work_DIR=$base/${work_dirs[$method]}
  output_DIR=$base/${out_dirs[$method]}
  truth=$tool/${truth_files[$method]}

  for region in 3UTR 5UTR CDS Introns; do
    python $tool/python_recovery_RBP_plot_final_main.py \
      -m $work_DIR/merged_motif_by_rbp_matrix_5mingenes_twocell_lines_${region}_nes.csv \
      -t $truth \
      -r $output_DIR/motif_best_rbp_rank_twocell_lines_${region}_nes.csv \
      -f $output_DIR/motif_recovery_curve_top100_twocell_lines_${region}_nes.png \
      -a $output_DIR/motif_auc_score_twocell_lines_${region}_nes.txt \
      --max_rank 100
  done
done
```

---

## Step 12 — Compare recovery curves across 5 methods

Compare all 5 motif scanning methods side-by-side for each mRNA region.

### Usage pattern

```bash
python $tool/compare_recovery_curves.py \
  -i <method1_rank>.csv <method2_rank>.csv ... <method5_rank>.csv \
  -o <output_figure>.png \
  --max_rank 100
```

### Full commands (one per region)

```bash
tool=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE
base=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE_isoforms/00_scRBP_isoforms
output_dir=$base/00_merged_5methods_output

clustered_dir=$base/1_clustered_output
singleton_dir=$base/2_singleton_output
arch_dir=$base/3_archetype_output
fimo_dir=$base/4_FIMO_output
homer_dir=$base/5_HOMER2_output

for region in 3UTR 5UTR CDS Introns; do
  python $tool/compare_recovery_curves.py \
    -i $clustered_dir/motif_best_rbp_rank_twocell_lines_${region}_nes.csv \
       $arch_dir/motif_best_rbp_rank_twocell_lines_${region}_nes.csv \
       $singleton_dir/motif_best_rbp_rank_twocell_lines_${region}_nes.csv \
       $fimo_dir/motif_best_rbp_rank_twocell_lines_${region}_nes.csv \
       $homer_dir/motif_best_rbp_rank_twocell_lines_${region}_nes.csv \
    -o $output_dir/recovery_curve_all_5methods_${region}.png \
    --max_rank 100
done
```

---

## Step 13 — Merge all 4 mRNA regions and run final comparison

### 13.1 Merge 4 regions per method

For each method, merge the best RBP ranks from 3UTR, 5UTR, CDS, and Introns into a single file, then compute recovery curve and AUC.

```bash
tool=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE
base=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE_isoforms/00_scRBP_isoforms
output_dir=$base/00_merged_5methods_output

declare -A rank_dirs=(
  ["clusteredMotifs"]="1_clustered_output"
  ["SingletonMotifs"]="2_singleton_output"
  ["archetypeMotifs"]="3_archetype_output"
  ["FIMO_Motifs"]="4_FIMO_output"
  ["HOMER_Motifs"]="5_HOMER2_output"
)

for motif_tag in clusteredMotifs SingletonMotifs archetypeMotifs FIMO_Motifs HOMER_Motifs; do
  rank_dir=$base/${rank_dirs[$motif_tag]}

  python $tool/merge_rbp_best_rank.py \
    --input_files \
      $rank_dir/motif_best_rbp_rank_twocell_lines_3UTR_nes.csv \
      $rank_dir/motif_best_rbp_rank_twocell_lines_5UTR_nes.csv \
      $rank_dir/motif_best_rbp_rank_twocell_lines_CDS_nes.csv \
      $rank_dir/motif_best_rbp_rank_twocell_lines_Introns_nes.csv \
    --output     $output_dir/merged_best_rank_5utr_3utr_cds_intron_${motif_tag}.csv \
    --fig_output $output_dir/merged_recovery_curve_5utr_3utr_cds_intron_${motif_tag}.png \
    --auc_output $output_dir/AUC_score_5utr_3utr_cds_intron_${motif_tag}.txt \
    --max_rank 100
done
```

### 13.2 Final 5-method comparison (merged across all regions)

```bash
tool=/mnt/isilon/gandal_lab/mayl/.../05_benchmark_ENCODE
output_dir=$base/00_merged_5methods_output

# PNG output
python $tool/compare_recovery_curves.py \
  -i $output_dir/merged_best_rank_5utr_3utr_cds_intron_clusteredMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_archetypeMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_SingletonMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_FIMO_Motifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_HOMER_Motifs.csv \
  -o $output_dir/recovery_curve_all_5methods_merged_4regions.png \
  --max_rank 100

# PDF output (publication quality)
python $tool/compare_recovery_curves_plotpdf.py \
  -i $output_dir/merged_best_rank_5utr_3utr_cds_intron_clusteredMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_archetypeMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_SingletonMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_FIMO_Motifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_HOMER_Motifs.csv \
  -o $output_dir/recovery_curve_all_5methods_merged_4regions.pdf \
  --max_rank 100

# Top-20 zoomed view
python $tool/compare_recovery_curves.py \
  -i $output_dir/merged_best_rank_5utr_3utr_cds_intron_clusteredMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_archetypeMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_SingletonMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_FIMO_Motifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_HOMER_Motifs.csv \
  -o $output_dir/recovery_curve_all_5methods_merged_4regions_top20.png \
  --max_rank 20
```

---


## Pipeline Flowchart

```
Step 7: Pruned results → motif-by-RBP matrix
            ↓
Step 8: Add RBP names (per cell line: K562, HepG2)
            ↓
Step 9: Merge two cell lines → single matrix
            ↓                (× 5 methods × 4 regions)
Step 10: Z-score normalization (row-wise)
            ↓
Step 11: Recovery curve + AUC (per method × region)
            ↓
Step 12: Compare 5 methods (per region)
            ↓
Step 13: Merge 4 regions → final 5-method comparison
```
