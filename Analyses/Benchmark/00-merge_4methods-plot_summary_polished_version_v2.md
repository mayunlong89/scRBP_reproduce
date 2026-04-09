# Benchmarking Five Motif Scanning Methods Against 223 ENCODE eCLIP Datasets

## Overview

This workflow benchmarks five RBP motif scanning strategies by evaluating their ability to recover known RBP-target interactions from ENCODE eCLIP data. Recovery curves and AUC scores are computed at both the **gene level** and **isoform level** across four mRNA domains.

### Methods Compared

| Method | Description |
|--------|-------------|
| **Clustered motifs** | Clustered motifs from scRBP |
| **Archetype motifs** | Archetype motifs from scRBP |
| **Singleton motifs** | Individual motifs run independently from scRBP |
| **FIMO** | Motif scanning via FIMO (MEME Suite) |
| **HOMER2** | Motif scanning via HOMER2 |

### Genomic Regions

All analyses are stratified by four mRNA domains: **3'UTR**, **5'UTR**, **CDS**, and **Introns**.

---

## Workflow

```
┌──────────────────────────────────────────────────────────────────────┐
│                     INPUT: 223 ENCODE eCLIP datasets                │
│              (RBP binding sites as ground truth)                     │
└──────────────────────────┬───────────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────────────┐
│          STEP 1: Per-Region Recovery Curve Comparison                │
│                                                                      │
│  For each mRNA domain (3'UTR, 5'UTR, CDS, Introns):                │
│    • Compute NES-based best RBP rank for each method                │
│    • Generate per-region recovery curves comparing all 5 methods    │
│                                                                      │
│  Script: compare_recovery_curves.py                                  │
│  Input:  motif_best_rbp_rank_twocell_lines_{region}_nes.csv         │
│          (one per method × region)                                   │
│  Output: recovery_curve_all_5methods_{region}.png                    │
└──────────────────────────┬───────────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────────────┐
│        STEP 2: Merge Four Genomic Regions Per Method                 │
│                                                                      │
│  For each of the 5 methods independently:                            │
│    • Merge best RBP ranks across 3'UTR, 5'UTR, CDS, and Introns   │
│    • Generate individual recovery curves and AUC scores              │
│                                                                      │
│  Script: merge_rbp_best_rank.py                                      │
│  Input:  motif_best_rbp_rank_twocell_lines_{region}_nes.csv         │
│          (4 regions per method)                                      │
│  Output: merged_best_rank_5utr_3utr_cds_intron_{method}.csv         │
│          AUC_score_5utr_3utr_cds_intron_{method}.txt                │
└──────────────────────────┬───────────────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────────────┐
│     STEP 3: Cross-Method Comparison (Merged Across Regions)          │
│                                                                      │
│  Compare all 5 methods using merged (4-region) rankings:             │
│    • Generate combined recovery curves (top 20 / 50 / 100 ranks)    │
│    • Export publication-quality PDF figures                           │
│                                                                      │
│  Script: compare_recovery_curves.py                                  │
│          compare_recovery_curves_plotpdf.py  (PDF output)            │
│  Input:  merged_best_rank_5utr_3utr_cds_intron_{method}.csv         │
│          (one per method)                                            │
│  Output: recovery_curve_all_5methods_merged4domains_top{N}.pdf       │
└──────────────────────────────────────────────────────────────────────┘
```

---

## I. Gene-Level AUC Analysis

### Step 1: Per-Region Recovery Curves

Compare all five methods within each genomic region independently.

```bash
# --- Configuration ---
tool=05_benchmark_ENCODE
arch_dir=z_archetype_motifs/output
clustered_dir=z_clustered_motifs/output
singleton_dir=z_individual_run_clusteredmotifs/output
fimo_output_DIR=02_fimo/00_merged_2batches/output
HOMER_output_DIR=01_HOMER/00_HOMER_merge_batches/output

# Run for each region: 3UTR, 5UTR, CDS, Introns
for region in 3UTR 5UTR CDS Introns; do
  python $tool/compare_recovery_curves.py \
    -i $clustered_dir/motif_best_rbp_rank_twocell_lines_${region}_nes.csv \
       $arch_dir/motif_best_rbp_rank_twocell_lines_${region}_nes.csv \
       $singleton_dir/motif_best_rbp_rank_twocell_lines_${region}_nes.csv \
       $fimo_output_DIR/motif_best_rbp_rank_twocell_lines_${region}_nes.csv \
       $HOMER_output_DIR/motif_best_rbp_rank_twocell_lines_${region}_nes.csv \
    -o $output_dir/recovery_curve_all_5methods_${region}.png \
    --max_rank 100
done
```

### Step 2: Merge Four Genomic Regions Per Method

Aggregate best RBP rankings across all four mRNA domains for each method.

```bash
output_dir=z_merged_four_features_compare

# Run for each method: Clustered, Singleton, Archetype, FIMO, HOMER
# Example for Clustered motifs:
python $tool/merge_rbp_best_rank.py \
    --input_files \
      $clustered_dir/motif_best_rbp_rank_twocell_lines_3UTR_nes.csv \
      $clustered_dir/motif_best_rbp_rank_twocell_lines_5UTR_nes.csv \
      $clustered_dir/motif_best_rbp_rank_twocell_lines_CDS_nes.csv \
      $clustered_dir/motif_best_rbp_rank_twocell_lines_Introns_nes.csv \
    --output $output_dir/merged_best_rank_5utr_3utr_cds_intron_clusteredMotifs.csv \
    --fig_output $output_dir/merged_recovery_curve_clusteredMotifs.png \
    --auc_output $output_dir/AUC_score_clusteredMotifs.txt \
    --max_rank 100

# Repeat analogously for: Singleton, Archetype, FIMO, HOMER
```

### Step 3: Cross-Method Comparison

Compare all five methods using merged region-aggregated rankings.

```bash
# Generate recovery curves at different rank cutoffs (top 20, 50, 100)
for max_rank in 20 50 100; do
  python $tool/compare_recovery_curves.py \
    -i $output_dir/merged_best_rank_5utr_3utr_cds_intron_clusteredMotifs.csv \
       $output_dir/merged_best_rank_5utr_3utr_cds_intron_archetypeMotifs.csv \
       $output_dir/merged_best_rank_5utr_3utr_cds_intron_SingletonMotifs.csv \
       $output_dir/merged_best_rank_5utr_3utr_cds_intron_FIMO_Motifs.csv \
       $output_dir/merged_best_rank_5utr_3utr_cds_intron_HOMER_Motifs.csv \
    -o $output_dir/recovery_curve_all_5methods_merged4domains_top${max_rank}.png \
    --max_rank $max_rank
done

# Export publication-quality PDF
python $tool/compare_recovery_curves_plotpdf.py \
  -i $output_dir/merged_best_rank_5utr_3utr_cds_intron_clusteredMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_archetypeMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_SingletonMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_FIMO_Motifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_HOMER_Motifs.csv \
  -o $output_dir/recovery_curve_all_5methods_merged4domains_top100.pdf \
  --max_rank 100
```

---

## II. Isoform-Level AUC Analysis

The same three-step workflow is applied at the **isoform level** to evaluate method performance on transcript-specific RBP binding.

### Step 1: Merge Four Genomic Regions Per Method

```bash
output_dir=05_benchmark_ENCODE_isoforms/00_scRBP_isoforms/00_merged_5methods_output

# Run merge_rbp_best_rank.py for each of the 5 methods
# (analogous to gene-level Step 2, using isoform-level input directories)
```

### Step 2: Cross-Method Comparison

```bash
# Compare all 5 methods at top 20 and top 100 rank cutoffs
python $tool/compare_recovery_curves.py \
  -i $output_dir/merged_best_rank_5utr_3utr_cds_intron_clusteredMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_archetypeMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_SingletonMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_FIMO_Motifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_HOMER_Motifs.csv \
  -o $output_dir/recovery_curve_all_5methods_merged4domains.png \
  --max_rank 100

# Export publication-quality PDF
python $tool/compare_recovery_curves_plotpdf.py \
  -i $output_dir/merged_best_rank_5utr_3utr_cds_intron_clusteredMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_archetypeMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_SingletonMotifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_FIMO_Motifs.csv \
     $output_dir/merged_best_rank_5utr_3utr_cds_intron_HOMER_Motifs.csv \
  -o $output_dir/recovery_curve_all_5methods_merged4domains_isoform.pdf \
  --max_rank 100
```

---

## Key Scripts

| Script | Purpose |
|--------|---------|
| `compare_recovery_curves.py` | Compare recovery curves across methods (PNG output) |
| `compare_recovery_curves_plotpdf.py` | Generate publication-quality PDF recovery curves |
| `merge_rbp_best_rank.py` | Merge best RBP ranks across four mRNA domains and compute AUC |
