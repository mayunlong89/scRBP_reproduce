# RBP Binding-Region Preference Across 30 Tissues

## Overview

This pipeline characterizes **RBP binding-region preferences** across 2 brain tissues and 28 human tissues from the Tabula Sapiens single-cell atlas. For each tissue, we run the full scRBP pipeline (infer → merge → modules → prune) using 30 independent GRNBoost2 seeds, then analyze which mRNA regions (3'UTR, 5'UTR, CDS, Introns) each RBP preferentially binds across tissues.

**Experimental design:**
- 2 brain tissues and 28 other human tissues, each with ~20K cells (geometric-sketched from Tabula Sapiens)
- 30 GRNBoost2 runs per tissue (seeds 1001–1030) for consensus network construction
- 4 mRNA regions per tissue: 3'UTR, 5'UTR, CDS, Introns
- Regulon filtering: min 10 genes per regulon

**Directory structure:**
```
28organs_single_cell_atlas/
├── loom_dir/                              # Input .loom files per tissue
├── z_GRNBoost2_{TISSUE}_30times/          # Per-tissue output directory
│   ├── z_GRNBoost2_{TISSUE}_TSP1_seed*_correlation_v7.tsv   # Step 1: 30 GRN runs
│   ├── merged_scRBP_GRNs_consensus.tsv                       # Step 2: merged consensus
│   ├── z_GRNBoost2_{TISSUE}_20Kcells_merged_*_modules.tsv    # Step 3: modules
│   ├── Results_final_*_top1500_{REGION}/                      # Step 4: pruned regulons
│   └── scRBP_allRegions_min10genes.symbol_default.gmt         # Step 6: merged GMT
```

---

## Step 1: GRN Inference (30 seeds per tissue)

Run `scRBP_infer_v7_cor.py` 30 times per tissue using a SLURM array job (seeds 1001–1030). Each seed produces an independent GRNBoost2 network with correlation information.

```bash
# Submit for each tissue
sbatch z_grn_array_GRN_for_28organs.sh Bladder
sbatch z_grn_array_GRN_for_28organs.sh Heart
# ... (repeat for all 28 tissues)
```

### z_grn_array_GRN_for_28organs.sh

```bash
#!/usr/bin/env bash
#SBATCH -J GRN_${ORGAN}
#SBATCH --mem=150G
#SBATCH --cpus-per-task=20
#SBATCH --time=10-00:00:00
#SBATCH --array=1001-1030
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail

ORGAN="$1"
SEED="${SLURM_ARRAY_TASK_ID}"

TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/Palmer_hardwick_short_read
LOOM_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/28organs_single_cell_atlas/loom_dir
OUT_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/28organs_single_cell_atlas/z_GRNBoost2_${ORGAN}_30times

srun python ${TOOL}/scRBP_infer_v7_cor.py \
  --matrix    ${LOOM_DIR}/${ORGAN}_TSP1_30_v2.loom \
  --rbp_list  ${TOOL}/rbp_list_616RBPs.tsv \
  --output    ${OUT_DIR}/z_GRNBoost2_${ORGAN}_TSP1_seed${SEED}_correlation_v7.tsv \
  --n_workers 10 \
  --batch_size 5 \
  --threshold 0.03 \
  --method    grnboost2 \
  --correlation True \
  --seed      ${SEED} \
  --log       ${OUT_DIR}/z_GRNBoost2_${ORGAN}_TSP1_seed${SEED}_correlation_v7.log
```

---

## Step 2: Merge 30 GRN Runs into Consensus Network

Aggregate the 30 independent GRNBoost2 outputs into a single consensus network. An edge is retained if it appears in at least 10 of 30 runs (≥30% present rate).

```bash
sbatch z_grn_merged_30times_for_28organs.sh Bladder
```

### z_grn_merged_30times_for_28organs.sh

```bash
#!/usr/bin/env bash
#SBATCH -J MERGE_GRN_${ORGAN}
#SBATCH --mem=200G
#SBATCH --cpus-per-task=1
#SBATCH --time=10-00:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail

ORGAN="$1"

# Paths
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data
DATA=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/28organs_single_cell_atlas/z_GRNBoost2_${ORGAN}_30times

SCRIPT="${TOOL}/merge_grnboost2_seeds_v2.py"
PATTERN="${DATA}/z_GRNBoost2_${ORGAN}_TSP1_seed*_correlation_v7_scRBP_GRNs.tsv"
OUTFILE="${DATA}/merged_scRBP_GRNs_consensus.tsv"

# Consensus parameters
CORR_THRESHOLD="0.0"
N_PRESENT="10"
PRESENT_RATE="0.3"

# Validation
mkdir -p logs
[[ -d "${DATA}" ]]   || { echo "[ERROR] Missing DATA dir: ${DATA}" >&2; exit 1; }
[[ -f "${SCRIPT}" ]] || { echo "[ERROR] Missing script: ${SCRIPT}" >&2; exit 1; }

echo "### merge_grnboost2_seeds_v2.py: ${ORGAN} ###" >&2
echo "[INFO] Pattern : ${PATTERN}"  >&2
echo "[INFO] Output  : ${OUTFILE}"  >&2

srun python "${SCRIPT}" \
  --pattern      "${PATTERN}" \
  --output       "${OUTFILE}" \
  --corr-threshold "${CORR_THRESHOLD}" \
  --n_present    "${N_PRESENT}" \
  --present_rate "${PRESENT_RATE}"
```

---

## Step 3: Extract RBP Regulatory Modules (6 strategies)

Apply `scRBP_merge_grn_modules.py` to extract RBP modules from the consensus network using 6 merging strategies (combinations of top-N filtering and percentile thresholds).

```bash
# Submit for all tissues
TISSUES=(Bladder Blood Bone_Marrow Ear Fat Heart Kidney Large_Intestine
         Lung Lymph_Node Mammary Muscle Ovary Pancreas Prostate
         Salivary_Gland Skin Small_Intestine Spleen Stomach Testis Thymus)

for T in "${TISSUES[@]}"; do
  sbatch z_grn_getModules_for_28organs.sh "$T"
done
```

### z_grn_getModules_for_28organs.sh

```bash
#!/usr/bin/env bash
#SBATCH -J MERGE_MODULES_${ORGAN}
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail

ORGAN="$1"

# Paths
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/Palmer_hardwick_short_read/z_03_regulon_construction
SCRIPT="${TOOL}/scRBP_merge_grn_modules.py"
DATA_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/28organs_single_cell_atlas/z_GRNBoost2_${ORGAN}_30times
INFILE="${DATA_DIR}/merged_scRBP_GRNs_consensus.tsv"
OUTFILE="${DATA_DIR}/z_GRNBoost2_${ORGAN}_20Kcells_merged_scRBP_GRNs_consensus.tsv"

# Module extraction parameters
IMPORTANCE_THRESHOLD="0.005"
TOP_N_LIST="5,10,50"
TARGET_TOP_N="50"
PERCENTILE="0.75,0.90"

# Validation
mkdir -p logs
[[ -f "${SCRIPT}" ]] || { echo "[ERROR] Missing script: ${SCRIPT}" >&2; exit 1; }
[[ -f "${INFILE}" ]] || { echo "[ERROR] Missing input: ${INFILE}" >&2; exit 1; }

echo "### scRBP_merge_grn_modules.py: ${ORGAN} ###" >&2

srun python "${SCRIPT}" \
  --input                "${INFILE}" \
  --importance_threshold "${IMPORTANCE_THRESHOLD}" \
  --top_n_list           "${TOP_N_LIST}" \
  --target_top_n         "${TARGET_TOP_N}" \
  --percentile           "${PERCENTILE}" \
  --output_merged        "${OUTFILE}" \
  --verbose
```

---

## Step 4: Motif-Based Regulon Pruning (4 mRNA regions)

Run `scRBP_prune_ctxscore_mt_v4.py` across all 4 mRNA regions for each tissue. A helper function `run_region()` avoids code duplication across regions.

```bash
for T in "${TISSUES[@]}"; do
  sbatch z_grn_getRegulons_for_28organs.sh "$T"
done
```

### z_grn_getRegulons_for_28organs.sh

```bash
#!/usr/bin/env bash
#SBATCH -J PRUNE_REGULON_${ORGAN}
#SBATCH --mem=100G
#SBATCH --cpus-per-task=32
#SBATCH --time=10-00:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail

ORGAN="$1"
mkdir -p logs

# ---- Paths ----
PRUNE_TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/01_gene_level_final_25044motifs/Cluster_Buster_matrix_3UTR_gene_final
PRUNE_SCRIPT="${PRUNE_TOOL}/scRBP_prune_ctxscore_mt_v4.py"

ANNO_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/01_gene_level_final_25044motifs
MOTIF_RBP_LINKS="${ANNO_DIR}/11_616RBPs_20746motifs_links_annotation.csv"

MODULES_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/28organs_single_cell_atlas/z_GRNBoost2_${ORGAN}_30times
RBP_TARGETS="${MODULES_DIR}/z_GRNBoost2_${ORGAN}_20Kcells_merged_scRBP_GRNs_consensus.tsv"

BASE_RANK="${ANNO_DIR}"

# ---- Pruning parameters ----
RANK_THRESHOLD="1500"
AUC_THRESHOLD="0.05"
MIN_GENES="20"
NES_THRESHOLD="3.0"
N_JOBS="8"
CHUNKSIZE="1"

# ---- Validation ----
[[ -f "${PRUNE_SCRIPT}" ]]  || { echo "[ERROR] Missing: ${PRUNE_SCRIPT}" >&2; exit 1; }
[[ -f "${MOTIF_RBP_LINKS}" ]] || { echo "[ERROR] Missing: ${MOTIF_RBP_LINKS}" >&2; exit 1; }
[[ -f "${RBP_TARGETS}" ]]  || { echo "[ERROR] Missing: ${RBP_TARGETS}" >&2; exit 1; }

# ---- Helper: run pruning for one region ----
run_region() {
  local REGION="$1"
  local RANK_DIR="${BASE_RANK}/Cluster_Buster_matrix_${REGION}_gene_final"
  local SAVE_DIR="${MODULES_DIR}/Results_final_9metabolictrait_all7strategies_20Kcells_top1500_${REGION}"
  local RANK_FEATHER="${RANK_DIR}/motif_gene_rank_formatted.feather"

  [[ -d "${RANK_DIR}" ]]     || { echo "[ERROR] Missing: ${RANK_DIR}" >&2; exit 1; }
  [[ -f "${RANK_FEATHER}" ]] || { echo "[ERROR] Missing: ${RANK_FEATHER}" >&2; exit 1; }
  mkdir -p "${SAVE_DIR}"

  echo "### scRBP prune: ${ORGAN} / ${REGION} ###" >&2

  srun python "${PRUNE_SCRIPT}" \
    --rbp_targets      "${RBP_TARGETS}" \
    --motif_rbp_links  "${MOTIF_RBP_LINKS}" \
    --motif_target_ranks "${RANK_FEATHER}" \
    --save_dir         "${SAVE_DIR}" \
    --rank_threshold   "${RANK_THRESHOLD}" \
    --auc_threshold    "${AUC_THRESHOLD}" \
    --min_genes        "${MIN_GENES}" \
    --nes_threshold    "${NES_THRESHOLD}" \
    --n_jobs           "${N_JOBS}" \
    --chunksize        "${CHUNKSIZE}"
}

# ---- Run all 4 regions ----
run_region "3UTR"
run_region "5UTR"
run_region "CDS"
run_region "Introns"
```

---

## Step 5: Convert Regulons to GMT Format

### 5.1 Single-tissue conversion (example)

Convert pruned regulon results to GMT format (gene-symbol and Entrez-ID), retaining only regulons with ≥10 genes. Gene mapping uses NCBI38.gene.loc (~96.5% mapping rate).

```bash
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/01_gene_level_final_25044motifs/Cluster_Buster_matrix_3UTR_gene_final
ANNO=/mnt/isilon/gandal_lab/mayl/01_GWAS_tools/MAGMA

python ${TOOL}/scRBP_convert_format_v2.py \
  --input      scRBP_prune_matched_results.csv \
  --out-symbol scRBP_prune_matched_results_min10genes.symbol_default.gmt \
  --out-entrez scRBP_prune_matched_results_min10genes.entrez_default.gmt \
  --map-custom ${ANNO}/NCBI38.gene.loc \
  --min_genes 10 \
  --sanitize-rbp \
  --drop-unmapped-genes \
  --drop-empty-sets
```

### 5.2 Batch conversion across all tissues and regions

```bash
bash z_grn_run_convert_all_tissues_regions.sh
```

#### z_grn_run_convert_all_tissues_regions.sh

```bash
#!/usr/bin/env bash
set -euo pipefail

# Paths
BASE_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/28organs_single_cell_atlas
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/01_gene_level_final_25044motifs/Cluster_Buster_matrix_3UTR_gene_final
ANNO=/mnt/isilon/gandal_lab/mayl/01_GWAS_tools/MAGMA
SCRIPT=${TOOL}/scRBP_convert_format_v2.py

N_OK=0
N_SKIP=0

# Loop over all tissue directories
for TISSUE_DIR in "${BASE_DIR}"/z_GRNBoost2_*_30times; do
  [[ -d "$TISSUE_DIR" ]] || continue
  echo "========================================"
  echo "Tissue: $(basename $TISSUE_DIR)"
  echo "========================================"

  # Find all region output directories under this tissue
  shopt -s nullglob
  REGION_DIRS=("$TISSUE_DIR"/Results_final_*_top1500_*)
  shopt -u nullglob

  if [[ ${#REGION_DIRS[@]} -eq 0 ]]; then
    echo "  [WARN] No region folders found, skipping"
    continue
  fi

  for RDIR in "${REGION_DIRS[@]}"; do
    INPUT_FILE="$RDIR/scRBP_prune_matched_results.csv"
    if [[ ! -f "$INPUT_FILE" ]]; then
      echo "  [SKIP] Missing: $INPUT_FILE"
      ((N_SKIP+=1))
      continue
    fi

    OUT_SYMBOL="$RDIR/scRBP_prune_matched_results_min10genes.symbol_default.gmt"
    OUT_ENTREZ="$RDIR/scRBP_prune_matched_results_min10genes.entrez_default.gmt"
    echo "  -> Region: $(basename $RDIR)"

    python "$SCRIPT" \
      --input      "$INPUT_FILE" \
      --out-symbol "$OUT_SYMBOL" \
      --out-entrez "$OUT_ENTREZ" \
      --map-custom "$ANNO/NCBI38.gene.loc" \
      --min_genes 10 \
      --sanitize-rbp \
      --drop-unmapped-genes \
      --drop-empty-sets

    ((N_OK+=1))
  done
done

echo "========================================"
echo "DONE. Converted: $N_OK | Skipped: $N_SKIP"
echo "========================================"
```

---

## Step 6: Merge 4-Region GMT Files Per Tissue

Concatenate the per-region GMT files into a single tissue-level GMT, appending region tags (e.g., `RBP1_3UTR`) to distinguish regulons from different mRNA regions.

```bash
bash merge_regions_per_tissue.sh
```

#### merge_regions_per_tissue.sh

```bash
#!/usr/bin/env bash
set -euo pipefail

BASE_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/28organs_single_cell_atlas

for TISSUE_DIR in ${BASE_DIR}/z_GRNBoost2_*_30times; do
  [[ -d "$TISSUE_DIR" ]] || continue
  echo "Processing: $(basename $TISSUE_DIR)"

  OUT_FILE="$TISSUE_DIR/scRBP_allRegions_min10genes.symbol_default.gmt"
  > "$OUT_FILE"  # clear previous output

  for REGION_DIR in "$TISSUE_DIR"/Results_final_*_top1500_*; do
    [[ -d "$REGION_DIR" ]] || continue
    INPUT_GMT="$REGION_DIR/scRBP_prune_matched_results_min10genes.symbol_default.gmt"
    [[ -f "$INPUT_GMT" ]] || { echo "  [SKIP] Missing: $INPUT_GMT"; continue; }

    # Auto-detect region tag from directory name
    if   [[ "$REGION_DIR" == *"3UTR"* ]];   then REGION_TAG="3UTR"
    elif [[ "$REGION_DIR" == *"5UTR"* ]];   then REGION_TAG="5UTR"
    elif [[ "$REGION_DIR" == *"CDS"* ]];    then REGION_TAG="CDS"
    elif [[ "$REGION_DIR" == *"Intron"* ]]; then REGION_TAG="Introns"
    else REGION_TAG="Unknown"; fi

    # Append region tag to regulon names (column 1) and merge
    awk -v tag="$REGION_TAG" 'BEGIN{FS=OFS="\t"}{$1=$1"_"tag; print}' \
      "$INPUT_GMT" >> "$OUT_FILE"
  done

  echo "  -> Created: $OUT_FILE"
done

echo "All tissues finished."
```

---

## Step 7: Compute Cross-Tissue RBP Region Composition

Parse the merged GMT files across all tissues to build a tissue × RBP × region presence matrix, then compute per-RBP region composition proportions (how often each RBP's regulon appears in each mRNA region across tissues).

```bash
python z_grn_scRBP_region_composition.py \
  --base_dir /mnt/isilon/.../28organs_single_cell_atlas \
  --top_n 60 \
  --out_prefix scRBP_30tissues_regionComposition
```

**Outputs:**
- `*.rbp_counts_and_composition.csv` — per-RBP counts and proportions across tissues
- `*.tissue_by_rbp_presence.csv` — full tissue × RBP × region binary presence matrix
- `*.plot1_top{N}_stacked_composition.png` — top-N RBP stacked bar plot
- `*.plot2_cluster_heatmap_composition_K{K}.png` — clustered heatmap

### z_grn_scRBP_region_composition.py

```python
#!/usr/bin/env python3
"""
Compute RBP binding-region composition across tissues from merged GMT files.

For each tissue, parses the merged GMT (containing regulons tagged with _3UTR,
_5UTR, _CDS, _Introns) to build a binary presence matrix. Aggregates across
tissues to compute per-RBP region composition proportions, then generates
summary tables and visualizations (stacked bar plot, clustered heatmap).
"""

import os
import re
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list

REGIONS = ["3UTR", "5UTR", "CDS", "Introns"]


def infer_tissue_name(tissue_dir: str) -> str:
    """Extract tissue name from directory: z_GRNBoost2_Liver_30times -> Liver."""
    base = os.path.basename(tissue_dir.rstrip("/"))
    m = re.match(r"z_GRNBoost2_(.+)_30times$", base)
    return m.group(1) if m else base


def parse_merged_gmt_presence(gmt_path: str) -> pd.DataFrame:
    """
    Parse merged GMT file to build per-RBP region presence table.

    Each line: setName<TAB>desc<TAB>gene1<TAB>gene2...
    setName ends with _3UTR / _5UTR / _CDS / _Introns.

    Returns DataFrame with columns: rbp, 3UTR, 5UTR, CDS, Introns (0/1).
    """
    seen = set()
    rows = []

    with open(gmt_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            set_name = parts[0]

            # Identify region from suffix
            region = None
            for r in REGIONS:
                if set_name.endswith("_" + r):
                    region = r
                    break
            if region is None:
                continue

            # Extract RBP name (remove region suffix and optional _regulon)
            rbp = set_name[: -(len(region) + 1)]
            rbp = re.sub(r"_regulon$", "", rbp)

            key = (rbp, region)
            if key in seen:
                continue
            seen.add(key)
            rows.append({"rbp": rbp, "region": region, "present": 1})

    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=["rbp"] + REGIONS)

    pres = (
        df.pivot_table(index="rbp", columns="region", values="present",
                       aggfunc="max", fill_value=0)
        .reindex(columns=REGIONS, fill_value=0)
        .reset_index()
    )
    return pres


def cluster_order(mat: np.ndarray, method: str = "average",
                  metric: str = "euclidean"):
    """Return row indices ordered by hierarchical clustering."""
    if mat.shape[0] <= 2:
        return np.arange(mat.shape[0])
    Z = linkage(mat, method=method, metric=metric)
    return leaves_list(Z)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base_dir", required=True,
                    help="Parent folder containing z_GRNBoost2_*_30times")
    ap.add_argument("--pattern", default="z_GRNBoost2_*_30times",
                    help="Glob pattern under base_dir")
    ap.add_argument("--gmt_name",
                    default="scRBP_allRegions_min10genes.symbol_default.gmt",
                    help="Merged GMT filename inside each tissue folder")
    ap.add_argument("--top_n", type=int, default=50,
                    help="Top N RBPs for stacked bar plot")
    ap.add_argument("--out_prefix", default="scRBP_regionComposition",
                    help="Output filename prefix")
    args = ap.parse_args()

    # ---- Discover tissue directories ----
    tissue_dirs = sorted(glob.glob(os.path.join(args.base_dir, args.pattern)))
    if not tissue_dirs:
        raise SystemExit(f"[ERROR] No tissue dirs found: {args.base_dir}/{args.pattern}")

    # ---- Collect per-tissue presence matrices ----
    all_tissues = []
    all_rbps = set()
    tissue_presence = {}

    for tdir in tissue_dirs:
        tissue = infer_tissue_name(tdir)
        gmt_path = os.path.join(tdir, args.gmt_name)
        if not os.path.isfile(gmt_path):
            print(f"[WARN] Missing GMT, skip: {gmt_path}")
            continue

        pres = parse_merged_gmt_presence(gmt_path)
        pres["tissue"] = tissue
        all_tissues.append(tissue)
        all_rbps.update(pres["rbp"].tolist())
        tissue_presence[tissue] = pres

    all_tissues = sorted(list(set(all_tissues)))
    all_rbps = sorted(list(all_rbps))

    if not all_tissues or not all_rbps:
        raise SystemExit("[ERROR] No data parsed. Check base_dir / gmt_name.")

    # ---- Build full tissue × RBP × region presence table ----
    full_rows = []
    for tissue in all_tissues:
        pres = tissue_presence[tissue].set_index("rbp")
        pres = pres.reindex(all_rbps, fill_value=0)[REGIONS].astype(int)
        pres["tissue"] = tissue
        pres["rbp"] = pres.index
        full_rows.append(pres.reset_index(drop=True))

    tissue_rbp = pd.concat(full_rows, ignore_index=True)

    # ---- Aggregate across tissues: region counts (0–N_tissues) ----
    rbp_counts = tissue_rbp.groupby("rbp")[REGIONS].sum().reset_index()
    rbp_counts.rename(columns={r: f"nTissues_{r}" for r in REGIONS}, inplace=True)

    rbp_counts["nTissues_anyRegion"] = tissue_rbp.groupby("rbp").apply(
        lambda df: int((df[REGIONS].sum(axis=1) > 0).sum())
    ).values

    nT = len(all_tissues)
    rbp_counts["n_tissues_total"] = nT

    # ---- Compute region composition proportions (sum to 1 per RBP) ----
    count_cols = [f"nTissues_{r}" for r in REGIONS]
    rbp_counts["total_region_occurrences"] = rbp_counts[count_cols].sum(axis=1)

    for r in REGIONS:
        rbp_counts[f"prop_{r}"] = (
            rbp_counts[f"nTissues_{r}"]
            / rbp_counts["total_region_occurrences"].replace({0: np.nan})
        )

    # Filter RBPs with zero occurrences
    rbp_counts_use = rbp_counts[rbp_counts["total_region_occurrences"] > 0].copy()

    # ---- Save tables ----
    out_counts = f"{args.out_prefix}.rbp_counts_and_composition.csv"
    out_tissue = f"{args.out_prefix}.tissue_by_rbp_presence.csv"
    rbp_counts_use.to_csv(out_counts, index=False)
    tissue_rbp.to_csv(out_tissue, index=False)
    print(f"[OK] Wrote: {out_counts}")
    print(f"[OK] Wrote: {out_tissue}")

    # ---- Plot 1: Top-N stacked bar (composition proportions) ----
    top = (rbp_counts_use
           .sort_values("total_region_occurrences", ascending=False)
           .head(args.top_n)
           .set_index("rbp")[[f"prop_{r}" for r in REGIONS]])

    ax = top.plot(kind="bar", stacked=True, figsize=(14, 5))
    ax.set_ylabel("Region composition (sums to 1)")
    ax.set_xlabel(f"RBP (top {args.top_n} by total region-occurrences)")
    ax.set_title("RBP binding-region composition across tissues")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}.plot1_top{args.top_n}_stacked.png", dpi=200)
    plt.close()

    # ---- Plot 2: Clustered heatmap (RBPs with ≥K occurrences) ----
    K = 5
    hm = (rbp_counts_use[rbp_counts_use["total_region_occurrences"] >= K]
          .set_index("rbp")[[f"prop_{r}" for r in REGIONS]]
          .fillna(0.0))

    order = cluster_order(hm.values, method="average", metric="euclidean")
    hm_ord = hm.iloc[order]

    fig, ax = plt.subplots(figsize=(6, max(6, 0.16 * hm_ord.shape[0])))
    im = ax.imshow(hm_ord.values, aspect="auto", interpolation="nearest")
    ax.set_yticks(np.arange(hm_ord.shape[0]))
    ax.set_yticklabels(hm_ord.index, fontsize=6)
    ax.set_xticks(np.arange(len(REGIONS)))
    ax.set_xticklabels(REGIONS)
    ax.set_title(f"Clustered heatmap (RBPs with ≥{K} region-occurrences)")
    plt.colorbar(im, ax=ax, label="Composition proportion")
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}.plot2_heatmap_K{K}.png", dpi=200)
    plt.close()

    print("[OK] All plots saved.")


if __name__ == "__main__":
    main()
```

---

## Step 8: Visualize RBP Region Preference with Hierarchical Clustering

Generate publication-quality stacked bar plots and dendrogram-annotated horizontal bar plots showing RBP region-binding preferences, with optional hierarchical clustering into k groups.

[z_grn_plot_rbp_region_pref_v2.py](https://github.com/mayunlong89/scRBP_reproduce/blob/main/Analyses/transcript_region_preference_30Tissues/z_grn_plot_rbp_region_pref_counts_props_v2.py)

### Usage

```bash
# Plot all RBPs with hierarchical clustering
python z_grn_plot_rbp_region_pref_v2.py \
  --stats_csv scRBP_30tissues_regionComposition.rbp_counts_and_composition.csv \
  --prefix scRBP_RBPRegionPref \
  --min_tissues 5 \
  --mode all \
  --n_clusters 4 \
  --cluster_method ward

# Plot balanced selection (top 50 per dominant region, 200 RBPs total)
python z_grn_plot_rbp_region_pref_v2.py \
  --stats_csv scRBP_30tissues_regionComposition.rbp_counts_and_composition.csv \
  --prefix scRBP_RBPRegionPref \
  --min_tissues 5 \
  --mode balanced \
  --n_per_region 50 \
  --n_clusters 4 \
  --cluster_method ward
```

### z_grn_plot_rbp_region_pref_v2.py

```python
#!/usr/bin/env python3
"""
Visualize RBP binding-region preferences with hierarchical clustering.

Reads the per-RBP region composition table (from Step 7), filters by minimum
tissue count, performs hierarchical clustering on region proportions, and
produces a dendrogram + stacked horizontal bar plot + cluster annotation strip.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram

REGIONS = ["3UTR", "5UTR", "CDS", "Introns"]

REGION_COLORS = {
    "3UTR":    "#0B6B3A",
    "5UTR":    "#A6D854",
    "CDS":     "#2CB34A",
    "Introns": "#FFB000",
}

CLUSTER_COLORS = {
    1: "#4E79A7",  # blue
    2: "#E15759",  # red
    3: "#B07AA1",  # purple
    4: "#9C755F",  # brown
}


# ============================================================
# Utility functions
# ============================================================

def _require_cols(df, cols, context):
    miss = [c for c in cols if c not in df.columns]
    if miss:
        raise SystemExit(f"[ERROR] Missing columns in {context}: {miss}")


def _save_fig(fig, out_base):
    fig.savefig(out_base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(out_base + ".pdf", bbox_inches="tight")
    plt.close(fig)


def add_dominant_region(df):
    """Assign each RBP its dominant mRNA region based on max proportion."""
    prop_cols = [f"prop_{r}" for r in REGIONS]
    _require_cols(df, prop_cols, "stats_csv")
    out = df.copy()
    out["dominant_region"] = (out[prop_cols].idxmax(axis=1)
                              .str.replace("prop_", "", regex=False))
    return out


def sort_by_preference(df, metric="prop"):
    """Sort RBPs by dominant region, then by region-specific score."""
    df = df.copy()
    count_cols = [f"nTissues_{r}" for r in REGIONS]
    _require_cols(df, count_cols, "stats_csv")

    if "total_region_occurrences" not in df.columns:
        df["total_region_occurrences"] = df[count_cols].sum(axis=1)

    parts = []
    for r in REGIONS:
        sub = df[df["dominant_region"] == r].copy()
        if sub.empty:
            continue
        sort_col = f"prop_{r}" if metric == "prop" else f"nTissues_{r}"
        sub = sub.sort_values([sort_col, "total_region_occurrences"],
                              ascending=[False, False])
        parts.append(sub)

    return pd.concat(parts, ignore_index=True) if parts else df.iloc[0:0].copy()


def select_rbplist(df_sorted, mode, top_n, n_per_region):
    """
    Select RBPs for plotting.

    Modes:
      - 'all': return all RBPs
      - 'top': return top_n RBPs overall
      - 'balanced': return n_per_region RBPs per dominant region
    """
    if mode == "all":
        return df_sorted.copy()
    elif mode == "top":
        return df_sorted.head(top_n).copy()
    elif mode == "balanced":
        parts = [df_sorted[df_sorted["dominant_region"] == r].head(n_per_region)
                 for r in REGIONS]
        return pd.concat(parts, ignore_index=True)
    raise ValueError("mode must be: all, top, or balanced")


# ============================================================
# Hierarchical clustering
# ============================================================

def hierarchical_cluster(df, n_clusters=4, method="ward"):
    """
    Cluster RBPs by region composition proportions.

    Returns:
      - df with 'cluster' and 'plot_order' columns, reordered by dendrogram
      - Z: linkage matrix for dendrogram plotting
    """
    prop_cols = [f"prop_{r}" for r in REGIONS]
    _require_cols(df, prop_cols, "hierarchical_cluster")

    x = df[prop_cols].values.astype(float)
    Z = linkage(x, method=method, metric="euclidean")

    # Get dendrogram leaf order
    dendro = dendrogram(Z, no_plot=True)
    leaves = dendro["leaves"]

    # Assign cluster labels
    raw_clusters = fcluster(Z, t=n_clusters, criterion="maxclust")
    out = df.copy()
    out["cluster_raw"] = raw_clusters

    # Relabel clusters by dominant region for interpretability
    cluster_info = []
    for cl in sorted(out["cluster_raw"].unique()):
        sub = out[out["cluster_raw"] == cl]
        means = sub[prop_cols].mean()
        dom_region = means.idxmax().replace("prop_", "")
        cluster_info.append((cl, dom_region, means.max()))

    region_rank = {r: i for i, r in enumerate(REGIONS)}
    cluster_info.sort(key=lambda x: (region_rank[x[1]], -x[2]))
    relabel = {old: i + 1 for i, (old, _, _) in enumerate(cluster_info)}

    out["cluster"] = out["cluster_raw"].map(relabel)
    out = out.iloc[leaves].copy().reset_index(drop=True)
    out["plot_order"] = np.arange(out.shape[0])

    return out, Z


def summarize_clusters(df):
    """Compute mean region proportions per cluster."""
    prop_cols = [f"prop_{r}" for r in REGIONS]
    rows = []
    for cl in sorted(df["cluster"].unique()):
        sub = df[df["cluster"] == cl]
        means = sub[prop_cols].mean()
        rows.append({
            "cluster": cl,
            "n_rbps": sub.shape[0],
            **{f"mean_{c}": means[c] for c in prop_cols},
            "dominant_region": means.idxmax().replace("prop_", "")
        })
    return pd.DataFrame(rows)


# ============================================================
# Plotting
# ============================================================

def plot_horizontal_clustered(df, Z, out_base, title=None):
    """
    Produce a three-panel horizontal figure:
      Left: dendrogram | Middle: stacked horizontal bars | Right: cluster strip
    """
    prop_cols = [f"prop_{r}" for r in REGIONS]
    _require_cols(df, ["rbp", "cluster"] + prop_cols, "plot")

    n = df.shape[0]
    fig_h = max(8, 0.22 * n)
    fig = plt.figure(figsize=(10, fig_h))

    gs = fig.add_gridspec(nrows=1, ncols=3,
                          width_ratios=[1.7, 5.8, 0.4], wspace=0.05)
    ax_d = fig.add_subplot(gs[0, 0])  # dendrogram
    ax_b = fig.add_subplot(gs[0, 1])  # stacked bars
    ax_c = fig.add_subplot(gs[0, 2])  # cluster strip

    # ---- Dendrogram ----
    dendrogram(Z, orientation="left", no_labels=True,
               color_threshold=None, above_threshold_color="black", ax=ax_d)
    ax_d.invert_yaxis()
    ax_d.set_xticks([])
    ax_d.set_yticks([])
    for s in ["top", "right", "bottom", "left"]:
        ax_d.spines[s].set_visible(False)

    # ---- Stacked horizontal bars ----
    y = np.arange(n)
    left = np.zeros(n)
    for r in REGIONS:
        vals = df[f"prop_{r}"].values
        ax_b.barh(y, vals, left=left, color=REGION_COLORS[r],
                  edgecolor="none", height=0.85, label=r)
        left += vals

    ax_b.set_ylim(-0.5, n - 0.5)
    ax_b.invert_yaxis()
    ax_b.set_xlim(0, 1.0)
    ax_b.set_xlabel("Proportion")
    ax_b.set_yticks(y)
    ax_b.set_yticklabels(df["rbp"].tolist(), fontsize=7)
    if title:
        ax_b.set_title(title, fontsize=13, pad=10)

    # Region legend
    handles, labels = ax_b.get_legend_handles_labels()
    label_map = {"3UTR": "3\u2032UTR", "5UTR": "5\u2032UTR",
                 "CDS": "CDS", "Introns": "Introns"}
    ax_b.legend(handles, [label_map.get(l, l) for l in labels],
                title="Region", frameon=False,
                bbox_to_anchor=(1.02, 1.0), loc="upper left")

    # ---- Cluster annotation strip ----
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(-0.5, n - 0.5)
    ax_c.invert_yaxis()

    cluster_vals = df["cluster"].tolist()
    for i, cl in enumerate(cluster_vals):
        ax_c.add_patch(Rectangle((0, i - 0.425), 1, 0.85,
                                 facecolor=CLUSTER_COLORS.get(cl, "#999999"),
                                 edgecolor="none"))

    # Cluster boundary lines
    for i in range(1, n):
        if cluster_vals[i] != cluster_vals[i - 1]:
            ax_b.axhline(i - 0.5, color="black", lw=1.0)
            ax_c.axhline(i - 0.5, color="white", lw=1.2)

    # Cluster labels
    start = 0
    for i in range(1, n + 1):
        if i == n or cluster_vals[i] != cluster_vals[i - 1]:
            cl = cluster_vals[start]
            mid = (start + i - 1) / 2
            ax_c.text(0.5, mid, f"C{cl}", ha="center", va="center",
                      rotation=90, fontsize=9, fontweight="bold", color="white")
            start = i

    ax_c.set_xticks([])
    ax_c.set_yticks([])
    ax_c.set_ylabel("Cluster", rotation=270, labelpad=15)
    for s in ["top", "right", "bottom", "left"]:
        ax_c.spines[s].set_visible(False)

    plt.tight_layout()
    _save_fig(fig, out_base)


# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stats_csv", required=True,
                    help="Per-RBP region composition CSV from Step 7")
    ap.add_argument("--prefix", default="scRBP_RBPRegionPref")
    ap.add_argument("--min_tissues", type=int, default=5,
                    help="Minimum tissues an RBP must appear in")
    ap.add_argument("--mode", choices=["all", "top", "balanced"], default="balanced",
                    help="RBP selection mode")
    ap.add_argument("--top_n", type=int, default=200,
                    help="N RBPs when --mode top")
    ap.add_argument("--n_per_region", type=int, default=50,
                    help="N RBPs per dominant region when --mode balanced")
    ap.add_argument("--n_clusters", type=int, default=4,
                    help="Number of hierarchical clusters")
    ap.add_argument("--cluster_method",
                    choices=["ward", "average", "complete", "single"],
                    default="ward")
    args = ap.parse_args()

    # ---- Load and filter ----
    df = pd.read_csv(args.stats_csv)
    count_cols = [f"nTissues_{r}" for r in REGIONS]
    prop_cols = [f"prop_{r}" for r in REGIONS]
    _require_cols(df, ["rbp", "nTissues_anyRegion"] + count_cols + prop_cols, "stats_csv")

    df = df[df["nTissues_anyRegion"] >= args.min_tissues].copy()
    if df.empty:
        raise SystemExit(f"[ERROR] No RBPs with nTissues_anyRegion >= {args.min_tissues}")

    df = add_dominant_region(df)
    if "total_region_occurrences" not in df.columns:
        df["total_region_occurrences"] = df[count_cols].sum(axis=1)

    # ---- Select and cluster ----
    df_sorted = sort_by_preference(df, metric="prop")
    df_sel = select_rbplist(df_sorted, args.mode, args.top_n, args.n_per_region)
    df_clu, Z = hierarchical_cluster(df_sel, args.n_clusters, args.cluster_method)
    cluster_summary = summarize_clusters(df_clu)

    # ---- Output tag ----
    tag = f"min{args.min_tissues}.{args.mode}"
    if args.mode == "top":
        tag += f".top{args.top_n}"
    if args.mode == "balanced":
        tag += f".perRegion{args.n_per_region}"
    tag += f".k{args.n_clusters}.{args.cluster_method}"

    # ---- Plot and save ----
    plot_horizontal_clustered(
        df=df_clu, Z=Z,
        out_base=f"{args.prefix}.proportion_hclust.{tag}",
        title=(f"RBP region composition across tissues "
               f"(k={args.n_clusters}, {args.cluster_method})")
    )

    df_clu.to_csv(f"{args.prefix}.proportion_hclust.{tag}.csv", index=False)
    cluster_summary.to_csv(f"{args.prefix}.cluster_summary.{tag}.csv", index=False)

    print("[OK] Saved:")
    print(f"  {args.prefix}.proportion_hclust.{tag}.png/.pdf")
    print(f"  {args.prefix}.proportion_hclust.{tag}.csv")
    print(f"  {args.prefix}.cluster_summary.{tag}.csv")


if __name__ == "__main__":
    main()
```
