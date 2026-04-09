#!/usr/bin/env bash
set -euo pipefail

# ===============================
# 1) Paths
# ===============================
BASE_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/28organs_single_cell_atlas

TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/01_gene_level_final_25044motifs/Cluster_Buster_matrix_3UTR_gene_final
ANNO=/mnt/isilon/gandal_lab/mayl/01_GWAS_tools/MAGMA
SCRIPT=$TOOL/scRBP_convert_format_v2.py

# 统计用
N_OK=0
N_SKIP=0

# ===============================
# 2) Loop tissues
# ===============================
for TISSUE_DIR in "${BASE_DIR}"/z_GRNBoost2_*_30times; do
  [[ -d "$TISSUE_DIR" ]] || continue
  echo "========================================"
  echo "Tissue: $TISSUE_DIR"
  echo "========================================"

  # 每个 tissue 下的 4 个 region 目录（你截图2那种 Results_final_*）
  shopt -s nullglob
  REGION_DIRS=("$TISSUE_DIR"/Results_final_*_top1500_*)
  shopt -u nullglob

  if [[ ${#REGION_DIRS[@]} -eq 0 ]]; then
    echo "⚠️  No region folders found under: $TISSUE_DIR"
    continue
  fi

  for RDIR in "${REGION_DIRS[@]}"; do
    INPUT_FILE="$RDIR/scRBP_prune_matched_results.csv"
    if [[ ! -f "$INPUT_FILE" ]]; then
      echo "⚠️  Skip (missing): $INPUT_FILE"
      ((N_SKIP+=1))
      continue
    fi

    # 用 region 名做输出前缀，避免覆盖
    REGION_NAME="$(basename "$RDIR")"
    OUT_SYMBOL="$RDIR/scRBP_prune_matched_results_min10genes.symbol_default.gmt"
    OUT_ENTREZ="$RDIR/scRBP_prune_matched_results_min10genes.entrez_default.gmt"

    echo "  -> Region: $REGION_NAME"
    echo "     Input : $INPUT_FILE"
    echo "     Output: $OUT_SYMBOL"
    echo "     Output: $OUT_ENTREZ"

    python "$SCRIPT" \
      --input "$INPUT_FILE" \
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
