#!/usr/bin/env bash
#SBATCH -J PRUNE_REGULON_${ORGAN}
#SBATCH --mem=100G
#SBATCH --cpus-per-task=32
#SBATCH --time=10-00:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail

# =========================
# 0) Args
# =========================
ORGAN="$1"

# =========================
# 1) Common paths
# =========================
mkdir -p logs

# scRBP prune script 所在目录（你之前一直用的那个工具目录）
PRUNE_TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/01_gene_level_final_25044motifs/Cluster_Buster_matrix_3UTR_gene_final
PRUNE_SCRIPT="${PRUNE_TOOL}/scRBP_prune_ctxscore_mt_v4.py"

# motif-RBP link annotation（对所有区域通用）
ANNO_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/01_gene_level_final_25044motifs
MOTIF_RBP_LINKS="${ANNO_DIR}/11_616RBPs_20746motifs_links_annotation.csv"

# GRN modules (Step3 输出所在目录)
MODULES_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/28organs_single_cell_atlas/z_GRNBoost2_${ORGAN}_30times
RBP_TARGETS="${MODULES_DIR}/z_GRNBoost2_${ORGAN}_20Kcells_merged_scRBP_GRNs_consensus.tsv"

# =========================
# 2) Params (keep same as your original)
# =========================
RANK_THRESHOLD="1500"
AUC_THRESHOLD="0.05"
MIN_GENES="20"
NES_THRESHOLD="3.0"
N_JOBS="8"
CHUNKSIZE="1"

# =========================
# 3) Checks
# =========================
[[ -f "${PRUNE_SCRIPT}" ]] || { echo "[ERROR] Missing prune script: ${PRUNE_SCRIPT}" >&2; exit 1; }
[[ -f "${MOTIF_RBP_LINKS}" ]] || { echo "[ERROR] Missing motif_rbp_links: ${MOTIF_RBP_LINKS}" >&2; exit 1; }
[[ -f "${RBP_TARGETS}" ]] || { echo "[ERROR] Missing rbp_targets: ${RBP_TARGETS}" >&2; exit 1; }

# =========================
# 4) Helper: run one region
# =========================
run_region() {
  local REGION="$1"          # e.g. 3UTR / 5UTR / CDS / Introns
  local RANK_DIR="$2"        # motif ranking directory for that region
  local SAVE_DIR="$3"        # output directory

  local RANK_FEATHER="${RANK_DIR}/motif_gene_rank_formatted.feather"

  echo "############# scRBP getRegulon: ${ORGAN}, '${REGION}' #############" >&2
  echo "[INFO] rbp_targets       : ${RBP_TARGETS}" >&2
  echo "[INFO] motif_rbp_links   : ${MOTIF_RBP_LINKS}" >&2
  echo "[INFO] motif_target_ranks: ${RANK_FEATHER}" >&2
  echo "[INFO] save_dir          : ${SAVE_DIR}" >&2

  [[ -d "${RANK_DIR}" ]] || { echo "[ERROR] Missing rank dir: ${RANK_DIR}" >&2; exit 1; }
  [[ -f "${RANK_FEATHER}" ]] || { echo "[ERROR] Missing rank feather: ${RANK_FEATHER}" >&2; exit 1; }
  mkdir -p "${SAVE_DIR}"

  srun python "${PRUNE_SCRIPT}" \
    --rbp_targets "${RBP_TARGETS}" \
    --motif_rbp_links "${MOTIF_RBP_LINKS}" \
    --motif_target_ranks "${RANK_FEATHER}" \
    --save_dir "${SAVE_DIR}" \
    --rank_threshold "${RANK_THRESHOLD}" \
    --auc_threshold "${AUC_THRESHOLD}" \
    --min_genes "${MIN_GENES}" \
    --nes_threshold "${NES_THRESHOLD}" \
    --n_jobs "${N_JOBS}" \
    --chunksize "${CHUNKSIZE}"
}

# =========================
# 5) Run 4 regions
# =========================
BASE_RANK=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/01_gene_level_final_25044motifs

run_region "3UTR"   "${BASE_RANK}/Cluster_Buster_matrix_3UTR_gene_final"    "${MODULES_DIR}/Results_final_9metabolictrait_all7strategies_20Kcells_top1500_3UTR"
run_region "5UTR"   "${BASE_RANK}/Cluster_Buster_matrix_5UTR_gene_final"    "${MODULES_DIR}/Results_final_9metabolictrait_all7strategies_20Kcells_top1500_5UTR"
run_region "CDS"    "${BASE_RANK}/Cluster_Buster_matrix_CDS_gene_final"     "${MODULES_DIR}/Results_final_9metabolictrait_all7strategies_20Kcells_top1500_CDS"
run_region "Introns" "${BASE_RANK}/Cluster_Buster_matrix_Introns_gene_final" "${MODULES_DIR}/Results_final_9metabolictrait_all7strategies_20Kcells_top1500_Introns"


