#!/usr/bin/env bash
#SBATCH -J MERGE_MODULES_${ORGAN}
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail

# =========================
# 0) Args
# =========================
ORGAN="$1"

# =========================
# 1) Paths
# =========================
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/Palmer_hardwick_short_read/z_03_regulon_construction
SCRIPT="${TOOL}/scRBP_merge_grn_modules.py"

# 输入/输出（建议放到 organ 对应目录里，避免不同 organ 覆盖）
DATA_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/28organs_single_cell_atlas/z_GRNBoost2_${ORGAN}_30times
INFILE="${DATA_DIR}/merged_scRBP_GRNs_consensus.tsv"

# 你原来的输出名我保留风格，但默认加上 organ 防止覆盖
OUTFILE="${DATA_DIR}/z_GRNBoost2_${ORGAN}_20Kcells_merged_scRBP_GRNs_consensus.tsv"

# =========================
# 2) Params
# =========================
IMPORTANCE_THRESHOLD="0.005"
TOP_N_LIST="5,10,50"
TARGET_TOP_N="50"
PERCENTILE="0.75,0.90"

# =========================
# 3) Checks
# =========================
mkdir -p logs
[[ -f "${SCRIPT}" ]] || { echo "[ERROR] Missing script: ${SCRIPT}" >&2; exit 1; }
[[ -f "${INFILE}" ]] || { echo "[ERROR] Missing input: ${INFILE}" >&2; exit 1; }

echo "############# scRBP_merge_grn_modules.py: ${ORGAN} #############" >&2
echo "[INFO] Input   : ${INFILE}"  >&2
echo "[INFO] Output  : ${OUTFILE}"  >&2
echo "[INFO] Params  : importance=${IMPORTANCE_THRESHOLD}, top_n_list=${TOP_N_LIST}, target_top_n=${TARGET_TOP_N}, percentile=${PERCENTILE}" >&2

# =========================
# 4) Run
# =========================
srun python "${SCRIPT}" \
  --input "${INFILE}" \
  --importance_threshold "${IMPORTANCE_THRESHOLD}" \
  --top_n_list "${TOP_N_LIST}" \
  --target_top_n "${TARGET_TOP_N}" \
  --percentile "${PERCENTILE}" \
  --output_merged "${OUTFILE}" \
  --verbose
  
  
