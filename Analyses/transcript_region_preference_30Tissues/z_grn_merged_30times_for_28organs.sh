#!/usr/bin/env bash
#SBATCH -J MERGE_GRN_${ORGAN}
#SBATCH --mem=200G
#SBATCH --cpus-per-task=1
#SBATCH --time=10-00:00:00
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
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data
DATA=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/28organs_single_cell_atlas/z_GRNBoost2_${ORGAN}_30times

SCRIPT="${TOOL}/merge_grnboost2_seeds_v2.py"
PATTERN="${DATA}/z_GRNBoost2_${ORGAN}_TSP1_seed*_correlation_v7_scRBP_GRNs.tsv"
OUTFILE="${DATA}/merged_scRBP_GRNs_consensus.tsv"

# =========================
# 2) Params
# =========================
CORR_THRESHOLD="0.0"
N_PRESENT="10"
PRESENT_RATE="0.3"

# =========================
# 3) Checks
# =========================
mkdir -p logs
[[ -d "${DATA}" ]]   || { echo "[ERROR] Missing DATA dir: ${DATA}" >&2; exit 1; }
[[ -f "${SCRIPT}" ]] || { echo "[ERROR] Missing script: ${SCRIPT}" >&2; exit 1; }

echo "############# merge_grnboost2_seeds_v2.py: ${ORGAN} #############" >&2
echo "[INFO] Pattern : ${PATTERN}"  >&2
echo "[INFO] Output  : ${OUTFILE}"  >&2
echo "[INFO] Params  : corr=${CORR_THRESHOLD}, n_present=${N_PRESENT}, present_rate=${PRESENT_RATE}" >&2

# =========================
# 4) Run
# =========================
srun python "${SCRIPT}" \
  --pattern "${PATTERN}" \
  --output "${OUTFILE}" \
  --corr-threshold "${CORR_THRESHOLD}" \
  --n_present "${N_PRESENT}" \
  --present_rate "${PRESENT_RATE}"
