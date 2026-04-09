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
LOOM_FILES=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/28organs_single_cell_atlas/loom_dir
GRN_FILES=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/28organs_single_cell_atlas/z_GRNBoost2_${ORGAN}_30times

srun python $TOOL/scRBP_infer_v7_cor.py \
  --matrix $LOOM_FILES/${ORGAN}_TSP1_30_v2.loom \
  --rbp_list $TOOL/rbp_list_616RBPs.tsv \
  --output $GRN_FILES/z_GRNBoost2_${ORGAN}_TSP1_seed${SEED}_correlation_v7.tsv \
  --n_workers 10 \
  --batch_size 5 \
  --threshold 0.03 \
  --method grnboost2 \
  --correlation True \
  --seed $SEED \
  --log $GRN_FILES/z_GRNBoost2_${ORGAN}_TSP1_seed${SEED}_correlation_v7.log

