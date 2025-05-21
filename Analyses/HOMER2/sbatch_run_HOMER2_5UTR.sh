#!/bin/bash
#SBATCH --job-name=5UTR_homer2_find   # Job name
#SBATCH --output=logs/homer2_find_5UTR.out  # Output log file for each array task
#SBATCH --error=logs/homer2_find_5UTR.err   # Error log file for each array task
#SBATCH --array=1-10             # Job array for 100 tasks
#SBATCH --ntasks=1                # Number of tasks per job
#SBATCH --cpus-per-task=1         # Number of CPU cores per task

# Load necessary modules
#module load homer/2.0  # Load the HOMER module, or adjust this line based on your environment

# 加载 conda 环境
source ~/.bashrc  # 如果 .bashrc 中包含 conda 初始化，这一步是必要的
conda activate pyscenic

# Set up paths and parameters
fasta_ref="/mnt/isilon/gandal_lab/mayl/reference/transcript_regions"
motif_dir="/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs/002_HOMER2_results/pwm_splites"  # Directory containing motif files
motif_file="${motif_dir}/motifs_part_${SLURM_ARRAY_TASK_ID}.pwm"
output_file="${motif_dir}/zzz_5UTR/homer2_5UTR_regions_targets_motifs_part_${SLURM_ARRAY_TASK_ID}.txt"

# Run homer2 find
homer2 find -i "$fasta_ref/5UTR_transcript_regions.fa" \
             -o "$output_file" \
             -m "$motif_file" \
             -mscore -p 1 -strand both
