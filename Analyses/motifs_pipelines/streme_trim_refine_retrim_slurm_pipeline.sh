#!/bin/bash
#SBATCH --job-name=fasta_processing   # 任务名称
#SBATCH --output=logs/%x_%A_%a.out    # 标准输出日志
#SBATCH --error=logs/%x_%A_%a.err     # 标准错误日志
#SBATCH --array=1-371             # 定义任务数组范围，总共371个任务，每次最多运行30个
#SBATCH --cpus-per-task=5            # 每个任务使用的CPU核数
#SBATCH --mem=10G                      # 每个任务使用的内存
#SBATCH --time=5:00:00                # 每个任务的时间限制

# 设置输入输出路径以及Snakemake路径
fasta_folder="$1"
output_prefix="$2"
snakemake_loc=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/01_snakemake_pipeline/

# 获取所有 .fasta 文件的列表
fasta_files=($(ls "$fasta_folder"/*.fasta))

# 计算当前任务处理的文件索引
file_idx=$((SLURM_ARRAY_TASK_ID - 1))

# 获取当前任务需要处理的 .fasta 文件
fasta_file=${fasta_files[$file_idx]}

# 提取文件名（不带路径和扩展名）
filename=$(basename "$fasta_file" .fasta)

# 确保输出目录存在
mkdir -p "$output_prefix"

echo "Processing $fasta_file..."

# Step 1: Find motifs using STREME
streme_dir="${output_prefix}/motifs/${filename}_motifs_dir"
mkdir -p "$streme_dir"
streme --oc "$streme_dir" --p "$fasta_file" --minw 4 --maxw 12 --thresh 0.05

# Step 2: Trim motifs
trimmed_meme="${output_prefix}/trimmed_motifs/${filename}_trimmed_meme_file.meme"
mkdir -p "$(dirname "$trimmed_meme")"
Rscript "$snakemake_loc/universalmotif.r" "$streme_dir/streme.txt" "$trimmed_meme"

# Step 3: Refine motifs using BaMMmotif
refine_dir="${output_prefix}/refined_motifs/${filename}_refine"
mkdir -p "$refine_dir"
BaMMmotif "$refine_dir" "$fasta_file" --PWMFile "$trimmed_meme" --EM -k 2 --advanceEM --saveBaMMs

# Step 4: Convert .ihbcp files to .meme files
refined_meme_dir="${refine_dir}/refine_meme"
mkdir -p "$refined_meme_dir"
bash "$snakemake_loc/convert_ihbcp_to_meme.sh" "$refine_dir" "$refined_meme_dir"

# Step 5: Merge MEME files into PWM format
merged_meme_prefix="${output_prefix}/merged_motifs/${filename}_merged_meme2pwm.pwm"
mkdir -p "$(dirname "$merged_meme_prefix")"
python "$snakemake_loc/merge_meme_files_into_PWM.py" "$refined_meme_dir" "$merged_meme_prefix"

# Step 6: Convert PWM file to MEME format using chen2meme
chen2meme "${merged_meme_prefix}" > "${merged_meme_prefix}_all_memes.meme"

# Step 7: Repeat trimming on the refined MEME file
Rscript "$snakemake_loc/universalmotif.r" "${merged_meme_prefix}_all_memes.meme" "${output_prefix}/${filename}_trimmed_refined_bamm2meme_file.meme"

echo "Completed processing for $fasta_file"

