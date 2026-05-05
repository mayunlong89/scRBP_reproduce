#!/bin/bash

# 检查是否提供了两个参数：输入 FASTA 文件夹和输出目录前缀
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <fasta_folder> <output_prefix>"
  exit 1
fi

# 获取输入的文件夹路径和输出前缀
fasta_folder="$1"
output_prefix="$2"

# 定义 Snakemake 工具文件路径
snakemake_loc=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/01_snakemake_pipeline/

# 确保输出目录存在
mkdir -p "$output_prefix"

# 遍历输入文件夹中的每个 .fasta 文件
for fasta_file in "$fasta_folder"/*.fasta; do
  # 提取文件名（不带路径和扩展名）
  filename=$(basename "$fasta_file" .fasta)
  
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
  BaMMmotif "$refine_dir" "$fasta_file" --PWMFile "$trimmed_meme" --EM -k 2  --advanceEM --extend 2 2 --saveBaMMs

  # Step 4: Convert .ihbcp files to .meme files
  refined_meme_dir="${refine_dir}/refine_meme"
  mkdir -p "$refined_meme_dir"
  bash "$snakemake_loc/convert_ihbcp_to_meme.sh" "$refine_dir" "$refined_meme_dir"


  # Step 5: Merge MEME files into PWM format
  merged_meme_prefix="${output_prefix}/merged_motifs/${filename}_merged_meme2pwm.pwm"
  mkdir -p "$(dirname "$merged_meme_prefix")"
  python "$snakemake_loc/merge_meme_files_into_PWM.py" "$refined_meme_dir" "$merged_meme_prefix"
  
  # 确保生成的 PWM 文件存在
  if [ ! -f "${merged_meme_prefix}" ]; then
    echo "Error: Failed to generate PWM file at ${merged_meme_prefix}"
    exit 1
  fi

  # Step 6: Convert PWM file to MEME format using chen2meme
  chen2meme "${merged_meme_prefix}" > "${merged_meme_prefix}_all_memes.meme"

	# 确保 MEME 文件成功生成
  if [ ! -f "${merged_meme_prefix}_all_memes.meme" ]; then
    echo "Error: Failed to generate MEME file at ${merged_meme_prefix}_all_memes.meme"
    exit 1
  fi

  # Step 7: Repeat trimming on the refined MEME file
  Rscript "$snakemake_loc/universalmotif.r" "${merged_meme_prefix}_all_memes.meme" "${output_prefix}/${filename}_trimmed_refined_bamm2meme_file.meme"

  echo "Completed processing for $fasta_file"
done
