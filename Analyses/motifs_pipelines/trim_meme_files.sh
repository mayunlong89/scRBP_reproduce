#!/bin/bash

# 检查是否传入了足够的参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

# 从命令行参数获取输入文件夹、输出文件夹路径和universal_loc路径
input_dir="$1"
output_dir="$2"
universal_loc=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/01_snakemake_pipeline/

# 创建输出文件夹（如果不存在）
mkdir -p "$output_dir"

# 循环处理输入文件夹中的每个 MEME 文件
for meme_file in "$input_dir"/*.meme; do
    # 获取文件名（不带路径和扩展名）
    filename=$(basename -- "$meme_file")
    name="${filename%.*}"

    # 设置输出文件路径
    output_trimmed_file="$output_dir/${name}_trimmed.meme"

    # 使用Rscript调用universalmotif.r进行trimming处理
    Rscript "$universal_loc/universalmotif.r" "$meme_file" "$output_trimmed_file"
done
 
echo "所有 MEME 文件已经过 trimming 处理并存储在 $output_dir 目录中"


