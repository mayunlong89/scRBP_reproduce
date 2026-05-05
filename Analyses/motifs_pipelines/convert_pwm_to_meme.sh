#!/bin/bash

# 检查是否传入了足够的参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

# 从命令行参数获取输入文件夹和输出文件夹路径
input_dir="$1"
output_dir="$2"

# 创建输出文件夹（如果不存在）
mkdir -p "$output_dir"

# 循环处理输入文件夹中的每个 PWM 文件
for pwm_file in "$input_dir"/*.pwm; do
    # 获取文件名（不带路径和扩展名）
    filename=$(basename -- "$pwm_file")
    name="${filename%.*}"

    # 设置输出文件路径
    output_file="$output_dir/${name}.meme"

    # 将 PWM 文件转换为 MEME 格式并保存
    chen2meme "$pwm_file" > "$output_file"
done

echo "所有 PWM 文件已转换为 MEME 格式并存储在 $output_dir 目录中"

