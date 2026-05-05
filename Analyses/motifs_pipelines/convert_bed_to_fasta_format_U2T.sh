#!/bin/bash

# 输入文件
input_file="motif_sequence_forSpliceAid_F_v3_rm_repeat"
# 输出文件
output_file="motif_sequence_U2T.fasta"

# 清空或创建输出文件
> "$output_file"

declare -A name_count

# 读取输入文件的每一行
while read -r line; do
  # 分割行，获取描述和序列
  name=$(echo "$line" | awk '{print $1}')
  sequence=$(echo "$line" | awk '{print $2}')
  
  # 将序列中的U替换为T
  sequence=$(echo "$sequence" | tr 'U' 'T')
  
  # 检查名字是否已存在并更新计数
  if [[ -v name_count["$name"] ]]; then
    name_count["$name"]=$((name_count["$name"] + 1))
    new_name="${name}.${name_count["$name"]}"
  else
    name_count["$name"]=1
    new_name="${name}.1"
  fi

  # 将描述和序列写入输出文件，格式化为FASTA
  echo ">$new_name" >> "$output_file"
  echo "$sequence" >> "$output_file"
done < "$input_file"

echo "Conversion to FASTA format completed. Output file: $output_file"
