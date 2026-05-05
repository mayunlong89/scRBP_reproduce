#!/bin/bash

# ================================
# Step 0: 设置参考注释文件路径
# ================================
ref="/mnt/isilon/gandal_lab/mayl/reference/gencode.v44.transcript_protein_coding_lncRNA.bed"

# ================================
# 批量处理所有 .txt 文件
# ================================
for file in *.txt; do
    base_name=$(basename "$file" .txt)
    final_bed="${base_name}_targets_isoforms_final.bed"

    # 如果最终输出已存在，则跳过
    if [[ -f "$final_bed" ]]; then
        echo "⏩ Skipping completed: $file"
        continue
    fi

    echo "🔄 Processing $file ..."

    # ------------------------------------------
    # Step 1: 提取 cluster_id + score
    # ------------------------------------------
    grep -v "MTF" "$file" | cut -f 4,5 > "${base_name}.score"

    # ------------------------------------------
    # Step 2: 提取 cluster_id + motif_id
    # ------------------------------------------
    cut -f 13 "$file" | grep -v "^-" | grep -v "^motif" | sed 's/__motif_/\t/g' > "${base_name}.ID"

    # ------------------------------------------
    # Step 3: 排序两个文件，为后续 join 做准备
    # ------------------------------------------
    sort -k1,1 "${base_name}.score" > "${base_name}.score.sorted"
    sort -k1,1 "${base_name}.ID"    > "${base_name}.ID.sorted"

    # ------------------------------------------
    # Step 4: 合并 motif_id + score（模拟 pandas merge left join）
    # ------------------------------------------
    awk 'NR==FNR {score[$1]=$2; next} {print $1 "\t" $2 "\t" score[$1]}' \
        "${base_name}.score.sorted" "${base_name}.ID.sorted" > "${base_name}.merged.tsv"

    # ------------------------------------------
    # Step 5: 提取 chr, start, end, motif_id, score（模拟 str.extract → BED）
    # ------------------------------------------
    awk -F '\t' '{
        split($1, b, "__");              # 去除 cluster 后缀
        split(b[1], a, /[:\-]/);         # 提取 chr:start-end
        print a[1] "\t" a[2] "\t" a[3] "\t" $2 "\t" $3
    }' "${base_name}.merged.tsv" | sort -u > "${base_name}.bed"

    # ------------------------------------------
    # Step 6: 使用 bedtools 注释基因信息
    # ------------------------------------------
    bedtools intersect -a "${base_name}.bed" -b "$ref" -wa -wb > "${base_name}_annotated.bed"

    # ------------------------------------------
    # Step 7: 提取 gene_name 字段
    # ------------------------------------------
    awk '{
    		split($9, iso, "|");
    		print $1, $2, $3, $4, $5, iso[1]
    		}' OFS="\t" "${base_name}_annotated.bed" > "${base_name}_annotated_clean.bed"
 

    # ------------------------------------------
    # Step 8: 对 motif-gene 保留最大 score（模拟 groupby + idxmax）
    # ------------------------------------------
    awk '
    {
        key = $4 OFS $6
        if (!(key in max) || $5 > max[key]) {
            max[key] = $5
            line[key] = $0
        }
    }
    END {
        for (k in line) print line[k]
    }' "${base_name}_annotated_clean.bed" > "$final_bed"

    # ------------------------------------------
    # Step 9: 清理所有中间文件
    # ------------------------------------------
    rm -f "${base_name}.score" \
          "${base_name}.ID" \
          "${base_name}.score.sorted" \
          "${base_name}.ID.sorted" \
          "${base_name}.merged.tsv" \
          "${base_name}.bed" \
          "${base_name}_annotated.bed" \
          "${base_name}_annotated_clean.bed"

    echo "✅ Finished processing $file."
done

echo "🎉 All .txt files have been processed!"
