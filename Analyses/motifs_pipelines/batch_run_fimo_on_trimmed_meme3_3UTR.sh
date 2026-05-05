#!/bin/bash
#SBATCH --job-name=multi_task_job
#SBATCH --output=logs/task_%A_%a.out
#SBATCH --error=logs/task_%A_%a.err
#SBATCH --time=20-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2
#SBATCH --exclusive


# 检查是否传入了足够的参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir_base>"
    exit 1
fi

# 从命令行参数获取输入文件夹、输出文件夹基础路径、参考序列路径和背景文件路径
input_dir="$1"
output_dir_base="$2"

# reference
ref_dir="/mnt/isilon/gandal_lab/mayl/reference/transcript_regions/"
bgfile_dir="/mnt/isilon/gandal_lab/mayl/reference/"

# 设置参考序列文件和背景文件路径
ref_fasta="$ref_dir/3UTR_transcript_regions.fa"
bg_file="$bgfile_dir/background.meme"

# 循环处理输入文件夹中的每个 trimmed MEME 文件
for meme_file in "$input_dir"/*_trimmed.meme; do
    # 获取文件名（不带路径和扩展名）
    filename=$(basename -- "$meme_file")
    name="${filename%.*}"

    # 设置输出文件夹路径
    output_dir="$output_dir_base/fimo_output_${name}"

    # 使用 srun 投放到集群节点上运行 fimo，并将任务放到后台
    srun --pty fimo --verbosity 2 --max-stored-scores 1000000 --thresh 0.001 --bgfile "$bg_file" --oc "$output_dir" "$meme_file" "$ref_fasta" &

done

# 等待所有后台任务完成
wait

echo "所有 trimmed MEME 文件的 FIMO 扫描任务已经提交至集群节点并运行完成。"

#------------------------------------------------------------------------------------------------------

# & 符号: 使得每个 srun 命令在后台运行，从而能够同时提交多个任务。
# wait 命令: 脚本中的 wait 命令会等待所有后台任务完成后再继续执行脚本后续的代码。这确保了所有的 FIMO 任务都完成后，脚本才会结束。
