


```bash

###------2025-04-15---

###----AUC---5utr HepG2
tool=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE
work_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions/z_ENCODE_HepG2_103RBPs_top0.05_5mingenes_adjustedMethods_none_5UTR
output_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions

python $tool/python_recovery_RBP_plot_final.py \
  -m $work_dir/motif_by_rbp_matrix_5mingenes.csv \
  -t $tool/encode_rbp_motif_links.csv \
  -r $output_dir/best_rank_output_5utr_HepG2.csv \
  -f $output_dir/recovery_curve_rank100motifs_5mingenes_5utr_HepG2.png \
  -a $output_dir/auc_5utr_HepG2.txt \
  --max_rank 100
  



###----AUC---5utr K562
tool=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE
work_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions/z_ENCODE_K562_120RBPs_top0.05_5mingenes_adjustedMethods_none_5UTR
output_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions

python $tool/python_recovery_RBP_plot_final.py \
  -m $work_dir/motif_by_rbp_matrix_5mingenes.csv \
  -t $tool/encode_rbp_motif_links.csv \
  -r $output_dir/best_rank_output_5utr_K562.csv \
  -f $output_dir/recovery_curve_rank100motifs_5mingenes_5utr_K562.png \
  -a $output_dir/auc_5utr_K562.txt \
  --max_rank 100
  
 
 
 ###----AUC---3utr HepG2
tool=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE
work_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions/z_ENCODE_HepG2_103RBPs_top0.05_5mingenes_adjustedMethods_none_3UTR
output_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions

python $tool/python_recovery_RBP_plot_final.py \
  -m $work_dir/motif_by_rbp_matrix_5mingenes.csv \
  -t $tool/encode_rbp_motif_links.csv \
  -r $output_dir/best_rank_output_3utr_HepG2.csv \
  -f $output_dir/recovery_curve_rank100motifs_5mingenes_3utr_HepG2.png \
  -a $output_dir/auc_3utr_HepG2.txt \
  --max_rank 100
  



###----AUC---3utr K562
tool=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE
work_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions/z_ENCODE_K562_120RBPs_top0.05_5mingenes_adjustedMethods_none_3UTR
output_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions

python $tool/python_recovery_RBP_plot_final.py \
  -m $work_dir/motif_by_rbp_matrix_5mingenes.csv \
  -t $tool/encode_rbp_motif_links.csv \
  -r $output_dir/best_rank_output_3utr_K562.csv \
  -f $output_dir/recovery_curve_rank100motifs_5mingenes_3utr_K562.png \
  -a $output_dir/auc_3utr_K562.txt \
  --max_rank 100
  
 
 
  ###----AUC---CDS HepG2
tool=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE
work_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions/z_ENCODE_HepG2_103RBPs_top0.05_5mingenes_adjustedMethods_none_CDS
output_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions

python $tool/python_recovery_RBP_plot_final.py \
  -m $work_dir/motif_by_rbp_matrix_5mingenes.csv \
  -t $tool/encode_rbp_motif_links.csv \
  -r $output_dir/best_rank_output_cds_HepG2.csv \
  -f $output_dir/recovery_curve_rank100motifs_5mingenes_cds_HepG2.png \
  -a $output_dir/auc_cds_HepG2.txt \
  --max_rank 100
  



###----AUC---CDS K562
tool=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE
work_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions/z_ENCODE_K562_120RBPs_top0.05_5mingenes_adjustedMethods_none_CDS
output_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions

python $tool/python_recovery_RBP_plot_final.py \
  -m $work_dir/motif_by_rbp_matrix_5mingenes.csv \
  -t $tool/encode_rbp_motif_links.csv \
  -r $output_dir/best_rank_output_cds_K562.csv \
  -f $output_dir/recovery_curve_rank100motifs_5mingenes_cds_K562.png \
  -a $output_dir/auc_cds_K562.txt \
  --max_rank 100
```

### python_recovery_RBP_plot_final.py
```python
#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import argparse
import sys
from sklearn.metrics import auc
def calculate_recovery_auc(rbp_best_rank: dict, max_rank: int = 100):
    """
    Input:
        rbp_best_rank: dictionary, format {RBP: BestRank}
        max_rank: cutoff rank threshold (e.g., 100)
    
    Output:
        recovery_curve: points on the recovery curve
        auc_score: area under the curve (normalized AUC between 0 and 1)
    """
    recovery_hits = np.zeros(max_rank)
    for rank in rbp_best_rank.values():
        if rank <= max_rank:
            recovery_hits[rank - 1:] += 1
    x = np.linspace(1, max_rank, max_rank)
    y = recovery_hits
    auc_score = auc(x / max_rank, y / y[-1]) if y[-1] > 0 else 0.0
    return y, auc_score
def main():
    parser = argparse.ArgumentParser(description="Compute and plot RBP recovery curve, and output AUC score")
    parser.add_argument("-m", "--matrix", required=True, help="Path to the input motif-by-RBP matrix CSV file")
    parser.add_argument("-t", "--truth", required=True, help="Path to the ground-truth RBP-motif mapping CSV file")
    parser.add_argument("-r", "--rank_output", default="rbp_best_motif_rank.csv", help="Output CSV file path for the best motif rank per RBP")
    parser.add_argument("-f", "--fig_output", default="rbp_recovery_curve_top100.png", help="Output PNG file path for the recovery curve plot")
    parser.add_argument("-a", "--auc_output", default="AUC_score.txt", help="Output file path for the AUC score")
    parser.add_argument("--max_rank", type=int, default=100, help="Maximum rank threshold for plotting (default: 100)")
    args = parser.parse_args()
    # Step 1: Read input data
    try:
        motif_rbp_df = pd.read_csv(args.matrix, index_col=0)
        truth_df = pd.read_csv(args.truth)
    except Exception as e:
        print(f"Failed to read file: {e}")
        sys.exit(1)
    # Step 2: Rank motifs for each RBP
    motif_rank_df = motif_rbp_df.rank(ascending=False, method="min")
    # Step 3: Get the best (minimum) rank of the true motif for each RBP
    rbp_best_rank = {}
    for _, row in tqdm(truth_df.iterrows(), total=truth_df.shape[0], desc="Ranking motifs"):
        rbp = row['RBP']
        motif = row['Motif']
        if rbp not in motif_rank_df.columns or motif not in motif_rank_df.index:
            continue
        rank = motif_rank_df.at[motif, rbp]
        if pd.isna(rank):
            continue
        rank = int(rank)
        if rbp not in rbp_best_rank:
            rbp_best_rank[rbp] = rank
        else:
            rbp_best_rank[rbp] = min(rbp_best_rank[rbp], rank)
    # Step 3.5: Save the best motif rank for each RBP
    rbp_rank_df = pd.DataFrame(list(rbp_best_rank.items()), columns=['RBP', 'BestRank'])
    rbp_rank_df.to_csv(args.rank_output, index=False)
    print(f"[✓] Best motif rank per RBP saved to: {args.rank_output}")
    # Step 4: Build recovery curve and compute AUC
    recovery_hits, auc_score = calculate_recovery_auc(rbp_best_rank, max_rank=args.max_rank)
    # Step 5: Save AUC
    with open(args.auc_output, 'w') as f:
        f.write(f"AUC_score (normalized) = {auc_score:.6f}\n")
    print(f"[✓] AUC score saved to: {args.auc_output}  value: {auc_score:.4f}")
    # Step 6: Plot
    plt.figure(figsize=(6, 5))
    plt.plot(range(1, args.max_rank + 1), recovery_hits, marker='o', markersize=3, label='YourMethod')
    plt.xlabel("Rank")
    plt.ylabel("Number of hits (Recovered RBPs)")
    plt.title("RBP Recovery Curve")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.fig_output, dpi=300)
    print(f"[✓] Recovery curve plot saved to: {args.fig_output}")
    plt.show()
if __name__ == "__main__":
    main()

```


### Merge all AUCs from two cell lines and four regions (3UTR, 5UTR, CDS, Introns)

```bash
#2025-04-15---merge best ranking for calculating AUC score

#Top 100 without 10 RBPs
tool=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE
work_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions
output_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions

python $tool/merge_rbp_best_rank.py \
    --input_files $work_dir/best_rank_output_5utr_HepG2.csv $work_dir/best_rank_output_5utr_K562.csv $work_dir/best_rank_output_3utr_HepG2.csv $work_dir/best_rank_output_3utr_K562.csv $work_dir/best_rank_output_cds_HepG2.csv $work_dir/best_rank_output_cds_K562.csv $work_dir/best_rank_output_Intron_HepG2.csv $work_dir/best_rank_output_Intron_K562.csv \
    --output merged_best_rank_5utr_3utr_cds_intron.csv \
    --fig_output merged_recovery_curve_5utr_3utr_cds_intron.png \
    --auc_output AUC_score_5utr_3utr_cds_intron.txt \
    --max_rank 100


#Top 100 with adding 10 RBPs
tool=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE
work_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions
output_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_test_four_regions

python $tool/merge_rbp_best_rank.py \
    --input_files $work_dir/best_rank_output_5utr_HepG2.csv $work_dir/best_rank_output_5utr_K562.csv $work_dir/best_rank_output_3utr_HepG2.csv $work_dir/best_rank_output_3utr_K562.csv $work_dir/best_rank_output_cds_HepG2.csv $work_dir/best_rank_output_cds_K562.csv $work_dir/best_rank_output_Intron_HepG2.csv $work_dir/best_rank_output_Intron_K562.csv $work_dir/best_rank_output_10RBPs_2.csv \
    --output merged_best_rank_5utr_3utr_cds_intron_10RBPs.csv \
    --fig_output merged_recovery_curve_5utr_3utr_cds_intron.png \
    --auc_output AUC_score_5utr_3utr_cds_intron.txt \
    --max_rank 100
```


### `merge_rbp_best_rank.py`

```python
#!/usr/bin/env python3

import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import auc

def load_best_rank(file_path):
    df = pd.read_csv(file_path)
    return dict(zip(df['RBP'], df['BestRank']))

def calculate_recovery_auc(rbp_best_rank: dict, max_rank: int = 100):
    """
    输入:
        rbp_best_rank: 字典，格式为 {RBP: BestRank}
        max_rank: 截止的排名阈值（如 100）
    
    输出:
        recovery_curve: Recovery 曲线上的点
        auc_score: 曲线下的面积 (normalized AUC between 0 and 1)
    """
    recovery_hits = np.zeros(max_rank)
    
    for rank in rbp_best_rank.values():
        if rank <= max_rank:
            recovery_hits[rank - 1:] += 1
    
    x = np.linspace(1, max_rank, max_rank)
    y = recovery_hits
    auc_score = auc(x / max_rank, y / y[-1]) if y[-1] > 0 else 0.0  # 防止除以0

    return y, auc_score

def main():
    parser = argparse.ArgumentParser(description="合并多个 RBP 最佳 motif 排名文件并绘制 recovery curve 和计算 AUC")
    parser.add_argument('--input_files', nargs='+', required=True, help="多个 RBP BestRank CSV 文件路径")
    parser.add_argument('--output', default='merged_best_rank.csv', help="输出合并后的 CSV 文件路径")
    parser.add_argument('--fig_output', default='merged_recovery_curve.png', help="合并后的 recovery 曲线图像路径")
    parser.add_argument('--auc_output', default='AUC_score.txt', help="输出 AUC 分数的 TXT 文件路径")
    parser.add_argument('--max_rank', type=int, default=100, help="最大排名阈值用于绘图和 AUC")

    args = parser.parse_args()

    # Step 1: 读取所有 RBP BestRank 字典
    all_ranks = [load_best_rank(f) for f in args.input_files]

    # Step 2: 合并，取最小 rank
    merged_rank = {}
    for rank_dict in all_ranks:
        for rbp, rank in rank_dict.items():
            if rbp not in merged_rank:
                merged_rank[rbp] = rank
            else:
                merged_rank[rbp] = min(merged_rank[rbp], rank)

    # Step 3: 保存合并结果
    merged_df = pd.DataFrame(list(merged_rank.items()), columns=['RBP', 'BestRank'])
    merged_df.to_csv(args.output, index=False)
    print(f"[✓] 合并结果保存至: {args.output}")

    # Step 4: 计算 recovery curve 和 AUC
    recovery_curve, auc_val = calculate_recovery_auc(merged_rank, max_rank=args.max_rank)
    print(f"[✓] Recovery AUC: {auc_val:.4f}")

    # Step 5: 保存 AUC 分数
    with open(args.auc_output, 'w') as f:
        f.write(f"AUC_score (normalized) = {auc_val:.6f}\n")
    print(f"[✓] AUC 分数保存至: {args.auc_output}")

    # Step 6: 绘图
    plt.figure(figsize=(6, 5))
    plt.plot(range(1, args.max_rank + 1), recovery_curve, marker='o', markersize=3, label='Merged')
    plt.xlabel("Rank")
    plt.ylabel("Number of hits (Recovered RBPs)")
    plt.title("Merged RBP Recovery Curve")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.fig_output, dpi=300)
    print(f"[✓] Recovery 曲线保存至: {args.fig_output}")
    plt.show()

if __name__ == "__main__":
    main()
```
