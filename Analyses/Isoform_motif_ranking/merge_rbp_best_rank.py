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

