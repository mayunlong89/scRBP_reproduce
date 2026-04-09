#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import argparse
import sys
from sklearn.metrics import auc

def calculate_recovery_auc(rbp_hit_ranks: dict, max_rank: int = 100):
    """
    输入:
        rbp_hit_ranks: 字典，格式为 {RBP: BestRank}
        max_rank: 截止的排名阈值（如 100）

    输出:
        recovery_curve: Recovery 曲线上的点（回收的 RBP 数）
        auc_score: 曲线下面积（normalized AUC between 0 and 1）
    """
    recovery_hits = np.zeros(max_rank)

    for rank in rbp_hit_ranks.values():
        if rank <= max_rank:
            recovery_hits[rank - 1:] += 1

    x = np.linspace(1, max_rank, max_rank)
    y = recovery_hits
    auc_score = auc(x / max_rank, y / y[-1]) if y[-1] > 0 else 0.0

    return y, auc_score

def main():
    parser = argparse.ArgumentParser(description="计算并绘制 RBP recovery curve，并输出 AUC 值")
    parser.add_argument("-m", "--matrix", required=True, help="输入的 motif-by-RBP 矩阵 CSV 文件路径")
    parser.add_argument("-t", "--truth", required=True, help="真实的 motif-RBP 对应关系 CSV 文件路径")
    parser.add_argument("-r", "--rank_output", default="rbp_best_rank.csv", help="输出 RBP 最佳排名的 CSV 文件路径")
    parser.add_argument("-f", "--fig_output", default="rbp_recovery_curve_top100.png", help="保存 recovery curve 图像的 PNG 文件路径")
    parser.add_argument("-a", "--auc_output", default="AUC_score.txt", help="保存 AUC 值的文件路径")
    parser.add_argument("--max_rank", type=int, default=100, help="最大排名阈值用于绘图（默认: 100）")

    args = parser.parse_args()

    # Step 1: 读取输入数据
    try:
        motif_rbp_df = pd.read_csv(args.matrix, index_col=0)
        truth_df = pd.read_csv(args.truth)
    except Exception as e:
        print(f"读取文件失败: {e}")
        sys.exit(1)

    # Step 2: 对每个 motif 的 RBP 进行排名
    rbp_rank_df = motif_rbp_df.rank(axis=1, ascending=False, method="min")

    # Step 3: 获取每个 RBP 被其真实 motif 命中的最小排名
    rbp_best_rank = {}

    for _, row in tqdm(truth_df.iterrows(), total=truth_df.shape[0], desc="Ranking RBPs"):
        rbp = row['RBP']
        motif = row['Motif']

        if motif not in rbp_rank_df.index or rbp not in rbp_rank_df.columns:
            continue

        rank = rbp_rank_df.at[motif, rbp]
        if pd.isna(rank):
            continue

        rank = int(rank)
        if rbp not in rbp_best_rank:
            rbp_best_rank[rbp] = rank
        else:
            rbp_best_rank[rbp] = min(rbp_best_rank[rbp], rank)

    # Step 3.5: 保存每个 RBP 的最佳 motif 排名
    rbp_rank_df = pd.DataFrame(list(rbp_best_rank.items()), columns=['RBP', 'BestRank'])
    rbp_rank_df.to_csv(args.rank_output, index=False)
    print(f"[✓] RBP 的最佳 motif 排名已保存至: {args.rank_output}")

    # Step 4: 构建 recovery curve 并计算 AUC
    recovery_hits, auc_score = calculate_recovery_auc(rbp_best_rank, max_rank=args.max_rank)

    # Step 5: 保存 AUC
    with open(args.auc_output, 'w') as f:
        f.write(f"AUC_score (normalized) = {auc_score:.6f}\n")
    print(f"[✓] AUC 分数保存至: {args.auc_output}  值为: {auc_score:.4f}")

    # Step 6: 绘图
    plt.figure(figsize=(6, 5))
    plt.plot(range(1, args.max_rank + 1), recovery_hits, marker='o', markersize=3, label='Motif-to-RBP Recovery')
    plt.xlabel("Rank")
    plt.ylabel("Number of hits (Recovered RBPs)")
    plt.title("RBP Recovery Curve")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.fig_output, dpi=300)
    print(f"[✓] Recovery curve 图已保存至: {args.fig_output}")
    plt.show()

if __name__ == "__main__":
    main()




