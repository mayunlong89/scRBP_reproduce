#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from sklearn.metrics import auc

def calculate_recovery_curve(best_ranks: pd.Series, max_rank: int = 100):
    """
    输入:
        best_ranks: pd.Series，RBP -> BestRank
        max_rank: 截止排名

    输出:
        recovery_curve: recovery 数组
        auc_score: normalized AUC
        recall: TP / total RBP
    """
    recovery_hits = np.zeros(max_rank)
    valid_ranks = best_ranks[best_ranks != 'NA'].astype(int)
    total_rbps = len(best_ranks)

    for rank in valid_ranks:
        if rank <= max_rank:
            recovery_hits[rank - 1:] += 1

    x = np.linspace(1, max_rank, max_rank)
    y = recovery_hits
    auc_score = auc(x / max_rank, y / y[-1]) if y[-1] > 0 else 0.0
    recall = (y[-1] / total_rbps) if total_rbps > 0 else 0.0

    return y, auc_score, recall

def main():
    parser = argparse.ArgumentParser(description="比较多个 motif_best_rank.csv 的 RBP recovery 曲线")
    parser.add_argument("-i", "--inputs", nargs="+", required=True, help="多个 motif_best_rank.csv 文件路径")
    parser.add_argument("-o", "--output", default="recovery_comparison.pdf", help="输出图像路径")
    parser.add_argument("--max_rank", type=int, default=100, help="最大排名 (默认: 100)")
    args = parser.parse_args()

    plt.figure(figsize=(7, 6))
    summary = []

    for input_file in args.inputs:
        df = pd.read_csv(input_file)
        recovery_curve, auc_score, recall = calculate_recovery_curve(df['BestRank'], args.max_rank)

        # 添加起点 (0, 0)
        x = [0] + list(range(1, args.max_rank + 1))
        y = [0] + list(recovery_curve)

        label = os.path.splitext(os.path.basename(input_file))[0] + f" (AUC={auc_score:.3f}, Recall={recall:.3f})"
        plt.plot(x, y, marker='o', markersize=2, label=label)

        summary.append({
            "File": os.path.basename(input_file),
            "AUC": round(auc_score, 4),
            "Recall": round(recall, 4)
        })

    # 图形美化
    plt.xlabel("Rank")
    plt.ylabel("Number of hits (Recovered RBPs)")
    plt.title("Comparison of RBP Recovery Curves")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"[✓] Recovery 曲线图已保存至: {args.output}")
    plt.show()

    # 保存 AUC + Recall 结果表
    summary_df = pd.DataFrame(summary)
    summary_file = os.path.splitext(args.output)[0] + "_metrics.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"[✓] Summary AUC + Recall 保存至: {summary_file}")

if __name__ == "__main__":
    main()

