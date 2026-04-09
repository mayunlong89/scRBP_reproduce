#!/usr/bin/env python3

import pandas as pd
import argparse
import sys

def main():
    # 设置参数
    parser = argparse.ArgumentParser(description="将 motif-RBP 表格转换为 motif × RBP 的矩阵")
    parser.add_argument("-i", "--input", required=True, help="输入的 CSV 文件路径（包含 Motif, RBP, AUC_Score 列）")
    parser.add_argument("-o", "--output", required=True, help="输出的矩阵 CSV 文件路径")

    args = parser.parse_args()

    # 读取文件
    try:
        df = pd.read_csv(args.input)
    except Exception as e:
        print(f"读取输入文件失败: {e}")
        sys.exit(1)

    # 检查所需列是否存在
    required_columns = {"Motif", "RBP", "AUC_Score"}
    if not required_columns.issubset(df.columns):
        print(f"错误: 输入文件缺少必要的列: {required_columns}")
        sys.exit(1)

    # 创建透视表
    nes_matrix = df.pivot(index="Motif", columns="RBP", values="AUC_Score")

    # 保存为输出文件
    try:
        nes_matrix.to_csv(args.output)
        print(f"转换完成，矩阵保存至: {args.output}")
    except Exception as e:
        print(f"保存输出文件失败: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()




