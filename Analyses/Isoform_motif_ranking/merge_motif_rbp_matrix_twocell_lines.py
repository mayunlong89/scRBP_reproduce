#!/usr/bin/env python3

import pandas as pd
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="合并两个 motif-by-RBP 矩阵，按 motif 行进行并集合并")
    parser.add_argument("-i1", "--input1", required=True, help="第一个 motif-by-RBP 矩阵的 CSV 文件路径")
    parser.add_argument("-i2", "--input2", required=True, help="第二个 motif-by-RBP 矩阵的 CSV 文件路径")
    parser.add_argument("-o", "--output", required=True, help="输出合并后的 CSV 文件路径")

    args = parser.parse_args()

    # 读取两个数据集
    try:
        df1 = pd.read_csv(args.input1, index_col=0)
        df2 = pd.read_csv(args.input2, index_col=0)
    except Exception as e:
        print(f"读取输入文件失败: {e}")
        sys.exit(1)

    # 按 index（motif）进行行合并，保留所有 motif（取并集）
    combined_df = pd.concat([df1, df2], axis=1, join='outer')

    # 保存输出文件
    try:
        combined_df.to_csv(args.output)
        print(f"合并后的矩阵已保存至: {args.output}")
    except Exception as e:
        print(f"保存输出文件失败: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

