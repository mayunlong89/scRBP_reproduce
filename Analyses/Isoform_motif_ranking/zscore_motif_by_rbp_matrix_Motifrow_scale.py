
import pandas as pd
import numpy as np
import argparse

def row_zscore_normalization(input_file, output_file):
    # 读取数据
    data = pd.read_csv(input_file)
    data.set_index('Motif', inplace=True)

    # 对每一行（motif）进行 z-score 标准化
    normalized_data = data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

    # 保存结果
    normalized_data.to_csv(output_file)
    print(f"按 motif 行 Z-score 标准化完成，结果已保存到 {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="对 motif-by-RBP 矩阵按 motif 行做 Z-score 标准化")
    parser.add_argument("-i", "--input", required=True, help="输入文件路径（CSV格式）")
    parser.add_argument("-o", "--output", required=True, help="输出文件路径（CSV格式）")

    args = parser.parse_args()
    row_zscore_normalization(args.input, args.output)

