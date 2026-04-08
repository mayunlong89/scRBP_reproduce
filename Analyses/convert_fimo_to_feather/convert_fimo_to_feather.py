import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import OrderedDict
import feather
from IPython.core.interactiveshell import InteractiveShell
import argparse
import os

InteractiveShell.ast_node_interactivity = "all"

def process_files(input_file, input_dir, output_dir):
    # Step 1: Read and process the fimo_target data
    contents = {}
    with open(os.path.join(input_dir, input_file)) as r:
        r.readline()  # Skip the header        
        line = r.readline()
        while line:
            arr = line.strip().split("\t")
            m_n = arr[3]
            g_n = arr[-2]
            rank = int(arr[-1])
            key = f"{m_n}\t{g_n}"
            if key in contents:
                contents[key] = min(rank, contents[key])
            else:
                contents[key] = rank
            line = r.readline()

    # Save the results to a txt file
    min_output_path = os.path.join(output_dir, "motif_gene_p_min.txt")
    with open(min_output_path, 'w') as w:
        for k, v in contents.items():
            w.write(f"{k}\t{v}\n")

    # Step 2: Read the saved minimum rank data and pivot
    df = pd.read_csv(min_output_path, header=None, sep="\t")
    df.columns = ["motif", "gene", "rank"]
    df_pivot = df.pivot(index='gene', columns='motif', values='rank')
    df_pivot = df_pivot.fillna(100000)  # Fill missing values

    # Save the pivoted results to a txt file
    find_output_path = os.path.join(output_dir, 'motif_gene_p_find.txt')
    df_pivot.to_csv(find_output_path, sep='\t', header=True, index=True)

    # Step 3: Rank the genes based on the motifs and save the results in .txt and .feather formats
    df_p_min = pd.read_csv(find_output_path, header=0, sep="\t", index_col=0)
    start_order = df_p_min.index.tolist()
    df_rank = pd.DataFrame()

    for col in df_p_min.columns:
        df_test_sort_b = df_p_min.sort_values(by=col, ascending=True)
        ind_index = OrderedDict()
        i = 1
        for ww in df_test_sort_b.index:
            ind_index[ww] = i
            i += 1
        rank_0 = [ind_index[ww] for ww in start_order]
        df_rank[col] = rank_0

    df_rank.index = start_order
    rank_output_path = os.path.join(output_dir, "motif_gene_rank.txt")
    df_rank.T.to_csv(rank_output_path, sep="\t", header=True, index=True)
    
    rank_output_path2 = os.path.join(output_dir, "motif_gene_rank.csv")
    df_rank.T.to_csv(rank_output_path2, sep=",", header=True, index=True)

    # Save the transposed rank DataFrame as a Feather file
    feather.write_dataframe(df_rank.T, os.path.join(output_dir, "motif_gene_rank.feather"))

def main():
    parser = argparse.ArgumentParser(description="Process fimo target data and generate output files.")
    parser.add_argument('--input_file', type=str, required=True, help="Name of the input file (e.g., fimo_best_site_narrowPeak_targets_rankedByPvalue_ranked.tsv).")
    parser.add_argument('--input_dir', type=str, required=True, help="Directory where the input file is located.")
    parser.add_argument('--output_dir', type=str, required=True, help="Directory where the output files will be saved.")
    
    args = parser.parse_args()
    
    # Process the files with the provided arguments
    process_files(args.input_file, args.input_dir, args.output_dir)

if __name__ == "__main__":
    main()

