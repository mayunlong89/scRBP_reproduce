"""
Optimized version of convert_fimo_to_feather_enhanced.py

Improvements over v1:
  1. No intermediate txt/csv files — everything stays in memory
  2. Chunked reading via pd.read_csv(chunksize=) for 100G+ files
  3. pandas groupby().min() replaces Python dict line-by-line accumulation
  4. df.rank(method='first') replaces per-column sort_values + OrderedDict
  5. Only outputs .feather by default (add --also_csv if needed)
  6. Uses float32/int32 to halve memory vs default float64
"""

import argparse
import os
import pandas as pd
import pyarrow.feather as feather


def process_files(input_file, input_dir, output_dir,
                  chunksize=5_000_000, fill_value=100000, also_csv=False):
    os.makedirs(output_dir, exist_ok=True)

    filepath = os.path.join(input_dir, input_file)

    # --- Determine column indices from the first data row ---
    first_row = pd.read_csv(filepath, sep="\t", nrows=1, header=None)
    ncol = first_row.shape[1]
    motif_col = 3          # column 4 (0-indexed: 3)
    gene_col  = ncol - 2   # second-to-last column
    rank_col  = ncol - 1   # last column

    print(f"Input: {filepath}")
    print(f"Detected {ncol} columns. "
          f"motif=col{motif_col}, gene=col{gene_col}, rank=col{rank_col}")

    # =====================================================================
    # Step 1: Chunked read + per-chunk groupby min
    #   - Only reads 3 columns (motif, gene, rank), ignores everything else
    #   - Each chunk is reduced to unique (motif, gene) pairs with min rank
    # =====================================================================
    chunk_results = []
    total_rows = 0

    reader = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        skiprows=1,          # skip the header line
        usecols=[motif_col, gene_col, rank_col],
        names=["motif", "gene", "rank"],
        dtype={"motif": "category", "gene": "str", "rank": "int32"},
        chunksize=chunksize,
    )

    for i, chunk in enumerate(reader):
        total_rows += len(chunk)
        # keep min rank per (motif, gene) within this chunk
        chunk = chunk.groupby(["motif", "gene"], observed=True,
                              as_index=False)["rank"].min()
        chunk_results.append(chunk)
        if (i + 1) % 10 == 0:
            print(f"  ... processed {total_rows:,} rows "
                  f"({len(chunk_results)} chunks accumulated)")

    print(f"Total rows read: {total_rows:,}")

    # =====================================================================
    # Step 2: Global min rank across all chunks
    # =====================================================================
    print("Merging chunks → global min per (motif, gene) ...")
    df = pd.concat(chunk_results, ignore_index=True)
    del chunk_results
    df = df.groupby(["motif", "gene"], observed=True,
                    as_index=False)["rank"].min()
    print(f"Unique (motif, gene) pairs: {len(df):,}")

    # =====================================================================
    # Step 3: Pivot → gene x motif matrix, fill missing with fill_value
    # =====================================================================
    print("Pivoting to gene x motif matrix ...")
    mat = df.pivot(index="gene", columns="motif", values="rank")
    del df
    mat = mat.fillna(fill_value).astype("float32")
    print(f"Matrix shape: {mat.shape[0]:,} genes x {mat.shape[1]:,} motifs")

    # =====================================================================
    # Step 4: Rank each motif column (ascending, ties → first occurrence)
    #   Equivalent to original: sort_values per col + OrderedDict numbering
    # =====================================================================
    print("Computing per-motif rankings ...")
    rank_mat = mat.rank(axis=0, method="first", ascending=True).astype("int32")
    del mat

    # =====================================================================
    # Step 5: Transpose to motif x gene → write .feather
    # =====================================================================
    print("Transposing and writing feather ...")
    out_df = rank_mat.T
    del rank_mat
    out_df.index.name = "motifs"
    out_df = out_df.reset_index()

    feather_path = os.path.join(output_dir, "motif_gene_rank.feather")
    feather.write_feather(out_df, feather_path)
    print(f"Wrote {feather_path}")

    if also_csv:
        csv_path = os.path.join(output_dir, "motif_gene_rank.csv")
        out_df.to_csv(csv_path, index=False)
        print(f"Wrote {csv_path}")

    print("Done.")


def main():
    parser = argparse.ArgumentParser(
        description="Convert BED-like motif ranking file to feather (optimized v2)."
    )
    parser.add_argument("--input_file", type=str, required=True,
                        help="Name of the input BED file.")
    parser.add_argument("--input_dir", type=str, required=True,
                        help="Directory where the input file is located.")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Directory for output files.")
    parser.add_argument("--chunksize", type=int, default=5_000_000,
                        help="Rows per chunk when reading (default: 5M).")
    parser.add_argument("--fill_value", type=int, default=100000,
                        help="Fill value for missing pairs (default: 100000).")
    parser.add_argument("--also_csv", action="store_true",
                        help="Also write a .csv output.")
    args = parser.parse_args()

    process_files(
        input_file=args.input_file,
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        chunksize=args.chunksize,
        fill_value=args.fill_value,
        also_csv=args.also_csv,
    )


if __name__ == "__main__":
    main()
