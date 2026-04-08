import numpy as np
import pandas as pd

# ── 1. 读取 bulk TPM 矩阵 ─────────────────────────────────────────────────
bulk_tpm = pd.read_csv('bulk_top5_V2_tpm_matrix_symbol.csv', index_col=0)
print(f"bulk_tpm loaded: {bulk_tpm.shape}")          # (n_replicates, n_genes)
print(f"cell lines: {bulk_tpm.index.unique().tolist()}")

# ── 2. 模拟函数 ───────────────────────────────────────────────────────────
def simulate_pseudocells_nb(
        bulk_tpm,
        n_cells_per_line = 200,
        target_counts    = 5000,
        dispersion       = 0.1,
        seed             = 42):

    np.random.seed(seed)
    pseudo_cells, labels = [], []

    for cell_line in bulk_tpm.index.unique():
        replicates = bulk_tpm.loc[cell_line]
        if replicates.ndim == 1:
            replicates = replicates.to_frame().T

        print(f"  {cell_line}: {len(replicates)} replicate(s) → {n_cells_per_line} pseudo-cells")

        for _ in range(n_cells_per_line):
            template  = replicates.sample(1).values.flatten()
            template  = np.maximum(template, 0)
            total_tpm = template.sum()
            if total_tpm == 0:
                continue

            mu      = template / total_tpm * target_counts
            n_param = 1.0 / dispersion
            p_param = np.clip(n_param / (n_param + mu + 1e-9), 1e-9, 1 - 1e-9)
            counts  = np.random.negative_binomial(n_param, p_param)

            pseudo_cells.append(counts)
            labels.append(cell_line)

    counts_df = pd.DataFrame(
        np.array(pseudo_cells, dtype=np.int32),
        columns=bulk_tpm.columns
    )
    return counts_df, labels

# ── 3. 运行模拟 ────────────────────────────────────────────────────────────
print("\nSimulating pseudo-cells...")
counts_df, labels = simulate_pseudocells_nb(bulk_tpm)

print(f"\ncounts_df.shape = {counts_df.shape}")   # (1600, ~58000)
print(f"labels sample   = {labels[:4]}")

# ── 4. 保存 ───────────────────────────────────────────────────────────────
counts_df.to_csv('pseudocell_counts_200cells.csv')
pd.Series(labels, name='cell_line').to_csv('pseudocell_labels_200cells.csv', index=False)
print("\nSaved: pseudocell_counts_200cells.csv, pseudocell_labels_200cells.csv")
