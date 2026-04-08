#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ENCODE Cell Line Complete Pipeline
====================================
整合了以下4个步骤：
  Step 1: build_bulk_matrix         — 读取 ENCODE RSEM TSV，构建 bulk TPM 矩阵，筛选代表性 replicates
  Step 2: convert_ensembl_to_symbol — Ensembl ID → Gene Symbol 转换（MyGene.info）
  Step 3: simulate_pseudocells_nb   — Negative Binomial 模拟 pseudo single cells
  Step 4: preprocess_pseudocells    — scRNA-seq 标准预处理，输出 AnnData (.h5ad)

运行方式:
  python encode_cell_line_pipeline.py encode_data --n_replicates 15 --prefix bulk_top5_V2

依赖:
  pip install numpy pandas scikit-learn matplotlib scanpy anndata requests
"""

import os
import glob
import argparse
import time
from collections import Counter

import numpy as np
import pandas as pd
import requests

# ══════════════════════════════════════════════════════════════
# 配置：目录名映射
# ══════════════════════════════════════════════════════════════

CELL_LINES_DIR = {
    'GM12878': 'GM12878',
    'K562':    'K562',
    'HCT116':  'HCT116',
    'MCF-7':   'MCF-7',
    'HepG2':   'HepG2',
    'IMR-90':  'IMR-90',
    'PC3':     'PC3',
    'Panc1':   'Panc1',
}


# ══════════════════════════════════════════════════════════════
# STEP 1: 构建 Bulk TPM 矩阵
# ══════════════════════════════════════════════════════════════

def read_rsem_gene_tsv(tsv_path):
    """读取 ENCODE RSEM gene quantification TSV，返回 TPM Series（已去版本号）"""
    df = pd.read_csv(tsv_path, sep='\t')

    tpm_col = next((c for c in df.columns if c.strip().upper() == 'TPM'), None)
    if tpm_col is None:
        raise ValueError(f"No TPM column. Got: {df.columns.tolist()}")

    gene_col = next((c for c in df.columns if c.strip().lower() == 'gene_id'), None)
    if gene_col is None:
        raise ValueError(f"No gene_id column. Got: {df.columns.tolist()}")

    gene_ids = df[gene_col].astype(str).str.replace(r'\.\d+$', '', regex=True)
    tpm = pd.Series(
        pd.to_numeric(df[tpm_col], errors='coerce').fillna(0).values,
        index=gene_ids
    )
    tpm = tpm.groupby(level=0).max()
    return tpm


def build_bulk_matrix(data_dir='encode_data', min_expressed=1000, min_file_kb=50):
    """
    遍历所有 cell line 子目录，读取 TSV 文件，构建 (replicates × genes) TPM 矩阵。

    Returns: bulk_df, sample_info
    """
    all_tpm = []
    sample_rows = []

    for label, dirname in CELL_LINES_DIR.items():
        pattern = os.path.join(data_dir, dirname, '*.tsv')
        tsvs = sorted(glob.glob(pattern))

        if not tsvs:
            print(f"[{label}] ⚠️  No TSV files at: {pattern}")
            continue

        print(f"\n[{label}]  {len(tsvs)} file(s):")
        for tsv in tsvs:
            kb = os.path.getsize(tsv) // 1024
            file_name = os.path.basename(tsv)
            replicate_id = os.path.splitext(file_name)[0]

            if kb < min_file_kb:
                print(f"  skip ({kb} KB too small): {file_name}")
                continue

            try:
                tpm = read_rsem_gene_tsv(tsv)
                n_expr = int((tpm > 0).sum())
                tpm_sum = float(tpm.sum())

                if n_expr < min_expressed:
                    print(f"  skip ({n_expr} expressed genes < {min_expressed}): {file_name}")
                    continue

                all_tpm.append(tpm)
                sample_rows.append({
                    'replicate_id': replicate_id,
                    'cell_line': label,
                    'file_name': file_name,
                    'file_path': tsv,
                    'expressed_genes': n_expr,
                    'tpm_sum': tpm_sum,
                })
                print(f"  ✓ {file_name:<45} expressed={n_expr:>6,}  TPM_sum={tpm_sum:>10,.0f}")

            except Exception as e:
                print(f"  ✗ {file_name}: {e}")

    if not all_tpm:
        print("\n❌ No data loaded!")
        return None, None

    sample_info = pd.DataFrame(sample_rows)

    bulk_df = pd.DataFrame(all_tpm).fillna(0)
    bulk_df = bulk_df.loc[:, (bulk_df > 0).any(axis=0)]
    bulk_df.index = sample_info['replicate_id'].tolist()

    print(f"\n{'='*60}")
    print(f"✅  bulk_df.shape = {bulk_df.shape}")
    print(f"    {bulk_df.shape[0]} replicates × {bulk_df.shape[1]} genes")
    print(f"\nReplicates per cell line:")
    for cl, cnt in sorted(Counter(sample_info['cell_line']).items()):
        warn = "  ⚠️  only 1 rep" if cnt < 2 else ""
        print(f"    {cl:<12}  {cnt}{warn}")

    return bulk_df, sample_info


def select_representative_replicates(bulk_df, sample_info, n_replicates=None, method='mean_corr'):
    """
    对每个 cell line 保留相关性最高的前 n_replicates 个 replicates。
    若 n_replicates=None 则保留全部。
    """
    if n_replicates is None or n_replicates <= 0:
        print("\nNo replicate subsampling applied.")
        sample_info = sample_info.copy()
        sample_info['selected'] = True
        sample_info['selection_score'] = np.nan
        return bulk_df, sample_info

    keep_ids = []
    out_rows = []

    print(f"\n{'='*60}")
    print(f"Selecting up to {n_replicates} representative replicates per cell line ...")

    for cl in sorted(sample_info['cell_line'].unique()):
        sub_info = sample_info[sample_info['cell_line'] == cl].copy()
        rep_ids = sub_info['replicate_id'].tolist()
        n = len(rep_ids)

        if n <= n_replicates:
            print(f"[{cl}] keep all {n} replicates (<= {n_replicates})")
            sub_info['selected'] = True
            sub_info['selection_score'] = np.nan
            keep_ids.extend(rep_ids)
            out_rows.append(sub_info)
            continue

        sub_expr = bulk_df.loc[rep_ids]
        X = np.log1p(sub_expr.values)

        if method == 'mean_corr':
            corr = np.corrcoef(X)
            if corr.shape[0] > 1:
                np.fill_diagonal(corr, np.nan)
                mean_corr = np.nanmean(corr, axis=1)
            else:
                mean_corr = np.array([1.0], dtype=float)
        else:
            raise ValueError(f"Unsupported method: {method}")

        sub_info['selection_score'] = mean_corr
        sub_info = sub_info.sort_values(
            by=['selection_score', 'expressed_genes', 'tpm_sum', 'replicate_id'],
            ascending=[False, False, False, True]
        )

        selected_ids = sub_info.head(n_replicates)['replicate_id'].tolist()
        sub_info['selected'] = sub_info['replicate_id'].isin(selected_ids)
        keep_ids.extend(selected_ids)
        out_rows.append(sub_info)

        print(
            f"[{cl}] {n} -> {len(selected_ids)} kept | "
            f"score range: "
            f"{sub_info[sub_info['selected']]['selection_score'].min():.4f} - "
            f"{sub_info[sub_info['selected']]['selection_score'].max():.4f}"
        )

    sample_info_scored = pd.concat(out_rows, axis=0, ignore_index=True)
    sample_info_sel = sample_info_scored[sample_info_scored['selected']].copy()
    bulk_df_sel = bulk_df.loc[sample_info_sel['replicate_id']].copy()

    print(f"\nSelected matrix shape: {bulk_df_sel.shape}")
    print("Selected replicates per cell line:")
    for cl, cnt in sorted(Counter(sample_info_sel['cell_line']).items()):
        print(f"    {cl:<12}  {cnt}")

    return bulk_df_sel, sample_info_scored


def qc_pca(bulk_df, sample_info, save_path='bulk_pca_qc.pdf'):
    """PCA sanity check，生成 PDF 图"""
    try:
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
    except ImportError:
        print("(sklearn/matplotlib not available, skip PCA plot)")
        return

    labels = sample_info.set_index('replicate_id').loc[bulk_df.index, 'cell_line'].tolist()

    X = np.log1p(bulk_df.values)
    X = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X)

    uniq = sorted(set(labels))
    cmap = {cl: cm.tab10(i / max(len(uniq), 1)) for i, cl in enumerate(uniq)}

    fig, ax = plt.subplots(figsize=(8, 6))
    for cl in uniq:
        idx = [i for i, l in enumerate(labels) if l == cl]
        ax.scatter(
            coords[idx, 0], coords[idx, 1],
            c=[cmap[cl]], s=120, label=cl,
            edgecolors='white', linewidth=0.5, zorder=3
        )
        cent = coords[idx].mean(axis=0)
        ax.annotate(
            cl, cent, fontsize=9, fontweight='bold',
            color=cmap[cl], ha='center',
            xytext=(0, 10), textcoords='offset points'
        )

    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    ax.set_title('Bulk RNA-seq PCA — sanity check')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\n📊 PCA saved: {save_path}")
    plt.close()


# ══════════════════════════════════════════════════════════════
# STEP 2: Ensembl ID → Gene Symbol 转换
# ══════════════════════════════════════════════════════════════

def ensg_to_symbol(ensg_ids, chunk_size=500, delay=0.3):
    """
    用 MyGene.info REST API 批量将 Ensembl ID 转换为 gene symbol。
    找不到的保留原 ENSG ID。
    """
    id_map = {}
    total = len(ensg_ids)
    ensg_ids = list(ensg_ids)

    for i in range(0, total, chunk_size):
        chunk = ensg_ids[i: i + chunk_size]
        pct = min(i + chunk_size, total)
        print(f"  querying {i}–{pct} / {total} ...", end=' ')

        try:
            r = requests.post(
                'https://mygene.info/v3/gene',
                data={
                    'ids':     ','.join(chunk),
                    'fields':  'symbol,name',
                    'species': 'human',
                },
                timeout=60,
            )
            r.raise_for_status()
            results = r.json()

            n_found = 0
            for item in results:
                ensg   = item.get('query', '')
                symbol = item.get('symbol', '')
                if symbol:
                    id_map[ensg] = symbol
                    n_found += 1
                else:
                    id_map[ensg] = ensg

            print(f"found {n_found}/{len(chunk)}")
            time.sleep(delay)

        except Exception as e:
            print(f"ERROR: {e}")
            for g in chunk:
                id_map[g] = g

    for g in ensg_ids:
        if g not in id_map:
            id_map[g] = g

    return id_map


def convert_ensembl_to_symbol(bulk, prefix):
    """
    将 bulk DataFrame 的列名从 Ensembl ID 转换为 gene symbol，
    保存转换后矩阵和 ID 映射表。
    """
    all_ensg = bulk.columns.tolist()
    ensg_cols  = [c for c in all_ensg if str(c).startswith('ENSG')]
    other_cols = [c for c in all_ensg if not str(c).startswith('ENSG')]

    if other_cols:
        print(f"  ⚠️  {len(other_cols)} non-ENSG columns (keep as-is): {other_cols[:5]}")
    print(f"  ENSG columns to convert: {len(ensg_cols)}")

    print("\nConverting Ensembl IDs → gene symbols via MyGene.info ...")
    id_map = ensg_to_symbol(ensg_cols)

    new_cols = [id_map.get(str(c), str(c)) for c in all_ensg]

    n_converted = sum(1 for c, n in zip(all_ensg, new_cols)
                      if str(c).startswith('ENSG') and not n.startswith('ENSG'))
    n_failed    = sum(1 for c, n in zip(all_ensg, new_cols)
                      if str(c).startswith('ENSG') and n.startswith('ENSG'))
    print(f"\nConversion result:")
    print(f"  ✅ converted to symbol : {n_converted}")
    print(f"  ⚠️  kept as ENSG (no symbol found): {n_failed}")

    bulk.columns = new_cols

    n_before = bulk.shape[1]
    bulk = bulk.T.groupby(level=0).mean().T
    n_after = bulk.shape[1]
    print(f"  Deduplicated: {n_before} → {n_after} columns")

    out_path = f'{prefix}_tpm_matrix_symbol.csv'
    bulk.to_csv(out_path)
    print(f"\n✅ Saved: {out_path}  shape: {bulk.shape}")

    map_df = pd.DataFrame([
        {'ensembl_id': ensg, 'symbol': id_map.get(ensg, ensg)}
        for ensg in ensg_cols
    ])
    map_df.to_csv('ensembl_to_symbol_map.csv', index=False)
    print(f"   ID map saved: ensembl_to_symbol_map.csv")

    return bulk


# ══════════════════════════════════════════════════════════════
# STEP 3: Negative Binomial Pseudo Single Cell 模拟
# ══════════════════════════════════════════════════════════════

def simulate_pseudocells_nb(
        bulk_tpm,
        n_cells_per_line=200,
        target_counts=5000,
        dispersion=0.1,
        seed=42):
    """
    对每个 cell line，从其 replicates 中随机抽取模板，
    用 Negative Binomial 分布模拟 pseudo single cell counts。
    """
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


# ══════════════════════════════════════════════════════════════
# STEP 4: scRNA-seq 标准预处理 → AnnData
# ══════════════════════════════════════════════════════════════

def preprocess_pseudocells(counts_df, labels, prefix):
    """
    完成 scRNA-seq 标准预处理：QC过滤、归一化、log1p、HVG、PCA。
    输出 AnnData (.h5ad) 文件。
    """
    import scanpy as sc
    import anndata as ad

    print("=" * 55)
    print("Step 4-1: Building AnnData")
    print("=" * 55)

    assert len(labels) == counts_df.shape[0], \
        f"Mismatch: {len(labels)} labels vs {counts_df.shape[0]} cells"

    adata = ad.AnnData(
        X   = counts_df.values.astype(np.float32),
        obs = pd.DataFrame(
                  {'cell_line': labels},
                  index=[f'cell_{i}' for i in range(len(labels))]
              ),
        var = pd.DataFrame(index=counts_df.columns)
    )
    print(f"  AnnData created: {adata}")

    print("\nStep 4-2: Preprocessing")

    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.filter_cells(adata, min_genes=200)
    print(f"  filter_genes: {n_genes_before} → {adata.n_vars} genes")
    print(f"  filter_cells: {n_cells_before} → {adata.n_obs} cells")

    adata.layers['counts'] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    print(f"  normalize_total: target_sum=10000")

    sc.pp.log1p(adata)
    print(f"  log1p done")

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes  = 2000,
        flavor       = 'seurat_v3',
        layer        = 'counts',
        span         = 0.3,
    )
    n_hvg = adata.var['highly_variable'].sum()
    print(f"  HVG selected: {n_hvg} genes")

    adata.layers['lognorm'] = adata.X.copy()
    print(f"  lognorm layer saved")

    sc.pp.scale(adata, max_value=10)
    print(f"  scale done (max_value=10)")

    print("\nStep 4-3: PCA (on HVG subset)")
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
    print(f"  PCA done, explained variance ratio (top5): "
          f"{adata.uns['pca']['variance_ratio'][:5].round(3).tolist()}")

    out_h5ad = f'{prefix}_pseudocells_preprocessed.h5ad'
    adata.write_h5ad(out_h5ad)
    print(f"\nSaved: {out_h5ad}")

    print(f"\n{'='*55}")
    print(f"✅  Preprocessing complete")
    print(f"    cells : {adata.n_obs}")
    print(f"    genes : {adata.n_vars}")
    print(f"    HVGs  : {n_hvg}")
    print(f"\n  adata.layers: {list(adata.layers.keys())}")
    print(f"  adata.obsm  : {list(adata.obsm.keys())}")
    print(f"  adata.var columns: {adata.var.columns.tolist()}")
    print(f"\n  Cell line distribution:")
    print(adata.obs['cell_line'].value_counts().sort_index().to_string())
    print("=" * 55)

    return adata


# ══════════════════════════════════════════════════════════════
# 参数解析
# ══════════════════════════════════════════════════════════════

def parse_args():
    parser = argparse.ArgumentParser(
        description="ENCODE Cell Line Complete Pipeline: "
                    "Bulk Matrix → Symbol Conversion → Pseudo Single Cells → Preprocessing"
    )
    parser.add_argument(
        'data_dir',
        nargs='?',
        default='encode_data',
        help='输入目录，包含各 cell line 子文件夹（default: encode_data）'
    )
    parser.add_argument(
        '--n_replicates',
        type=int,
        default=None,
        help='每个 cell line 保留的最多 replicates 数（按相关性筛选；default: 保留全部）'
    )
    parser.add_argument(
        '--prefix',
        default='bulk_top5_V2',
        help='输出文件前缀（default: bulk_top5_V2）'
    )
    parser.add_argument(
        '--min_expressed',
        type=int,
        default=1000,
        help='过滤文件：至少表达基因数（default: 1000）'
    )
    parser.add_argument(
        '--min_file_kb',
        type=int,
        default=50,
        help='过滤文件：最小文件大小 KB（default: 50）'
    )
    parser.add_argument(
        '--n_cells_per_line',
        type=int,
        default=200,
        help='每个 cell line 模拟 pseudo cell 数量（default: 200）'
    )
    parser.add_argument(
        '--target_counts',
        type=int,
        default=5000,
        help='每个 pseudo cell 的目标 UMI counts（default: 5000）'
    )
    parser.add_argument(
        '--dispersion',
        type=float,
        default=0.1,
        help='Negative Binomial dispersion 参数（default: 0.1）'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='随机数种子，保证结果可重现（default: 42）'
    )
    parser.add_argument(
        '--skip_symbol_conversion',
        action='store_true',
        help='跳过 Ensembl ID → Symbol 转换（离线环境使用）'
    )
    return parser.parse_args()


# ══════════════════════════════════════════════════════════════
# 主程序
# ══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    args = parse_args()

    # ─────────────────────────────────────────────────
    # STEP 1: 构建 Bulk TPM 矩阵
    # ─────────────────────────────────────────────────
    print("\n" + "═"*65)
    print("  STEP 1: Building Bulk TPM Matrix")
    print("═"*65)

    bulk_df, sample_info = build_bulk_matrix(
        data_dir=args.data_dir,
        min_expressed=args.min_expressed,
        min_file_kb=args.min_file_kb,
    )

    if bulk_df is None:
        raise SystemExit("❌ No data loaded. Please check your data_dir.")

    # 筛选代表性 replicates
    bulk_df_sel, sample_info_scored = select_representative_replicates(
        bulk_df=bulk_df,
        sample_info=sample_info,
        n_replicates=args.n_replicates,
        method='mean_corr',
    )

    selected_info = sample_info_scored[sample_info_scored['selected']].copy()

    # 保存矩阵（index 改为 cell_line，方便后续步骤）
    bulk_df_to_save = bulk_df_sel.copy()
    bulk_df_to_save.index = selected_info.set_index('replicate_id').loc[
        bulk_df_sel.index, 'cell_line'
    ].tolist()
    bulk_df_to_save.index.name = 'cell_line'

    matrix_out = f'{args.prefix}_tpm_matrix.csv'
    bulk_df_to_save.to_csv(matrix_out)
    sample_info_scored.to_csv(f'{args.prefix}_replicate_selection.tsv', sep='\t', index=False)
    selected_info.to_csv(f'{args.prefix}_selected_replicates.tsv', sep='\t', index=False)

    print(f"\n💾 Saved: {matrix_out}  {bulk_df_to_save.shape}")
    print(f"💾 Saved: {args.prefix}_replicate_selection.tsv")
    print(f"💾 Saved: {args.prefix}_selected_replicates.tsv")

    # PCA sanity check
    try:
        qc_pca(bulk_df_sel, selected_info, save_path=f'{args.prefix}_pca_qc.pdf')
    except Exception as e:
        print(f"(PCA skipped: {e})")

    # ─────────────────────────────────────────────────
    # STEP 2: Ensembl ID → Gene Symbol 转换
    # ─────────────────────────────────────────────────
    print("\n" + "═"*65)
    print("  STEP 2: Converting Ensembl ID → Gene Symbol")
    print("═"*65)

    if args.skip_symbol_conversion:
        print("  (skipped by --skip_symbol_conversion)")
        bulk_symbol = bulk_df_to_save.copy()
        symbol_matrix_path = matrix_out
    else:
        bulk_symbol = convert_ensembl_to_symbol(bulk_df_to_save.copy(), args.prefix)
        symbol_matrix_path = f'{args.prefix}_tpm_matrix_symbol.csv'

    # ─────────────────────────────────────────────────
    # STEP 3: Negative Binomial Pseudo Single Cell 模拟
    # ─────────────────────────────────────────────────
    print("\n" + "═"*65)
    print("  STEP 3: Simulating Pseudo Single Cells (Negative Binomial)")
    print("═"*65)

    print(f"\nbulk_tpm loaded: {bulk_symbol.shape}")
    print(f"cell lines: {bulk_symbol.index.unique().tolist()}")
    print(f"\nSimulating pseudo-cells...")

    counts_df, labels = simulate_pseudocells_nb(
        bulk_symbol,
        n_cells_per_line=args.n_cells_per_line,
        target_counts=args.target_counts,
        dispersion=args.dispersion,
        seed=args.seed,
    )

    print(f"\ncounts_df.shape = {counts_df.shape}")
    print(f"labels sample   = {labels[:4]}")

    counts_out = f'{args.prefix}_pseudocell_counts_{args.n_cells_per_line}cells.csv'
    labels_out = f'{args.prefix}_pseudocell_labels_{args.n_cells_per_line}cells.csv'

    counts_df.to_csv(counts_out)
    pd.Series(labels, name='cell_line').to_csv(labels_out, index=False)
    print(f"\n💾 Saved: {counts_out}")
    print(f"💾 Saved: {labels_out}")

    # ─────────────────────────────────────────────────
    # STEP 4: scRNA-seq 标准预处理 → AnnData
    # ─────────────────────────────────────────────────
    print("\n" + "═"*65)
    print("  STEP 4: Preprocessing Pseudo Cells → AnnData")
    print("═"*65)

    try:
        adata = preprocess_pseudocells(counts_df, labels, args.prefix)
    except ImportError as e:
        print(f"\n⚠️  scanpy/anndata not available, Step 4 skipped: {e}")
        print("   Install with: pip install scanpy anndata")

    print("\n" + "═"*65)
    print("  ✅  ALL STEPS COMPLETE")
    print("═"*65)
