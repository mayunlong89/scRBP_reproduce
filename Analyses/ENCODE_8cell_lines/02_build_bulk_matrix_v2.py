#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
build_bulk_matrix_v3.py
=======================
针对 ENCODE RSEM gene quantification 格式:
  列: gene_id | transcript_id(s) | length | effective_length |
       expected_count | TPM | FPKM | posterior_mean_count | ...

处理要点:
  1. gene_id 去版本号: ENSG00000XXX.14 → ENSG00000XXX
  2. 可选: 用 MyGene.info API 把 Ensembl ID 转成 gene symbol
  3. 构建 (replicates × genes) 矩阵
  4. 新增: 对每个 cell line 自动筛选代表性 replicates
     - 参数: --n_replicates
     - 若某个 cell line replicates 不足 N，则全部保留
     - 若超过 N，则保留该 cell line 内平均相关性最高的前 N 个 replicates
"""

import os
import glob
import argparse
from collections import Counter

import numpy as np
import pandas as pd


# ── 目录名映射（与硬盘实际目录一致）─────────────────────
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


# ════════════════════════════════════════════════════
# 核心读取函数
# ════════════════════════════════════════════════════

def read_rsem_gene_tsv(tsv_path, id_type='ensembl'):
    """
    读取 ENCODE RSEM gene quantification TSV。

    Parameters
    ----------
    tsv_path : str
    id_type  : 'ensembl'  → 保留 ENSG ID（去版本号）
               'symbol'   → 转成 gene symbol（需要联网调用 MyGene.info）

    Returns
    -------
    tpm : pd.Series
        index=gene_id, values=TPM (float)
    """
    df = pd.read_csv(tsv_path, sep='\t')

    # ── 找 TPM 列（大小写不敏感）──────────────────────
    tpm_col = next((c for c in df.columns if c.strip().upper() == 'TPM'), None)
    if tpm_col is None:
        raise ValueError(f"No TPM column. Got: {df.columns.tolist()}")

    # ── 找 gene_id 列 ──────────────────────────────────
    gene_col = next((c for c in df.columns if c.strip().lower() == 'gene_id'), None)
    if gene_col is None:
        raise ValueError(f"No gene_id column. Got: {df.columns.tolist()}")

    # ── 去 Ensembl 版本号 ─────────────────────────────
    gene_ids = df[gene_col].astype(str).str.replace(r'\.\d+$', '', regex=True)

    tpm = pd.Series(
        pd.to_numeric(df[tpm_col], errors='coerce').fillna(0).values,
        index=gene_ids
    )

    # 重复ID取最大值（极少）
    tpm = tpm.groupby(level=0).max()

    return tpm


# ════════════════════════════════════════════════════
# 可选：Ensembl ID → Gene Symbol 转换
# ════════════════════════════════════════════════════

def convert_ensembl_to_symbol(ensembl_ids, chunk_size=500):
    """
    用 MyGene.info REST API 批量转换 Ensembl ID → gene symbol。
    不需要安装任何额外包，只需要 requests。

    Parameters
    ----------
    ensembl_ids : list of str
        ENSG格式，已去版本号
    chunk_size : int
        每次API请求的基因数量

    Returns
    -------
    id_map : dict
        {ENSG_id: symbol}，找不到的保留原ENSG ID
    """
    import requests

    id_map = {}
    ids = list(ensembl_ids)

    print(f"  Converting {len(ids)} Ensembl IDs to symbols via MyGene.info ...")
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i + chunk_size]
        payload = {
            'ids': ','.join(chunk),
            'fields': 'symbol',
            'species': 'human',
        }
        try:
            r = requests.post('https://mygene.info/v3/gene', data=payload, timeout=30)
            r.raise_for_status()
            for item in r.json():
                ensg = item.get('query', '')
                sym = item.get('symbol', '')
                id_map[ensg] = sym if sym else ensg
        except Exception as e:
            print(f"  Warning: MyGene.info chunk {i}-{i+chunk_size} failed: {e}")
            for g in chunk:
                id_map[g] = g

    for g in ids:
        if g not in id_map:
            id_map[g] = g

    n_converted = sum(1 for v in id_map.values() if not v.startswith('ENSG'))
    print(f"  Converted {n_converted}/{len(ids)} IDs to symbols")
    return id_map


# ════════════════════════════════════════════════════
# 主函数：构建 bulk 矩阵
# ════════════════════════════════════════════════════

def build_bulk_matrix(
    data_dir='encode_data',
    convert_to_symbol=False,
    min_expressed=1000,
    min_file_kb=50,
):
    """
    Returns
    -------
    bulk_df : pd.DataFrame
        shape = (n_replicates, n_genes)
        index = replicate_id
        columns = gene IDs
    sample_info : pd.DataFrame
        每个 replicate 的元信息
        columns = [replicate_id, cell_line, file_name, file_path, expressed_genes, tpm_sum]
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

    # ── 对齐基因（fillna=0）──────────────────────────
    bulk_df = pd.DataFrame(all_tpm).fillna(0)
    bulk_df = bulk_df.loc[:, (bulk_df > 0).any(axis=0)]   # 去全零列
    bulk_df.index = sample_info['replicate_id'].tolist()

    # ── 可选：转换基因ID ──────────────────────────────
    if convert_to_symbol:
        id_map = convert_ensembl_to_symbol(bulk_df.columns.tolist())
        new_cols = [id_map.get(c, c) for c in bulk_df.columns]
        bulk_df.columns = new_cols
        # 重复symbol合并（取均值）
        bulk_df = bulk_df.T.groupby(level=0).mean().T

    # ── 汇总 ─────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"✅  bulk_df.shape = {bulk_df.shape}")
    print(f"    {bulk_df.shape[0]} replicates × {bulk_df.shape[1]} genes")

    print(f"\nReplicates per cell line:")
    for cl, cnt in sorted(Counter(sample_info['cell_line']).items()):
        warn = "  ⚠️  only 1 rep" if cnt < 2 else ""
        print(f"    {cl:<12}  {cnt}{warn}")

    print(f"\nFirst 5×5 preview:")
    print(bulk_df.iloc[:5, :5].to_string())

    return bulk_df, sample_info


# ════════════════════════════════════════════════════
# 选择代表性 replicates
# ════════════════════════════════════════════════════

def select_representative_replicates(
    bulk_df,
    sample_info,
    n_replicates=None,
    method='mean_corr',
):
    """
    对每个 cell line 选择代表性 replicates。

    Parameters
    ----------
    bulk_df : pd.DataFrame
        index = replicate_id
    sample_info : pd.DataFrame
        必须包含 [replicate_id, cell_line]
    n_replicates : int or None
        每个 cell line 最多保留多少个 replicate。
        若为 None，则不做筛选。
    method : str
        目前支持:
        - 'mean_corr': 按 replicate 在该 cell line 内的平均 Pearson 相关性排序，取前 N

    Returns
    -------
    bulk_df_sel : pd.DataFrame
    sample_info_sel : pd.DataFrame
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
            # 排除对角线自身相关性，避免人为抬高
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
            f"score range among kept: "
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


# ════════════════════════════════════════════════════
# PCA sanity check
# ════════════════════════════════════════════════════

def qc_pca(bulk_df, sample_info, save_path='bulk_pca_qc.pdf'):
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

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
    plt.show()


# ════════════════════════════════════════════════════
# main
# ════════════════════════════════════════════════════

def parse_args():
    parser = argparse.ArgumentParser(
        description="Build bulk TPM matrix from ENCODE RSEM gene quantification TSVs "
                    "and optionally keep top-N representative replicates per cell line."
    )
    parser.add_argument(
        'data_dir',
        nargs='?',
        default='encode_data',
        help='Input directory containing per-cell-line folders (default: encode_data)'
    )
    parser.add_argument(
        '--convert_to_symbol',
        action='store_true',
        help='Convert Ensembl IDs to gene symbols using MyGene.info'
    )
    parser.add_argument(
        '--min_expressed',
        type=int,
        default=1000,
        help='Minimum number of expressed genes (>0 TPM) to keep a file (default: 1000)'
    )
    parser.add_argument(
        '--min_file_kb',
        type=int,
        default=50,
        help='Minimum file size in KB to keep a file (default: 50)'
    )
    parser.add_argument(
        '--n_replicates',
        type=int,
        default=None,
        help='Maximum number of representative replicates to keep per cell line. '
             'If a cell line has fewer replicates, keep all. '
             'If omitted, keep all replicates.'
    )
    parser.add_argument(
        '--prefix',
        default='bulk',
        help='Output prefix (default: bulk)'
    )
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    # ── 步骤1: 构建矩阵 ───────────────────────────────
    bulk_df, sample_info = build_bulk_matrix(
        data_dir=args.data_dir,
        convert_to_symbol=args.convert_to_symbol,
        min_expressed=args.min_expressed,
        min_file_kb=args.min_file_kb,
    )

    if bulk_df is None:
        raise SystemExit(1)

    # ── 步骤2: 筛选代表性 replicates ─────────────────
    bulk_df_sel, sample_info_scored = select_representative_replicates(
        bulk_df=bulk_df,
        sample_info=sample_info,
        n_replicates=args.n_replicates,
        method='mean_corr',
    )

    selected_info = sample_info_scored[sample_info_scored['selected']].copy()

    # ── 步骤3: 保存结果 ─────────────────────────────
       # ── 步骤3: 保存结果 ─────────────────────────────
    matrix_out = f'{args.prefix}_tpm_matrix.csv'
    all_info_out = f'{args.prefix}_replicate_selection.tsv'
    selected_info_out = f'{args.prefix}_selected_replicates.tsv'
    pca_out = f'{args.prefix}_pca_qc.pdf'

    # 保存一个适合 pseudo-cell 模拟的 matrix：
    # 行仍然是 replicate，但 index 改成 cell_line（允许重复）
    bulk_df_to_save = bulk_df_sel.copy()
    bulk_df_to_save.index = selected_info.set_index('replicate_id').loc[bulk_df_sel.index, 'cell_line'].tolist()
    bulk_df_to_save.index.name = 'cell_line'

    bulk_df_to_save.to_csv(matrix_out)
    sample_info_scored.to_csv(all_info_out, sep='\t', index=False)
    selected_info.to_csv(selected_info_out, sep='\t', index=False)

    print(f"\n💾 Saved matrix: {matrix_out}  {bulk_df_to_save.shape}")
    print(f"💾 Saved replicate table (all): {all_info_out}")
    print(f"💾 Saved replicate table (selected only): {selected_info_out}")

    # ── 步骤4: PCA sanity check ─────────────────────
    try:
        qc_pca(bulk_df_sel, selected_info, save_path=pca_out)
    except ImportError:
        print("(sklearn/matplotlib not available, skip PCA plot)")
