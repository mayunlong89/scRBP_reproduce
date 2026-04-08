"""
build_bulk_matrix_v2.py
=======================
针对 ENCODE RSEM gene quantification 格式:
  列: gene_id | transcript_id(s) | length | effective_length |
       expected_count | TPM | FPKM | posterior_mean_count | ...

处理要点:
  1. gene_id 去版本号: ENSG00000XXX.14 → ENSG00000XXX
  2. 可选: 用 MyGene.info API 把 Ensembl ID 转成 gene symbol
  3. 构建 (replicates × genes) 矩阵
"""

import pandas as pd
import numpy as np
import glob, os, re
from collections import Counter

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
    tpm : pd.Series  index=gene_id, values=TPM (float)
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
    # ENSG00000223972.5  →  ENSG00000223972
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
    ensembl_ids : list of str (ENSG格式，已去版本号)
    chunk_size  : 每次API请求的基因数量

    Returns
    -------
    id_map : dict  {ENSG_id: symbol}  找不到的保留原ENSG ID
    """
    import requests
    id_map = {}
    ids = list(ensembl_ids)

    print(f"  Converting {len(ids)} Ensembl IDs to symbols via MyGene.info ...")
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i + chunk_size]
        payload = {
            'ids':    ','.join(chunk),
            'fields': 'symbol',
            'species': 'human',
        }
        try:
            r = requests.post('https://mygene.info/v3/gene',
                              data=payload, timeout=30)
            r.raise_for_status()
            for item in r.json():
                ensg = item.get('query', '')
                sym  = item.get('symbol', '')
                id_map[ensg] = sym if sym else ensg   # fallback to ENSG if no symbol
        except Exception as e:
            print(f"  Warning: MyGene.info chunk {i}–{i+chunk_size} failed: {e}")
            for g in chunk:
                id_map[g] = g   # fallback

    # 未找到的保留ENSG
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
        data_dir          = 'encode_data',
        convert_to_symbol = False,   # True = 转成 gene symbol（需联网）
        min_expressed     = 1000,    # 过滤：至少有这么多基因表达 > 0
        min_file_kb       = 50,      # 过滤：文件太小跳过
):
    """
    Returns
    -------
    bulk_df : DataFrame  (n_replicates × n_genes)
              index   = cell_line 标签
              columns = ENSG ID (或 gene symbol if convert_to_symbol=True)
    labels  : list[str]
    """
    all_tpm, labels = [], []

    for label, dirname in CELL_LINES_DIR.items():
        pattern = os.path.join(data_dir, dirname, '*.tsv')
        tsvs    = sorted(glob.glob(pattern))

        if not tsvs:
            print(f"[{label}] ⚠️  No TSV files at: {pattern}")
            continue

        print(f"\n[{label}]  {len(tsvs)} file(s):")
        for tsv in tsvs:
            kb = os.path.getsize(tsv) // 1024
            if kb < min_file_kb:
                print(f"  skip ({kb} KB too small): {os.path.basename(tsv)}")
                continue
            try:
                tpm      = read_rsem_gene_tsv(tsv)
                n_expr   = (tpm > 0).sum()
                if n_expr < min_expressed:
                    print(f"  skip ({n_expr} expressed genes < {min_expressed}): "
                          f"{os.path.basename(tsv)}")
                    continue
                all_tpm.append(tpm)
                labels.append(label)
                print(f"  ✓ {os.path.basename(tsv):<45} "
                      f"expressed={n_expr:>6,}  TPM_sum={tpm.sum():>10,.0f}")
            except Exception as e:
                print(f"  ✗ {os.path.basename(tsv)}: {e}")

    if not all_tpm:
        print("\n❌ No data loaded!")
        return None, None

    # ── 对齐基因（fillna=0）──────────────────────────
    bulk_df = pd.DataFrame(all_tpm).fillna(0)
    bulk_df = bulk_df.loc[:, (bulk_df > 0).any(axis=0)]   # 去全零列
    bulk_df.index = pd.Index(labels, name='cell_line')

    # ── 可选：转换基因ID ──────────────────────────────
    if convert_to_symbol:
        id_map  = convert_ensembl_to_symbol(bulk_df.columns.tolist())
        new_cols = [id_map.get(c, c) for c in bulk_df.columns]
        bulk_df.columns = new_cols
        # 重复symbol合并（取均值）
        bulk_df = bulk_df.T.groupby(level=0).mean().T

    # ── 汇总 ─────────────────────────────────────────
    print(f"\n{'='*55}")
    print(f"✅  bulk_df.shape = {bulk_df.shape}")
    print(f"    {bulk_df.shape[0]} replicates  ×  {bulk_df.shape[1]} genes")
    print(f"\nReplicates per cell line:")
    for cl, cnt in sorted(Counter(labels).items()):
        warn = "  ⚠️  only 1 rep" if cnt < 2 else ""
        print(f"    {cl:<12}  {cnt}{warn}")

    # ── 快速 sanity：显示前5行5列 ────────────────────
    print(f"\nFirst 5×5 preview:")
    print(bulk_df.iloc[:5, :5].to_string())

    return bulk_df, labels


# ════════════════════════════════════════════════════
# Sanity check PCA
# ════════════════════════════════════════════════════

def qc_pca(bulk_df, labels, save_path='bulk_pca_qc.pdf'):
    from sklearn.decomposition import PCA
    from sklearn.preprocessing  import StandardScaler
    import matplotlib.pyplot as plt
    import matplotlib.cm       as cm

    X      = np.log1p(bulk_df.values)
    X      = StandardScaler().fit_transform(X)
    coords = PCA(n_components=2, random_state=42).fit_transform(X)

    uniq   = sorted(set(labels))
    cmap   = {cl: cm.tab10(i / len(uniq)) for i, cl in enumerate(uniq)}

    fig, ax = plt.subplots(figsize=(8, 6))
    for cl in uniq:
        idx = [i for i, l in enumerate(labels) if l == cl]
        ax.scatter(coords[idx, 0], coords[idx, 1],
                   c=[cmap[cl]], s=120, label=cl,
                   edgecolors='white', linewidth=0.5, zorder=3)
        cent = coords[idx].mean(axis=0)
        ax.annotate(cl, cent, fontsize=9, fontweight='bold',
                    color=cmap[cl], ha='center',
                    xytext=(0, 10), textcoords='offset points')

    ax.set_xlabel('PC1'); ax.set_ylabel('PC2')
    ax.set_title('Bulk RNA-seq PCA — sanity check\n(8 cell lines should separate clearly)')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\n📊 PCA saved: {save_path}")
    plt.show()


# ════════════════════════════════════════════════════
if __name__ == '__main__':
    import sys
    data_dir = sys.argv[1] if len(sys.argv) > 1 else 'encode_data'

    # ── 步骤1: 构建矩阵（保留 Ensembl ID）────────────
    bulk_df, labels = build_bulk_matrix(
        data_dir          = data_dir,
        convert_to_symbol = False,   # 先保持ENSG ID，后面用gene list匹配时再转
        min_expressed     = 1000,
        min_file_kb       = 50,
    )

    if bulk_df is None:
        sys.exit(1)

    # ── 步骤2: 保存 ───────────────────────────────────
    bulk_df.to_csv('bulk_tpm_matrix.csv')
    print(f"\n💾 Saved: bulk_tpm_matrix.csv  {bulk_df.shape}")

    # ── 步骤3: PCA sanity check ──────────────────────
    try:
        qc_pca(bulk_df, labels)
    except ImportError:
        print("(sklearn/matplotlib not available, skip PCA plot)")
