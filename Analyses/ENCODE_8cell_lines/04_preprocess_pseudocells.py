"""
preprocess_pseudocells.py
=========================
读取 pseudocell_count_200cells.csv + pseudocell_labels_gene_symbol_200cells.csv
完成 scRNA-seq 标准预处理，输出 AnnData 对象
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

# ══════════════════════════════════════════════════
# 1. 读取数据
# ══════════════════════════════════════════════════
print("=" * 55)
print("Step 1: Loading data")
print("=" * 55)

counts_df = pd.read_csv('pseudocell_counts_200cells.csv', index_col=0)
labels_df = pd.read_csv('pseudocell_labels_200cells.csv')

# 兼容labels文件的各种列名
label_col = labels_df.columns[0]
labels    = labels_df[label_col].tolist()

print(f"  counts_df : {counts_df.shape}  (cells × genes)")
print(f"  labels    : {len(labels)} cells, {len(set(labels))} cell lines")
print(f"  cell lines: {sorted(set(labels))}")

assert len(labels) == counts_df.shape[0], \
    f"Mismatch: {len(labels)} labels vs {counts_df.shape[0]} cells"

# ══════════════════════════════════════════════════
# 2. 构建 AnnData
# ══════════════════════════════════════════════════
print("\nStep 2: Building AnnData")

adata = ad.AnnData(
    X   = counts_df.values.astype(np.float32),
    obs = pd.DataFrame(
              {'cell_line': labels},
              index=[f'cell_{i}' for i in range(len(labels))]
          ),
    var = pd.DataFrame(index=counts_df.columns)
)
print(f"  AnnData created: {adata}")

# ══════════════════════════════════════════════════
# 3. 预处理
# ══════════════════════════════════════════════════
print("\nStep 3: Preprocessing")

# 3-1. 基本质控过滤
n_cells_before = adata.n_obs
n_genes_before = adata.n_vars
sc.pp.filter_genes(adata, min_cells=10)       # 至少10个细胞中表达
sc.pp.filter_cells(adata, min_genes=200)      # 细胞至少200个基因表达
print(f"  filter_genes: {n_genes_before} → {adata.n_vars} genes")
print(f"  filter_cells: {n_cells_before} → {adata.n_obs} cells")

# 3-2. 保存原始counts（用于seurat_v3 HVG）
adata.layers['counts'] = adata.X.copy()

# 3-3. 文库大小归一化
sc.pp.normalize_total(adata, target_sum=1e4)
print(f"  normalize_total: target_sum=10000")

# 3-4. log1p 变换
sc.pp.log1p(adata)
print(f"  log1p done")

# 3-5. 高变基因（HVG baseline arm）
#      seurat_v3 需要原始counts，从layer读取
sc.pp.highly_variable_genes(
    adata,
    n_top_genes  = 2000,
    flavor       = 'seurat_v3',
    layer        = 'counts',   # ← 用raw counts计算variance
    span         = 0.3,
)
n_hvg = adata.var['highly_variable'].sum()
print(f"  HVG selected: {n_hvg} genes")

# 3-6. 保存 log-normalized 层（用于 RBP/TF regulon打分）
adata.layers['lognorm'] = adata.X.copy()
print(f"  lognorm layer saved")

# 3-7. Scale（用于PCA，不影响lognorm层）
sc.pp.scale(adata, max_value=10)
print(f"  scale done (max_value=10)")

# ══════════════════════════════════════════════════
# 4. 快速 PCA + UMAP 可视化（可选sanity check）
# ══════════════════════════════════════════════════
print("\nStep 4: PCA (on HVG subset)")
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
print(f"  PCA done, explained variance ratio (top5): "
      f"{adata.uns['pca']['variance_ratio'][:5].round(3).tolist()}")

# ══════════════════════════════════════════════════
# 5. 保存
# ══════════════════════════════════════════════════
print("\nStep 5: Saving")
adata.write_h5ad('pseudocells_preprocessed_200cells.h5ad')
print(f"  Saved: pseudocells_preprocessed_200cells.h5ad")

# ══════════════════════════════════════════════════
# 6. 汇总报告
# ══════════════════════════════════════════════════
print(f"\n{'='*55}")
print(f"✅  Preprocessing complete")
print(f"    cells : {adata.n_obs}")
print(f"    genes : {adata.n_vars}")
print(f"    HVGs  : {n_hvg}")
print(f"\n  adata.layers available:")
for k in adata.layers:
    print(f"    '{k}' — shape {adata.layers[k].shape}")
print(f"\n  adata.obsm:")
for k in adata.obsm:
    print(f"    '{k}' — shape {adata.obsm[k].shape}")
print(f"\n  adata.var columns: {adata.var.columns.tolist()}")
print(f"\n  Cell line distribution:")
print(adata.obs['cell_line'].value_counts().sort_index().to_string())
print("="*55)
