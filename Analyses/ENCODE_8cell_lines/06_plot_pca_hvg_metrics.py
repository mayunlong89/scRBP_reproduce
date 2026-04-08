#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_pca_hvg_metrics.py
=======================
基于 h5ad 文件中的 top HVGs 进行:
1) PCA
2) KMeans clustering
3) ARI / NMI evaluation
4) PCA scatter plot

输入:
- h5ad 文件
- obs 中必须包含: cell_line
- var 中必须包含: highly_variable

输出:
- PCA 图 (PDF / PNG)
- metrics.txt
- metrics.csv
"""

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.cluster import KMeans


CELL_LINE_COLORS = {
    'GM12878': '#E41A1C',   # red
    'K562':    '#FF7F00',   # orange
    'HCT116':  '#4DAF4A',   # green
    'MCF-7':   '#F781BF',   # pink
    'MCF7':    '#F781BF',
    'HepG2':   '#377EB8',   # blue
    'IMR-90':  '#984EA3',   # purple
    'IMR90':   '#984EA3',
    'PC3':     '#FF69B4',   # hot pink
    'Panc1':   '#8B4513',   # brown
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run PCA on HVGs from h5ad and compute ARI/NMI"
    )
    parser.add_argument(
        "--h5ad",
        required=True,
        help="Input h5ad file"
    )
    parser.add_argument(
        "--outprefix",
        required=True,
        help="Output prefix"
    )
    parser.add_argument(
        "--n_comps",
        type=int,
        default=50,
        help="Number of PCs to compute (default: 50)"
    )
    parser.add_argument(
        "--random_state",
        type=int,
        default=42,
        help="Random seed for KMeans (default: 42)"
    )
    parser.add_argument(
        "--kmeans_n_init",
        type=int,
        default=50,
        help="n_init for KMeans (default: 50)"
    )
    parser.add_argument(
        "--title",
        default="Top Variable genes",
        help="Plot title"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # ══════════════════════════════════════════════════
    # 1. 读取数据
    # ══════════════════════════════════════════════════
    print(f"\n[1/4] Loading h5ad: {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    print(f"Loaded: {adata}")

    if "cell_line" not in adata.obs.columns:
        raise ValueError("adata.obs must contain column 'cell_line'.")

    if "highly_variable" not in adata.var.columns:
        raise ValueError("adata.var must contain column 'highly_variable'.")

    # ══════════════════════════════════════════════════
    # 2. 提取 HVG 子集 + PCA
    # ══════════════════════════════════════════════════
    print("\n[2/4] Running PCA on HVGs ...")
    hvg_mask = adata.var["highly_variable"].astype(bool).values
    n_hvg = int(hvg_mask.sum())

    if n_hvg == 0:
        raise ValueError("No highly variable genes found in adata.var['highly_variable'].")

    adata_hvg = adata[:, hvg_mask].copy()
    print(f"  HVG subset: {adata_hvg.shape}")

    sc.tl.pca(adata_hvg, n_comps=args.n_comps, svd_solver='arpack')

    # ══════════════════════════════════════════════════
    # 3. KMeans + ARI + NMI
    # ══════════════════════════════════════════════════
    print("\n[3/4] Running KMeans + metrics ...")
    labels = adata_hvg.obs["cell_line"].astype(str).values
    n_clusters = len(np.unique(labels))

    km = KMeans(
        n_clusters=n_clusters,
        n_init=args.kmeans_n_init,
        random_state=args.random_state
    )
    pred = km.fit_predict(adata_hvg.obsm["X_pca"])

    ari = adjusted_rand_score(labels, pred)
    nmi = normalized_mutual_info_score(labels, pred)

    print(f"  Number of cell lines = {n_clusters}")
    print(f"  ARI = {ari:.4f}")
    print(f"  NMI = {nmi:.4f}")

    coords = adata_hvg.obsm["X_pca"][:, :2]

    # 保存 metrics
    metrics_txt = f"{args.outprefix}.metrics.txt"
    with open(metrics_txt, "w") as f:
        f.write(f"ARI\t{ari:.6f}\n")
        f.write(f"NMI\t{nmi:.6f}\n")
    print(f"Saved metrics: {metrics_txt}")

    metrics_csv = f"{args.outprefix}.metrics.csv"
    pd.DataFrame([{
        "method": "HVG",
        "input_h5ad": args.h5ad,
        "n_cells": adata_hvg.n_obs,
        "n_features": adata_hvg.n_vars,
        "ARI": ari,
        "NMI": nmi,
        "random_state": args.random_state,
        "kmeans_n_init": args.kmeans_n_init
    }]).to_csv(metrics_csv, index=False)
    print(f"Saved metrics csv: {metrics_csv}")

    # 保存 PCA 坐标
    coords_csv = f"{args.outprefix}.pca_coordinates.csv"
    pd.DataFrame({
        "cell": adata_hvg.obs_names.astype(str),
        "cell_line": labels,
        "kmeans_cluster": pred,
        "PC1": coords[:, 0],
        "PC2": coords[:, 1]
    }).to_csv(coords_csv, index=False)
    print(f"Saved PCA coordinates: {coords_csv}")

    # ══════════════════════════════════════════════════
    # 4. 画图
    # ══════════════════════════════════════════════════
    print("\n[4/4] Plotting PCA ...")
    fig, ax = plt.subplots(figsize=(6, 5))

    cell_lines = sorted(np.unique(labels))

    for cl in cell_lines:
        mask = labels == cl
        color = CELL_LINE_COLORS.get(cl, '#888888')

        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            c=color,
            s=6,
            alpha=0.6,
            linewidths=0,
            rasterized=True,
        )

        cx, cy = coords[mask, 0].mean(), coords[mask, 1].mean()
        ax.text(
            cx, cy, cl,
            fontsize=9,
            fontweight='bold',
            color=color,
            ha='center',
            va='center',
            path_effects=[pe.withStroke(linewidth=2.5, foreground='white')],
        )

    metric_text = f'ARI: {ari:.2f}\nNMI: {nmi:.2f}'
    ax.text(
        0.03, 0.97, metric_text,
        transform=ax.transAxes,
        fontsize=11,
        fontweight='bold',
        va='top',
        ha='left',
        color='black',
    )

    ax.set_title(args.title, fontsize=12, fontweight='bold', pad=8)

    var_ratio = adata_hvg.uns['pca']['variance_ratio']
    ax.set_xlabel(f'PC1 ({var_ratio[0]*100:.1f}%)', fontsize=10)
    ax.set_ylabel(f'PC2 ({var_ratio[1]*100:.1f}%)', fontsize=10)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.2)
    ax.spines['bottom'].set_linewidth(1.2)
    ax.tick_params(labelsize=8)

    plt.tight_layout()

    pdf_file = f"{args.outprefix}.pdf"
    png_file = f"{args.outprefix}.png"
    plt.savefig(pdf_file, dpi=300, bbox_inches='tight')
    plt.savefig(png_file, dpi=200, bbox_inches='tight')
    plt.close()

    print(f"Saved: {pdf_file}")
    print(f"Saved: {png_file}")
    print("\n✅ Done.")


if __name__ == "__main__":
    main()
