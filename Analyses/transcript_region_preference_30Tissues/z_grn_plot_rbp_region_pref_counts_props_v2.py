#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram

REGIONS = ["3UTR", "5UTR", "CDS", "Introns"]

REGION_COLORS = {
    "3UTR": "#0B6B3A",
    "5UTR": "#A6D854",
    "CDS": "#2CB34A",
    "Introns": "#FFB000",
}

CLUSTER_COLORS = {
  1: "#4E79A7",  # blue
  2: "#E15759",  # red
  3: "#B07AA1",  # purple
  4: "#9C755F",  # brown
}



def _require_cols(df: pd.DataFrame, cols: list, context: str):
    miss = [c for c in cols if c not in df.columns]
    if miss:
        raise SystemExit(f"[ERROR] Missing columns in {context}: {miss}")


def _save_fig(fig, out_base):
    fig.savefig(out_base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(out_base + ".pdf", bbox_inches="tight")
    plt.close(fig)


def add_dominant_region(df: pd.DataFrame) -> pd.DataFrame:
    prop_cols = [f"prop_{r}" for r in REGIONS]
    _require_cols(df, prop_cols, "stats_csv (need prop_*)")
    out = df.copy()
    out["dominant_region"] = out[prop_cols].idxmax(axis=1).str.replace("prop_", "", regex=False)
    return out


def sort_by_preference(df: pd.DataFrame, metric: str = "prop") -> pd.DataFrame:
    if metric not in ("prop", "count"):
        raise ValueError("metric must be 'prop' or 'count'")

    df = df.copy()
    count_cols = [f"nTissues_{r}" for r in REGIONS]
    _require_cols(df, count_cols, "stats_csv (need nTissues_*)")

    if "total_region_occurrences" not in df.columns:
        df["total_region_occurrences"] = df[count_cols].sum(axis=1)

    parts = []
    for r in REGIONS:
        sub = df[df["dominant_region"] == r].copy()
        if sub.empty:
            continue

        if metric == "prop":
            sub = sub.sort_values([f"prop_{r}", "total_region_occurrences"], ascending=[False, False])
        else:
            sub = sub.sort_values([f"nTissues_{r}", "total_region_occurrences"], ascending=[False, False])

        parts.append(sub)

    return pd.concat(parts, ignore_index=True) if parts else df.iloc[0:0].copy()


def select_rbplist(df_sorted: pd.DataFrame, mode: str, top_n: int, n_per_region: int) -> pd.DataFrame:
    mode = mode.lower()
    if mode == "all":
        return df_sorted.copy()
    if mode == "top":
        return df_sorted.head(top_n).copy()
    if mode == "balanced":
        parts = []
        for r in REGIONS:
            parts.append(df_sorted[df_sorted["dominant_region"] == r].head(n_per_region).copy())
        return pd.concat(parts, ignore_index=True)
    raise ValueError("mode must be one of: all, top, balanced")


def hierarchical_cluster(df: pd.DataFrame, n_clusters: int = 4, method: str = "ward"):
    prop_cols = [f"prop_{r}" for r in REGIONS]
    _require_cols(df, prop_cols, "hierarchical_cluster")

    x = df[prop_cols].values.astype(float)
    Z = linkage(x, method=method, metric="euclidean")

    # dendrogram order
    dendro = dendrogram(Z, no_plot=True)
    leaves = dendro["leaves"]

    # cluster assignment
    raw_clusters = fcluster(Z, t=n_clusters, criterion="maxclust")

    out = df.copy()
    out["cluster_raw"] = raw_clusters

    # reorder cluster labels by mean profile for interpretability
    cluster_info = []
    for cl in sorted(out["cluster_raw"].unique()):
        sub = out[out["cluster_raw"] == cl]
        means = sub[prop_cols].mean()
        dom_region = means.idxmax().replace("prop_", "")
        dom_value = means.max()
        cluster_info.append((cl, dom_region, dom_value))

    region_rank = {r: i for i, r in enumerate(REGIONS)}
    cluster_info = sorted(cluster_info, key=lambda x: (region_rank[x[1]], -x[2]))
    relabel = {old: i + 1 for i, (old, _, _) in enumerate(cluster_info)}

    out["cluster"] = out["cluster_raw"].map(relabel)

    # reorder rows by dendrogram
    out = out.iloc[leaves].copy().reset_index(drop=True)
    out["plot_order"] = np.arange(out.shape[0])

    return out, Z


def summarize_clusters(df: pd.DataFrame) -> pd.DataFrame:
    prop_cols = [f"prop_{r}" for r in REGIONS]
    rows = []
    for cl in sorted(df["cluster"].unique()):
        sub = df[df["cluster"] == cl]
        means = sub[prop_cols].mean()
        rows.append({
            "cluster": cl,
            "n_rbps": sub.shape[0],
            "mean_prop_3UTR": means["prop_3UTR"],
            "mean_prop_5UTR": means["prop_5UTR"],
            "mean_prop_CDS": means["prop_CDS"],
            "mean_prop_Introns": means["prop_Introns"],
            "dominant_region": means.idxmax().replace("prop_", "")
        })
    return pd.DataFrame(rows)


def plot_horizontal_clustered(df: pd.DataFrame, Z, out_base: str, title: str = None):
    """
    Layout:
      dendrogram | rbp names + stacked bars | cluster strip
    """
    prop_cols = [f"prop_{r}" for r in REGIONS]
    _require_cols(df, ["rbp", "cluster"] + prop_cols, "plot_horizontal_clustered")

    n = df.shape[0]
    fig_h = max(8, 0.22 * n)
    fig = plt.figure(figsize=(10, fig_h))

    gs = fig.add_gridspec(
        nrows=1, ncols=3,
        width_ratios=[1.7, 5.8, 0.4],
        wspace=0.05
    )

    ax_d = fig.add_subplot(gs[0, 0])   # dendrogram
    ax_b = fig.add_subplot(gs[0, 1])   # bars + names
    ax_c = fig.add_subplot(gs[0, 2])   # cluster strip

    # ----- left dendrogram -----
    dendrogram(
        Z,
        orientation="left",
        no_labels=True,
        color_threshold=None,
        above_threshold_color="black",
        ax=ax_d
    )
    ax_d.invert_yaxis()
    ax_d.set_xticks([])
    ax_d.set_yticks([])
    for s in ["top", "right", "bottom", "left"]:
        ax_d.spines[s].set_visible(False)

    # ----- middle stacked horizontal bars -----
    y = np.arange(n)
    left = np.zeros(n)

    for r in REGIONS:
        vals = df[f"prop_{r}"].values
        ax_b.barh(
            y, vals, left=left,
            color=REGION_COLORS[r],
            edgecolor="none",
            height=0.85,
            label=r
        )
        left += vals

    ax_b.set_ylim(-0.5, n - 0.5)
    ax_b.invert_yaxis()
    ax_b.set_xlim(0, 1.0)
    ax_b.set_xlabel("Proportion")
    ax_b.set_yticks(y)
    ax_b.set_yticklabels(df["rbp"].tolist(), fontsize=7)

    if title is not None:
        ax_b.set_title(title, fontsize=13, pad=10)

    # region legend
    handles, labels = ax_b.get_legend_handles_labels()
    label_map = {"3UTR": "3′UTR", "5UTR": "5′UTR", "CDS": "CDS", "Introns": "Introns"}
    labels = [label_map.get(x, x) for x in labels]
    ax_b.legend(handles, labels, title="Region", frameon=False,
                bbox_to_anchor=(1.02, 1.00), loc="upper left")

    # light horizontal separators
    for yy in np.arange(-0.5, n, 1):
        ax_b.axhline(yy, color="white", lw=0.35, zorder=0)

    # ----- right cluster strip -----
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(-0.5, n - 0.5)
    ax_c.invert_yaxis()

    for i, cl in enumerate(df["cluster"].tolist()):
        ax_c.add_patch(Rectangle((0, i - 0.425), 1, 0.85,
                                 facecolor=CLUSTER_COLORS.get(cl, "#999999"),
                                 edgecolor="none"))

    # boundary lines between clusters
    cluster_vals = df["cluster"].tolist()
    for i in range(1, n):
        if cluster_vals[i] != cluster_vals[i - 1]:
            ax_b.axhline(i - 0.5, color="black", lw=1.0)
            ax_c.axhline(i - 0.5, color="white", lw=1.2)

    # add cluster labels on strip center for each block
    start = 0
    for i in range(1, n + 1):
        if i == n or cluster_vals[i] != cluster_vals[i - 1]:
            cl = cluster_vals[start]
            mid = (start + i - 1) / 2
            ax_c.text(0.5, mid, f"C{cl}", ha="center", va="center",
                      rotation=90, fontsize=9, fontweight="bold", color="white")
            start = i

    ax_c.set_xticks([])
    ax_c.set_yticks([])
    ax_c.set_ylabel("Cluster", rotation=270, labelpad=15)
    for s in ["top", "right", "bottom", "left"]:
        ax_c.spines[s].set_visible(False)

    plt.tight_layout()
    _save_fig(fig, out_base)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats_csv", required=True)
    ap.add_argument("--prefix", default="scRBP_RBPRegionPref")
    ap.add_argument("--min_tissues", type=int, default=5)

    ap.add_argument("--mode", choices=["all", "top", "balanced"], default="balanced")
    ap.add_argument("--top_n", type=int, default=200)
    ap.add_argument("--n_per_region", type=int, default=50)

    ap.add_argument("--n_clusters", type=int, default=4)
    ap.add_argument("--cluster_method", choices=["ward", "average", "complete", "single"], default="ward")

    args = ap.parse_args()

    df = pd.read_csv(args.stats_csv)

    count_cols = [f"nTissues_{r}" for r in REGIONS]
    prop_cols = [f"prop_{r}" for r in REGIONS]

    _require_cols(df, ["rbp", "nTissues_anyRegion"] + count_cols + prop_cols, "stats_csv")

    df = df[df["nTissues_anyRegion"] >= args.min_tissues].copy()
    if df.empty:
        raise SystemExit(f"[ERROR] No RBPs left after filtering nTissues_anyRegion >= {args.min_tissues}")

    df = add_dominant_region(df)

    if "total_region_occurrences" not in df.columns:
        df["total_region_occurrences"] = df[count_cols].sum(axis=1)

    # choose by prop preference, because you want proportion-based clusters
    df_sorted = sort_by_preference(df, metric="prop")
    df_sel = select_rbplist(df_sorted, args.mode, args.top_n, args.n_per_region)

    # hierarchical clustering
    df_clu, Z = hierarchical_cluster(
        df_sel,
        n_clusters=args.n_clusters,
        method=args.cluster_method
    )

    cluster_summary = summarize_clusters(df_clu)

    tag = f"min{args.min_tissues}.{args.mode}"
    if args.mode == "top":
        tag += f".top{args.top_n}"
    if args.mode == "balanced":
        tag += f".perRegion{args.n_per_region}"
    tag += f".k{args.n_clusters}.{args.cluster_method}"

    plot_horizontal_clustered(
        df=df_clu,
        Z=Z,
        out_base=f"{args.prefix}.proportion_hclust.{tag}",
        title=f"RBP transcript-region composition across tissues (hierarchical clustering, k={args.n_clusters})"
    )

    df_clu.to_csv(f"{args.prefix}.proportion_hclust.{tag}.csv", index=False)
    cluster_summary.to_csv(f"{args.prefix}.cluster_summary.{tag}.csv", index=False)

    print("[OK] Saved:")
    print(f"  {args.prefix}.proportion_hclust.{tag}.png/.pdf")
    print(f"  {args.prefix}.proportion_hclust.{tag}.csv")
    print(f"  {args.prefix}.cluster_summary.{tag}.csv")


if __name__ == "__main__":
    main()
