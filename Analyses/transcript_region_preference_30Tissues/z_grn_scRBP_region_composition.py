#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import linkage, leaves_list


REGIONS = ["3UTR", "5UTR", "CDS", "Introns"]


def infer_tissue_name(tissue_dir: str) -> str:
    """z_GRNBoost2_Liver_30times -> Liver"""
    base = os.path.basename(tissue_dir.rstrip("/"))
    m = re.match(r"z_GRNBoost2_(.+)_30times$", base)
    return m.group(1) if m else base


def parse_merged_gmt_presence(gmt_path: str) -> pd.DataFrame:
    """
    Parse scRBP_allRegions_min10genes.symbol_default.gmt
    Each line: setName \t desc \t gene1...
    setName endswith _3UTR/_5UTR/_CDS/_Introns

    Return a table with columns: rbp, region, present(0/1)
    Presence is 1 if region line exists for that rbp (regardless of #genes).
    """
    seen = set()
    rows = []

    with open(gmt_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            set_name = parts[0]

            region = None
            for r in REGIONS:
                if set_name.endswith("_" + r):
                    region = r
                    break
            if region is None:
                continue

            # rbp name: drop _REGION suffix then drop optional _regulon
            rbp = set_name[: -(len(region) + 1)]
            rbp = re.sub(r"_regulon$", "", rbp)

            key = (rbp, region)
            if key in seen:
                continue
            seen.add(key)
            rows.append({"rbp": rbp, "region": region, "present": 1})

    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=["rbp"] + REGIONS)

    pres = (
        df.pivot_table(index="rbp", columns="region", values="present", aggfunc="max", fill_value=0)
        .reindex(columns=REGIONS, fill_value=0)
        .reset_index()
    )
    return pres


def cluster_order(mat: np.ndarray, method: str = "average", metric: str = "euclidean"):
    """Return row order indices for clustering."""
    if mat.shape[0] <= 2:
        return np.arange(mat.shape[0])
    Z = linkage(mat, method=method, metric=metric)
    return leaves_list(Z)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base_dir", required=True, help="Parent folder containing z_GRNBoost2_*_30times")
    ap.add_argument("--pattern", default="z_GRNBoost2_*_30times", help="Glob under base_dir")
    ap.add_argument(
        "--gmt_name",
        default="scRBP_allRegions_min10genes.symbol_default.gmt",
        help="Merged GMT name inside each tissue folder",
    )
    ap.add_argument("--top_n", type=int, default=50, help="Top N RBPs to show in stacked barplot")
    ap.add_argument("--out_prefix", default="scRBP_regionComposition", help="Prefix for outputs")
    args = ap.parse_args()

    tissue_dirs = sorted(glob.glob(os.path.join(args.base_dir, args.pattern)))
    if not tissue_dirs:
        raise SystemExit(f"[ERROR] No tissue dirs found: {args.base_dir}/{args.pattern}")

    # --------- collect per tissue presence ---------
    all_tissues = []
    all_rbps = set()
    tissue_presence = {}

    for tdir in tissue_dirs:
        tissue = infer_tissue_name(tdir)
        gmt_path = os.path.join(tdir, args.gmt_name)
        if not os.path.isfile(gmt_path):
            print(f"[WARN] missing GMT, skip: {gmt_path}")
            continue

        pres = parse_merged_gmt_presence(gmt_path)  # columns: rbp + regions (0/1)
        pres["tissue"] = tissue
        all_tissues.append(tissue)
        all_rbps.update(pres["rbp"].tolist())
        tissue_presence[tissue] = pres

    all_tissues = sorted(list(set(all_tissues)))
    all_rbps = sorted(list(all_rbps))

    if not all_tissues or not all_rbps:
        raise SystemExit("[ERROR] Parsed nothing. Check base_dir/gmt_name.")

    # --------- build tissue × rbp full matrix (0/1 presence per region) ---------
    # Long -> wide full
    full_rows = []
    for tissue in all_tissues:
        pres = tissue_presence[tissue].set_index("rbp")
        pres = pres.reindex(all_rbps, fill_value=0)
        pres = pres[REGIONS].astype(int)
        pres["tissue"] = tissue
        pres["rbp"] = pres.index
        full_rows.append(pres.reset_index(drop=True))

    tissue_rbp = pd.concat(full_rows, ignore_index=True)  # tissue, rbp, region cols

    # --------- across tissues: counts (0-30) ---------
    rbp_counts = tissue_rbp.groupby("rbp")[REGIONS].sum().reset_index()
    rbp_counts.rename(columns={r: f"nTissues_{r}" for r in REGIONS}, inplace=True)

    rbp_counts["nTissues_anyRegion"] = tissue_rbp.groupby("rbp").apply(
        lambda df: int((df[REGIONS].sum(axis=1) > 0).sum())
    ).values

    nT = len(all_tissues)
    rbp_counts["n_tissues_total"] = nT

    # --------- composition proportions: sum to 1 (for RBPs that appear at least once) ---------
    count_cols = [f"nTissues_{r}" for r in REGIONS]
    rbp_counts["total_region_occurrences"] = rbp_counts[count_cols].sum(axis=1)

    for r in REGIONS:
        rbp_counts[f"prop_{r}"] = rbp_counts[f"nTissues_{r}"] / rbp_counts["total_region_occurrences"].replace({0: np.nan})

    # filter RBPs that never appear in any region across all tissues (rare, but safe)
    rbp_counts_use = rbp_counts[rbp_counts["total_region_occurrences"] > 0].copy()

    # --------- save tables ---------
    out_counts = f"{args.out_prefix}.rbp_counts_and_composition.csv"
    out_tissue_rbp = f"{args.out_prefix}.tissue_by_rbp_presence.csv"
    rbp_counts_use.to_csv(out_counts, index=False)
    tissue_rbp.to_csv(out_tissue_rbp, index=False)
    print(f"[OK] wrote: {out_counts}")
    print(f"[OK] wrote: {out_tissue_rbp}")

    # ============================
    # Plot 1: Top-N stacked bar (composition)
    # ============================
    top = rbp_counts_use.sort_values("total_region_occurrences", ascending=False).head(args.top_n).copy()
    top = top.set_index("rbp")[[f"prop_{r}" for r in REGIONS]]

    ax = top.plot(kind="bar", stacked=True, figsize=(14, 5))
    ax.set_ylabel("Region composition (sums to 1)")
    ax.set_xlabel(f"RBP (top {args.top_n} by total region-occurrences across tissues)")
    ax.set_title("RBP binding-region composition across tissues (composition)")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}.plot1_top{args.top_n}_stacked_composition.png", dpi=200)
    plt.close()
    print(f"[OK] saved: {args.out_prefix}.plot1_top{args.top_n}_stacked_composition.png")

    # ============================
    # Plot 2: Cluster heatmap (composition)
    # ============================
    # Use more RBPs for clustering, e.g. those with at least K total occurrences
    # Here: keep RBPs that appear in >= 5 region-occurrences across tissues (adjustable if you want)
    K = 5
    hm = rbp_counts_use[rbp_counts_use["total_region_occurrences"] >= K].copy()
    hm = hm.set_index("rbp")[[f"prop_{r}" for r in REGIONS]].fillna(0.0)

    mat = hm.values
    order = cluster_order(mat, method="average", metric="euclidean")
    hm_ord = hm.iloc[order]

    fig, ax = plt.subplots(figsize=(6, max(6, 0.16 * hm_ord.shape[0])))
    im = ax.imshow(hm_ord.values, aspect="auto", interpolation="nearest")

    ax.set_yticks(np.arange(hm_ord.shape[0]))
    ax.set_yticklabels(hm_ord.index, fontsize=6)
    ax.set_xticks(np.arange(len(REGIONS)))
    ax.set_xticklabels(REGIONS)

    ax.set_title(f"Clustered heatmap of region composition (RBPs with total_occurrences ≥ {K})")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Composition proportion")

    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}.plot2_cluster_heatmap_composition_K{K}.png", dpi=200)
    plt.close()
    print(f"[OK] saved: {args.out_prefix}.plot2_cluster_heatmap_composition_K{K}.png")


if __name__ == "__main__":
    main()
