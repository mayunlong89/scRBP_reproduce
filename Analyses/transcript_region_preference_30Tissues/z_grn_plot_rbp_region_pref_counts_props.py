#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REGIONS = ["3UTR", "5UTR", "CDS", "Introns"]

REGION_COLORS = {
    "CDS": "#2CB34A",
    "3UTR": "#0B6B3A",
    "5UTR": "#A6D854",
    "Introns": "#FFB000",
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
    """
    df_sorted: already sorted & contains dominant_region.
    mode:
      - 'all': return df_sorted (no truncation)
      - 'top': return head(top_n)
      - 'balanced': take n_per_region from each dominant_region (in REGIONS order)
    """
    mode = mode.lower()
    if mode == "all":
        return df_sorted

    if mode == "top":
        return df_sorted.head(top_n).copy()

    if mode == "balanced":
        parts = []
        for r in REGIONS:
            sub = df_sorted[df_sorted["dominant_region"] == r].head(n_per_region).copy()
            parts.append(sub)
        out = pd.concat(parts, ignore_index=True)
        return out

    raise ValueError("mode must be one of: all, top, balanced")


def plot_stacked(df: pd.DataFrame, value_cols: list, title: str, ylabel: str, out_base: str):
    plot_df = df.set_index("rbp")[value_cols].copy()
    plot_df.columns = REGIONS

    fig_w = max(12, 0.22 * plot_df.shape[0])
    fig, ax = plt.subplots(figsize=(fig_w, 5))

    plot_df.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        color=[REGION_COLORS[r] for r in REGIONS],
        width=0.85,
    )

    ax.set_title(title, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xlabel(f"RBP (grouped by dominant region)", fontsize=11)
    plt.xticks(rotation=90, fontsize=8)

    legend_labels = ["3′UTR", "5′UTR", "CDS", "Introns"]
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles, legend_labels, title="Region", frameon=False,
              bbox_to_anchor=(1.01, 1.0), loc="upper left")

    plt.tight_layout()
    _save_fig(fig, out_base)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats_csv", required=True)
    ap.add_argument("--prefix", default="scRBP_RBPRegionPref")
    ap.add_argument("--min_tissues", type=int, default=5)

    # NEW: selection mode
    ap.add_argument("--mode", choices=["all", "top", "balanced"], default="balanced",
                    help="How to select RBPs for plotting")
    ap.add_argument("--top_n", type=int, default=200,
                    help="Used when --mode top (take top_n after sorting)")
    ap.add_argument("--n_per_region", type=int, default=50,
                    help="Used when --mode balanced (take N per dominant region)")

    args = ap.parse_args()

    df = pd.read_csv(args.stats_csv)

    count_cols = [f"nTissues_{r}" for r in REGIONS]
    prop_cols = [f"prop_{r}" for r in REGIONS]
    _require_cols(df, ["rbp", "nTissues_anyRegion"], "stats_csv")
    _require_cols(df, count_cols + prop_cols, "stats_csv")

    # filter
    df = df[df["nTissues_anyRegion"] >= args.min_tissues].copy()
    if df.empty:
        raise SystemExit(f"[ERROR] No RBPs left after filtering nTissues_anyRegion >= {args.min_tissues}")

    df = add_dominant_region(df)

    if "total_region_occurrences" not in df.columns:
        df["total_region_occurrences"] = df[count_cols].sum(axis=1)

    # sort for each plot type
    df_counts_sorted = sort_by_preference(df, metric="count")
    df_props_sorted = sort_by_preference(df, metric="prop")

    # select RBPs
    sel_counts = select_rbplist(df_counts_sorted, args.mode, args.top_n, args.n_per_region)
    sel_props = select_rbplist(df_props_sorted, args.mode, args.top_n, args.n_per_region)

    # labels for filenames
    tag = f"min{args.min_tissues}.{args.mode}"
    if args.mode == "top":
        tag += f".top{args.top_n}"
    if args.mode == "balanced":
        tag += f".perRegion{args.n_per_region}"

    # counts plot
    plot_stacked(
        df=sel_counts,
        value_cols=count_cols,
        title=f"Across 30 tissues: region enrichment counts per RBP (≥{args.min_tissues} tissues; {args.mode})",
        ylabel="# tissues enriched per region (0–30); stacked total (0–120)",
        out_base=f"{args.prefix}.counts_stacked.{tag}",
    )

    # proportion plot
    plot_stacked(
        df=sel_props,
        value_cols=prop_cols,
        title=f"Across 30 tissues: region composition per RBP (≥{args.min_tissues} tissues; {args.mode})",
        ylabel="Composition proportion (3′UTR + 5′UTR + CDS + Introns = 1)",
        out_base=f"{args.prefix}.proportion_stacked.{tag}",
    )

    # save the selected tables
    sel_counts.to_csv(f"{args.prefix}.selected_for_counts.{tag}.csv", index=False)
    sel_props.to_csv(f"{args.prefix}.selected_for_props.{tag}.csv", index=False)

    print("[OK] Saved PNG+PDF:")
    print(f"  {args.prefix}.counts_stacked.{tag}.png/.pdf")
    print(f"  {args.prefix}.proportion_stacked.{tag}.png/.pdf")
    print("[OK] Saved selected tables:")
    print(f"  {args.prefix}.selected_for_counts.{tag}.csv")
    print(f"  {args.prefix}.selected_for_props.{tag}.csv")


if __name__ == "__main__":
    main()
