#!/usr/bin/env python3
# =============================================================================
# Supplementary Fig. 14d: distribution of cognate-RBP self-ranks across regulons.
#
# For each of 1,702 regulons, target genes were scored (top-1000 cutoff AUROC)
# against all 496 RBP perturbations; the self-rank is the percentile of the
# cognate perturbation among all 496 (0.5 = chance). Real regulons are compared
# with a size-matched shuffled-target null.
#
#   - overlaid histograms: shuffled (grey) vs real (blue)
#   - median lines, and a one-sided Wilcoxon signed-rank test of real vs chance
#
# Input : sectionA/sectionA_topN.tsv   (rows with cutoff == "1000")
#         columns used: self_rank, self_rank_shuffled
# Output: sectionA/figures/panel_d_selfrank_hist_standalone.{png,pdf}
# Run   : srun --partition=defq --mem=4G --time=00:15:00 python3 49_panel_d_selfrank_hist.py
# =============================================================================
import csv
import numpy as np
from scipy.stats import wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- paths -----------------------------------------------------------------
BASE = "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/00_scRBP_perturb-seq_results/results_regulon_validation_K562_20260709"
TSV  = f"{BASE}/sectionA/sectionA_topN.tsv"
OUT  = f"{BASE}/sectionA/figures/panel_d_selfrank_hist_standalone"
CUTOFF = "1000"   # which top-N cutoff block to read

XLABEL = "self-rank: cognate perturbation's percentile among all 496 RBPs"

# ---- palette ---------------------------------------------------------------
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"
BLUE, GREY, RED = "#2a78d6", "#bdbbb2", "#e34948"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": "#4d4d4d", "axes.linewidth": .8, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False, "grid.color": GRID,
    "grid.linewidth": .5, "pdf.fonttype": 42, "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
})

# ---- load ------------------------------------------------------------------
rows = [r for r in csv.DictReader(open(TSV), delimiter="\t") if r["cutoff"] == CUTOFF]
sr  = np.array([float(r["self_rank"]) for r in rows])            # real regulons
sr0 = np.array([float(r["self_rank_shuffled"]) for r in rows])   # shuffled-target null
n = len(sr)
p_w = wilcoxon(sr - 0.5, alternative="greater").pvalue          # real vs chance (0.5)

# ---- draw ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.6, 4.0))
bins = np.linspace(0, 1, 31)
ax.hist(sr0, bins=bins, color=GREY, alpha=.80, label=f"shuffled targets (median {np.median(sr0):.2f})", zorder=2)
ax.hist(sr,  bins=bins, color=BLUE, alpha=.72, label=f"real regulons (median {np.median(sr):.2f})", zorder=3)
ax.axvline(0.5, color=MUTED, ls="--", lw=1, zorder=1)            # chance
ax.axvline(np.median(sr), color=BLUE, lw=2.2, zorder=4)          # real median
ax.axvspan(0.9, 1.0, color=RED, alpha=.08, zorder=0)            # top decile

ax.set_xlabel(XLABEL)
ax.set_ylabel("regulons")
ax.legend(fontsize=8, frameon=False, loc="upper left")
ax.text(.98, .97, f"n = {n:,} regulons\nWilcoxon P = {p_w:.0e}",
        transform=ax.transAxes, ha="right", va="top", fontsize=8.5, style="italic")

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=350, bbox_inches="tight")
print(f"wrote {OUT}.png / .pdf   (n={n}, median real={np.median(sr):.3f}, Wilcoxon P={p_w:.1e})")
