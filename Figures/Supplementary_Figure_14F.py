#!/usr/bin/env python3
# =============================================================================
# Supplementary Fig. 14f: canonical RBPs recovered by their own K562 perturbation.
#
# Cognate self-rank for the 15 canonical RBPs with the highest values. For each
# RBP, targets were aggregated across its cell-type- and transcript-region-
# specific regulons and scored (top-1000 cutoff AUROC) against all 496 RBP
# perturbations; self-rank is the percentile of the cognate perturbation
# (0.5 = chance). Lollipops are coloured by RNA-processing function; disease
# drivers are shown in italic.
#
# Input : sectionA/perRBP_selfrank_top1000.tsv   (columns: RBP, self_rank)
# Output: sectionA/figures/panel_f_canonical_rbps_standalone.{png,pdf}
# Run   : srun --partition=defq --mem=4G --time=00:15:00 python3 51_panel_f_canonical_rbps.py
# =============================================================================
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ---- paths -----------------------------------------------------------------
BASE = "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/00_scRBP_perturb-seq_results/results_regulon_validation_K562_20260709"
TSV  = f"{BASE}/sectionA/perRBP_selfrank_top1000.tsv"
OUT  = f"{BASE}/sectionA/figures/panel_f_canonical_rbps_standalone"

TITLE = "Top canonical RBPs validated by their own K562 perturbation"
XLABEL = "cognate self-rank  (per RBP, aggregated targets)"

# ---- curated 15 canonical RBPs (top self-rank among canonical factors) ------
# Displayed top -> bottom.
ORDER = ["PABPC1", "FBL", "DDX3X", "SFPQ", "DDX6", "SF3B1", "LSM6", "TRA2B",
         "RBM17", "PPIG", "SRSF1", "TARDBP", "DDX5", "NUDT21", "TAF15"]
# RNA-processing function -> colour
SPL, STB, HEL = "#2a78d6", "#1baf7a", "#eb6834"
CLASS = {
    "SF3B1": SPL, "SFPQ": SPL, "TRA2B": SPL, "SRSF1": SPL, "RBM17": SPL,
    "PPIG": SPL, "LSM6": SPL, "TAF15": SPL,
    "PABPC1": STB, "DDX6": STB, "TARDBP": STB, "NUDT21": STB,
    "DDX3X": HEL, "DDX5": HEL, "FBL": HEL,
}
DISEASE = {"SF3B1", "TARDBP", "TAF15", "DDX3X"}   # recurrent disease drivers (italic)

# ---- palette ---------------------------------------------------------------
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": "#4d4d4d", "axes.linewidth": .8, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False, "grid.color": GRID,
    "grid.linewidth": .5, "pdf.fonttype": 42, "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
})

# ---- load self-ranks -------------------------------------------------------
sr = {r["RBP"]: float(r["self_rank"]) for r in csv.DictReader(open(TSV), delimiter="\t")}

# ---- draw ------------------------------------------------------------------
lst = ORDER[::-1]                       # bottom -> top for plotting
yv = np.arange(len(lst))
fig, ax = plt.subplots(figsize=(7.0, 6.0))
for i, g in enumerate(lst):
    v = sr[g]; col = CLASS[g]
    ax.plot([0.5, v], [i, i], color=col, lw=2.4, alpha=.55, zorder=2)   # stem from chance
    ax.plot(v, i, "o", ms=9, color=col, mec="white", mew=.9, zorder=3)  # marker

ax.axvline(0.5, color=MUTED, ls="--", lw=1, zorder=1)                   # chance = 0.5
ax.set_yticks(yv); ax.set_yticklabels(lst, fontsize=10)
for t, g in zip(ax.get_yticklabels(), lst):
    t.set_color(CLASS[g]); t.set_weight("bold")
    if g in DISEASE:
        t.set_style("italic")
ax.set_xlim(0.48, 1.02); ax.set_ylim(-.6, len(lst) - .4)
ax.set_xlabel(XLABEL)
ax.set_title(TITLE, weight="bold", loc="left")

legend = [Patch(color=SPL, label="splicing / spliceosome"),
          Patch(color=STB, label="mRNA stability / decay"),
          Patch(color=HEL, label="helicase / rRNA")]
ax.legend(handles=legend, fontsize=8, frameon=False, loc="lower right",
          title="Function", title_fontsize=8.2)

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=350, bbox_inches="tight")
print(f"wrote {OUT}.png / .pdf")
