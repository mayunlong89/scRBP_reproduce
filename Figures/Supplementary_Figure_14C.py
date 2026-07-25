#!/usr/bin/env python3
# =============================================================================
# Supplementary Fig. 14c: RBP-regulon targets respond preferentially to their
# cognate RBP.
#
# Specificity = pooled enrichment of a regulon's targets among the responsive
# genes of its OWN (cognate) RBP, divided by the enrichment against non-cognate
# RBPs. Because the same target set is scored against matched vs mismatched
# perturbations, expression / generic-responsiveness biases cancel.
#
#   - one forest panel, hematopoietic lineages on the y-axis (all in black)
#   - four response-gene definitions per lineage, encoded by a Reds ramp
#     (top-100 darkest -> DE lightest)
#   - points = pooled specificity ratio ; bars = 95% bootstrap CI (over regulons)
#   - text   = BH-FDR from a 5,000x RBP-identity permutation test ("ns" if > 0.05)
#
# Input : comprehensive/comprehensive_3background_wDE.tsv  (background == "specificity")
# Output: figures/panel_c_specificity_standalone.{png,pdf}
# Run   : srun --partition=defq --mem=4G --time=00:15:00 python3 48_panel_c_specificity.py
# =============================================================================
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- paths -----------------------------------------------------------------
BASE = "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/00_scRBP_perturb-seq_results/results_regulon_validation_K562_20260709"
TSV  = f"{BASE}/comprehensive/comprehensive_3background_wDE.tsv"
OUT  = f"{BASE}/comprehensive/figures/panel_c_specificity_standalone"

# ---- figure text -----------------------------------------------------------
TITLE  = "Validation from K562 perturb-seq data (matched vs. non-cognate RBPs)"
XLABEL = "Specificity\n(matched / non-cognate RBP)"

# ---- palette ---------------------------------------------------------------
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"
Reds = plt.cm.Reds
# (reference, y-offset, marker size, line width, colour)
CUT = [
    ("100", 0.27, 6.5, 2.6, Reds(0.90)),
    ("200", 0.09, 5.5, 2.2, Reds(0.68)),
    ("500", -0.09, 4.6, 1.9, Reds(0.48)),
    ("DE",  -0.27, 3.9, 1.6, Reds(0.30)),
]
CUTLAB = {"100": "top-100", "200": "top-200", "500": "top-500", "DE": "DE (FDR<0.05)"}

ORDER = ["All hematopoietic", "Erythroid_Mega", "Myeloid", "HSPC_Stromal", "B", "T_NK"]
PRETTY = {"All hematopoietic": "All hematopoietic", "Erythroid_Mega": "Erythroid / Mega",
          "Myeloid": "Myeloid", "HSPC_Stromal": "HSPC / Stromal", "B": "B", "T_NK": "T / NK"}

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": "#4d4d4d", "axes.linewidth": .8, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "font.size": 9.5,
    "axes.spines.top": False, "axes.spines.right": False, "grid.color": GRID,
    "grid.linewidth": .5, "pdf.fonttype": 42, "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
})

# ---- load specificity rows -------------------------------------------------
comp = {(r["reference"], r["group"]): r
        for r in csv.DictReader(open(TSV), delimiter="\t")
        if r["background"] == "specificity"}

ypos = {g: i for i, g in enumerate(ORDER[::-1])}

# ---- draw ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.2, 4.4))
xlo = min(float(comp[(ref, g)]["ci_lo"]) for ref, *_ in CUT for g in ORDER)
xhi = max(float(comp[(ref, g)]["ci_hi"]) for ref, *_ in CUT for g in ORDER)

for g in ORDER:
    y = ypos[g]
    for ref, dy, ms, lw, col in CUT:
        r = comp[(ref, g)]
        v, lo, hi, fdr = float(r["value"]), float(r["ci_lo"]), float(r["ci_hi"]), float(r["FDR"])
        ax.plot([lo, hi], [y + dy, y + dy], color=col, lw=lw, solid_capstyle="round", zorder=3)
        ax.plot(v, y + dy, "o", ms=ms, color=col, mec="white", mew=.7, zorder=4)
        # FDR text label (value if significant, otherwise "ns")
        label = f"FDR = {fdr:.2g}" if fdr < 0.05 else "ns"
        ax.text(hi + (xhi - xlo) * 0.015, y + dy, label, va="center", ha="left",
                fontsize=6.0, color=INK2)

ax.axvline(1, color=MUTED, ls="--", lw=1, zorder=1)                      # ratio = 1 (no specificity)
ax.axhspan(ypos["All hematopoietic"] - .45, ypos["All hematopoietic"] + .45,
           color="#00000008", zorder=0)
ax.set_xlim(min(0.95, xlo - (xhi - xlo) * 0.05), xhi + (xhi - xlo) * 0.30)  # extra room for FDR text
ax.set_ylim(-.6, len(ORDER) - .4)
ax.set_yticks([ypos[g] for g in ORDER])
ax.set_yticklabels([PRETTY[g] for g in ORDER])
for t, g in zip(ax.get_yticklabels(), ORDER):
    t.set_color(INK)
    t.set_weight("bold" if g == "All hematopoietic" else "normal")
ax.set_xlabel(XLABEL)
ax.grid(axis="x", zorder=0)
ax.tick_params(axis="y", length=0)
ax.set_title(TITLE, weight="bold", fontsize=10, loc="left", pad=8)

for ref, _, ms, _, col in CUT:
    ax.plot([], [], "o", color=col, ms=ms, label=CUTLAB[ref])
ax.legend(fontsize=7.4, frameon=False, loc="lower right",
          title="Reference set", title_fontsize=7.6)

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=350, bbox_inches="tight")
print(f"wrote {OUT}.png / .pdf")
