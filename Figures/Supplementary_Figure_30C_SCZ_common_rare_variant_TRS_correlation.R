#!/usr/bin/env Rscript
## Panel c — SCZ: common vs rare-variant TRS correlation.
## One point per regulon x cell-type pair (id = RBP_region_celltype) tested in SCZ.
##   x = SCZ common-variant TRS,  y = SCZ rare-variant TRS
##   color = significance class:
##       both significant  (red)    common+rare both significant
##       rare-only         (green)  rare-variant model only
##       common-only       (orange) common-variant model only
##       not significant   (grey)
## Rare-variant significance = triple test (FDR_RAS/RGS/TRS all < 0.1);
## common significance = curated common-variant significant list.
## Both-significant and rare-only points are labelled RBP_region (ctAbbr).
## Reads panelC_SCZ_common_vs_rare_DATA.csv (same folder). Run:
##   Rscript panelC_SCZ_common_vs_rare_scatter.R
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(ggrepel)
})

## ------------------------- CONFIG (edit freely) -------------------------
here <- tryCatch(dirname(normalizePath(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) getwd())
if (length(here) == 0 || is.na(here)) here <- getwd()
data_csv <- file.path(here, "panelC_SCZ_common_vs_rare_DATA.csv")
out_png  <- file.path(here, "figures", "panelC_SCZ_common_vs_rare_R.png")
out_pdf  <- file.path(here, "figures", "panelC_SCZ_common_vs_rare_R.pdf")

COLORS <- c(ns = "#CFCFCF", common_only = "#E69F00",
            rare_only = "#1B9E77", both = "#C0392B")
LABELS <- c(ns = "not significant", common_only = "common-only significant",
            rare_only = "rare-only significant", both = "both significant")
SIZES  <- c(ns = 0.7, common_only = 2.0, rare_only = 2.6, both = 2.8)
ALPHAS <- c(ns = 0.35, common_only = 0.85, rare_only = 1.0, both = 1.0)
draw_order <- c("ns", "common_only", "rare_only", "both")
rho_txt <- "rho = 0.24\nP = 2e-108"          # Spearman (ASCII rho for base R fonts)
## ------------------------------------------------------------------------

pretty_region <- function(r) {
  r <- sub("^3UTR$", "3'UTR", r); r <- sub("^5UTR$", "5'UTR", r)
  sub("^Introns$", "Intron", r)
}
CT_ABBR <- c("Excitatory Neuron" = "ExN", "Inhibitory Neuron" = "InN",
             "Newborn Neuron" = "NbN", "IPC" = "IPC", "Div" = "Div")

g <- read_csv(data_csv, show_col_types = FALSE)
g$sig_class <- factor(g$sig_class, levels = draw_order)
g <- g %>% arrange(sig_class)                          # ns first, both drawn on top

## label both-significant + rare-only points
lab <- g %>% filter(sig_class %in% c("both", "rare_only"))
ct_short <- ifelse(lab$cell_type %in% names(CT_ABBR),
                   CT_ABBR[lab$cell_type], lab$cell_type)
lab$txt <- sprintf("%s_%s (%s)", lab$RBP, pretty_region(lab$region), ct_short)

p <- ggplot(g, aes(z_common, z_rare)) +
  geom_hline(yintercept = 0, color = "#DDDDDD", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "#DDDDDD", linewidth = 0.4) +
  geom_point(aes(color = sig_class, size = sig_class, alpha = sig_class),
             stroke = 0.25) +
  geom_point(data = lab, shape = 21, size = 2.9, fill = NA,
             color = "#111111", stroke = 0.4) +
  ggrepel::geom_text_repel(data = lab, aes(label = txt), size = 2.7,
             color = "#111111", min.segment.length = 0,
             segment.size = 0.3, segment.color = "#999999", box.padding = 0.7,
             point.padding = 0.3, force = 12, force_pull = 0.2,
             max.overlaps = Inf, max.iter = 100000, seed = 3) +
  scale_color_manual(values = COLORS, labels = LABELS, name = NULL,
                     breaks = c("both", "rare_only", "common_only", "ns")) +
  scale_size_manual(values = SIZES, guide = "none") +
  scale_alpha_manual(values = ALPHAS, guide = "none") +
  annotate("text", x = min(g$z_common), y = max(g$z_rare),
           label = rho_txt, hjust = 0, vjust = 1, size = 4, color = "#222222") +
  labs(title = "SCZ: common vs rare-variant TRS correlation",
       x = "SCZ TRS (common variant)", y = "SCZ TRS (rare variant)") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#F2F2F2", linewidth = 0.3),
        legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.background = element_blank(), legend.text = element_text(size = 9),
        plot.title = element_text(face = "bold", size = 13))

ggsave(out_png, p, width = 8.4, height = 7.2, dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 8.4, height = 7.2, bg = "white")
cat("saved:", out_png, "\n"); cat("saved:", out_pdf, "\n")
cat("labelled:", nrow(lab), "points\n")
