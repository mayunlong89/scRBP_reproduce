#!/usr/bin/env Rscript
## Panel b — SCZ-relevant regulons based on rare variants (dot plot).
## Rows = regulon (RBP + transcript region, "RBP_region"); columns = cell type.
## One dot per significant regulon x cell-type association (triple-significant:
## FDR_RAS < 0.1, FDR_RGS < 0.1 and FDR_TRS < 0.1).
##   dot color = TRS   (trait-regulon score, viridis)
##   dot size  = -log10(min FDR across RAS/RGS/TRS)
## A region color strip is drawn to the left of each regulon label.
## 16 significant regulon x cell-type pairs, 11 regulons, 9 RBPs, 5 cell types.
## Reads panelB_SCZ_rare_dotplot_DATA.csv (same folder). Run:
##   Rscript panelB_SCZ_rare_dotplot.R
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(viridis)
})

## ------------------------- CONFIG (edit freely) -------------------------
here <- tryCatch(dirname(normalizePath(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) getwd())
if (length(here) == 0 || is.na(here)) here <- getwd()
data_csv <- file.path(here, "panelB_SCZ_rare_dotplot_DATA.csv")
out_png  <- file.path(here, "figures", "panelB_SCZ_rare_dotplot_R.png")
out_pdf  <- file.path(here, "figures", "panelB_SCZ_rare_dotplot_R.pdf")

z_min <- 3; z_max <- 7                                   # TRS color range
size_range <- c(2.0, 6.5)                                # bubble radius (-log10 FDR)
CT_SHORT <- c("Newborn Neuron" = "NbN", "IPC" = "IPC", "Div" = "Div",
              "Excitatory Neuron" = "ExN", "Inhibitory Neuron" = "InN")
ct_order <- c("NbN", "IPC", "Div", "ExN", "InN")         # progenitor -> neuron
REGION_COLS <- c("5'UTR" = "#F4A6A0", "CDS" = "#1B9E77",
                 "Intron" = "#8DD3A8", "3'UTR" = "#8C7CB8")
## ------------------------------------------------------------------------

g <- read_csv(data_csv, show_col_types = FALSE)
g$ct <- factor(CT_SHORT[g$cell_type], levels = ct_order)

## regulon row order: by breadth (#cell types) then max TRS; broad/strong at top
lev <- g %>% group_by(regulon_lab) %>%
  summarise(nb = n_distinct(cell_type), mz = max(z_TRS),
            region = region[1], .groups = "drop") %>%
  arrange(nb, mz)
g$regulon_lab <- factor(g$regulon_lab, levels = lev$regulon_lab)
g$z_clip <- pmin(g$z_TRS, z_max)

## region color strip: one tile per regulon, placed just left of the plotting area
strip <- lev %>% mutate(regulon_lab = factor(regulon_lab, levels = lev$regulon_lab))
x_strip <- 0.35                                          # left of first cell-type column (x=1)

p <- ggplot(g, aes(x = ct, y = regulon_lab)) +
  ## region color strip + points
  geom_tile(data = strip, aes(x = x_strip, y = regulon_lab, fill = region),
            width = 0.28, height = 0.72, inherit.aes = FALSE) +
  geom_point(aes(size = neglog10_minq, color = z_clip), stroke = 0.3) +
  scale_color_viridis_c(name = "TRS", limits = c(z_min, z_max),
                        oob = scales::squish) +
  scale_size_continuous(name = expression(-log[10]*"(FDR)"), range = size_range) +
  scale_fill_manual(values = REGION_COLS, name = "Regions",
                    breaks = c("5'UTR", "CDS", "Intron", "3'UTR")) +
  scale_x_discrete(limits = ct_order, expand = expansion(add = c(0.9, 0.6))) +
  labs(title = "SCZ-relevant regulons based on rare variants",
       x = NULL, y = NULL) +
  guides(fill = guide_legend(override.aes = list(size = 0), order = 3)) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 10, color = "#222222"),
    axis.text.y = element_text(size = 9,  color = "#222222"),
    panel.grid.major = element_line(color = "#F0F0F0", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    legend.key.height = unit(0.7, "lines")
  )

ggsave(out_png, p, width = 6.0, height = 4.6, dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 6.0, height = 4.6, bg = "white")
cat("saved:", out_png, "\n"); cat("saved:", out_pdf, "\n")
