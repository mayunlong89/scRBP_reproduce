#!/usr/bin/env Rscript
## Fig. 7a — ASD rare-variant risk converges on neuronal and progenitor regulons.
## Dot plot of every cell type-resolved RBP regulon significantly associated with
## ASD rare-variant burden. Rows = regulon (RBP + transcript region, "RBP_region"),
## columns = cell type. A dot is drawn only where the association is triple-
## significant (FDR_RAS < 0.1, FDR_RGS < 0.1 and FDR_TRS < 0.1).
##   dot color = TRS  (trait-regulon score, viridis)
##   dot size  = -log10(min FDR across RAS/RGS/TRS)
## A region color strip sits to the left of every regulon label.
## 164 significant regulon x cell-type associations, 71 regulons, 54 RBPs, 14 cell types.
## Reads panelA_ASD_regulon_dotplot_DATA.csv (same folder). Run:
##   Rscript panelA_ASD_regulon_dotplot.R
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(viridis)
})

## ------------------------- CONFIG (edit freely) -------------------------
here <- tryCatch(dirname(normalizePath(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) getwd())
if (length(here) == 0 || is.na(here)) here <- getwd()
data_csv <- file.path(here, "panelA_ASD_regulon_dotplot_DATA.csv")
out_png  <- file.path(here, "figures", "panelA_ASD_regulon_dotplot_R.png")
out_pdf  <- file.path(here, "figures", "panelA_ASD_regulon_dotplot_R.pdf")

z_min <- 3; z_max <- 7                                   # TRS color range
size_range <- c(1.2, 4.5)                                # bubble radius (-log10 FDR)
## cell-type column order (progenitor/dividing -> neuron -> glia -> vascular)
CT_SHORT <- c("Newborn Neuron"="NbN","IPC"="IPC","Excitatory Neuron"="ExN",
              "RG.div"="RG.div","Div"="Div","Endothelial"="Endo","RG"="RG",
              "Astrocyte"="AS","Microglia"="MG","Pericyte"="Peri",
              "Inhibitory Neuron"="InN","OPC.div"="OPC.div","OPC"="OPC",
              "Oligodendrocyte"="ODC")
ct_order <- c("NbN","IPC","ExN","RG.div","Div","Endo","RG","AS","MG","Peri",
              "InN","OPC.div","OPC","ODC")
REGION_COLS <- c("5'UTR"="#F4A6A0","CDS"="#1B9E77","Intron"="#8DD3A8","3'UTR"="#8C7CB8")
## ------------------------------------------------------------------------

pretty_region <- function(r) {
  r <- sub("^3UTR$","3'UTR",r); r <- sub("^5UTR$","5'UTR",r); sub("^Introns$","Intron",r)
}

g <- read_csv(data_csv, show_col_types = FALSE)
g$ct  <- factor(CT_SHORT[g$cell_type], levels = ct_order)
g$reg <- sprintf("%s_%s", g$RBP, pretty_region(g$region))
g$region_p <- pretty_region(g$region)
g$z_clip <- pmin(g$z_TRS, z_max)
ct_rank <- setNames(seq_along(ct_order), ct_order)

## regulon row order: by barycenter of the cell types it hits (few crossings),
## so neuronal/progenitor-heavy regulons cluster near the top columns.
ord <- g %>% mutate(cr = ct_rank[as.character(ct)]) %>%
  group_by(reg) %>%
  summarise(bary = mean(cr), mz = max(z_TRS), region_p = region_p[1], .groups = "drop") %>%
  arrange(desc(bary), desc(mz))
g$reg <- factor(g$reg, levels = ord$reg)               # top row = first level -> drawn at top

## region color strip: one tile per regulon, just left of the first column
strip <- ord %>% mutate(reg = factor(reg, levels = ord$reg))
x_strip <- 0.30

p <- ggplot(g, aes(x = ct, y = reg)) +
  geom_tile(data = strip, aes(x = x_strip, y = reg, fill = region_p),
            width = 0.30, height = 0.75, inherit.aes = FALSE) +
  geom_point(aes(size = neglog10_minq, color = z_clip), stroke = 0.2) +
  scale_color_viridis_c(name = "TRS", limits = c(z_min, z_max), oob = scales::squish) +
  scale_size_continuous(name = expression(-log[10]*"(FDR)"), range = size_range) +
  scale_fill_manual(values = REGION_COLS, name = "Regions",
                    breaks = c("5'UTR","CDS","Intron","3'UTR")) +
  scale_x_discrete(limits = ct_order, expand = expansion(add = c(0.9, 0.6))) +
  labs(title = "ASD rare-variant risk converges on neuronal and progenitor regulons",
       x = NULL, y = NULL) +
  guides(fill = guide_legend(override.aes = list(size = 0), order = 3)) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5, color = "#222222"),
    axis.text.y = element_text(size = 5.6, color = "#222222"),
    panel.grid.major = element_line(color = "#F2F2F2", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 9.5),
    legend.key.height = unit(0.65, "lines")
  )

ggsave(out_png, p, width = 6.2, height = 10.5, dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 6.2, height = 10.5, bg = "white")
cat("saved:", out_png, "\n"); cat("saved:", out_pdf, "\n")
