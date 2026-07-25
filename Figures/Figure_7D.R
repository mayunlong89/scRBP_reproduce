#!/usr/bin/env Rscript
## Fig. 7d — SCZ rare-variant regulons converge predominantly on excitatory neurons.
## Bipartite arc diagram.
##   left  nodes = regulons (RBP + region, "RBP_region")   [11 regulons]
##   right nodes = cell types                              [5 cell types]
##   each arc    = one significant regulon x cell-type association [16 links]
##       arc COLOUR   = transcript region  (node dot coloured the same)
##       arc WIDTH    = TRS (signal strength)
##   right-node SIZE  = number of regulons converging on that cell type
## Only triple-significant associations are shown (FDR_RAS/RGS/TRS all < 0.1).
## RBFOX1_CDS is flagged red -> Fig. 7f,g.
## Reads panelD_SCZ_convergence_DATA.csv (same folder). Run:
##   Rscript panelD_SCZ_convergence_arc.R
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr)
})

## ------------------------- CONFIG (edit freely) -------------------------
here <- tryCatch(dirname(normalizePath(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) getwd())
if (length(here) == 0 || is.na(here)) here <- getwd()
data_csv <- file.path(here, "panelD_SCZ_convergence_DATA.csv")
out_png  <- file.path(here, "figures", "panelD_SCZ_convergence_arc_R.png")
out_pdf  <- file.path(here, "figures", "panelD_SCZ_convergence_arc_R.pdf")

CT_SHORT <- c("Newborn Neuron"="NbN","IPC"="IPC","Div"="Div",
              "Excitatory Neuron"="ExN","Inhibitory Neuron"="InN")
ct_order <- c("NbN","IPC","Div","ExN","InN")             # progenitor -> neuron
REGION_COLS <- c("5'UTR"="#E69F00","3'UTR"="#0072B2","CDS"="#009E73")
highlight_regulon <- "RBFOX1_CDS"; hi_col <- "#D7191C"
curvature <- -0.22
lw_range  <- c(0.35, 1.7)
## ------------------------------------------------------------------------

g <- read_csv(data_csv, show_col_types = FALSE)
g$ct  <- factor(CT_SHORT[g$cell_type], levels = ct_order)
g$reg <- g$regulon_lab

## right column: cell-type nodes, size = # regulons converging
ct_nodes <- g %>% group_by(ct) %>%
  summarise(n_reg = n_distinct(reg), .groups = "drop") %>%
  arrange(match(ct, ct_order))
ct_nodes$y <- seq(1, 0, length.out = nrow(ct_nodes))
ct_y <- setNames(ct_nodes$y, as.character(ct_nodes$ct))

## left column: regulon nodes, ordered by barycenter of targets (few crossings)
reg_bary <- g %>% mutate(ty = ct_y[as.character(ct)]) %>%
  group_by(reg) %>%
  summarise(bary = mean(ty), mz = max(z_TRS), region = region[1], .groups = "drop") %>%
  arrange(desc(bary), desc(mz))
reg_bary$y <- seq(1, 0, length.out = nrow(reg_bary))
reg_y <- setNames(reg_bary$y, reg_bary$reg)

xL <- 0; xR <- 1
links <- g %>% mutate(x = xL, xend = xR, y = reg_y[reg], yend = ct_y[as.character(ct)])

reg_df <- reg_bary; reg_df$x <- xL
reg_df$col  <- ifelse(reg_df$reg == highlight_regulon, hi_col, "#333333")
reg_df$face <- ifelse(reg_df$reg == highlight_regulon, "bold", "plain")
ct_df <- ct_nodes; ct_df$x <- xR

p <- ggplot() +
  geom_curve(data = links,
             aes(x = x, y = y, xend = xend, yend = yend, color = region, linewidth = z_TRS),
             curvature = curvature, alpha = 0.75, lineend = "round") +
  geom_point(data = reg_df, aes(x, y, fill = region), shape = 21,
             size = 2.8, stroke = 0.4, color = "white") +
  geom_text(data = reg_df, aes(x, y, label = reg), hjust = 1, nudge_x = -0.035,
            size = 3.0, color = reg_df$col, fontface = reg_df$face) +
  geom_point(data = ct_df, aes(x, y, size = n_reg), color = "#444444") +
  geom_text(data = ct_df, aes(x, y, label = sprintf("%s (%d)", ct, n_reg)),
            hjust = 0, nudge_x = 0.06, size = 3.4, fontface = "bold") +
  scale_color_manual(values = REGION_COLS, name = "Region",
                     breaks = c("5'UTR","3'UTR","CDS")) +
  scale_fill_manual(values = REGION_COLS, guide = "none") +
  scale_linewidth(range = lw_range, name = expression(TRS~z), breaks = c(3, 4)) +
  scale_size_area(max_size = 11, guide = "none") +
  scale_x_continuous(limits = c(-0.52, 1.42)) +
  scale_y_continuous(limits = c(-0.08, 1.08)) +
  labs(title = "SCZ rare-variant regulons converge predominantly on ExN",
       subtitle = "16 significant regulon -> cell-type links | 11 regulons | 5 cell types") +
  guides(color = guide_legend(order = 1, override.aes = list(linewidth = 2)),
         linewidth = guide_legend(order = 2)) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 9, color = "#666666", hjust = 0.5),
        legend.position = "bottom", legend.box = "horizontal",
        legend.key.size = unit(0.8, "lines"),
        plot.margin = margin(6, 6, 6, 6))

ggsave(out_png, p, width = 6.0, height = 5.0, dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 6.0, height = 5.0, bg = "white")
cat("saved:", out_png, "\n"); cat("saved:", out_pdf, "\n")
