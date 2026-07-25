#!/usr/bin/env Rscript
## Fig. 7e — Shared ASD-SCZ rare-variant regulons are concentrated in ExN.
## Dumbbell / paired-lollipop plot of the regulons that reach rare-variant
## significance in BOTH ASD and SCZ in at least one shared cell type.
##   one row per shared regulon (RBP + region)
##   two dots  = rare-variant TRS in ASD (green) and SCZ (purple), evaluated at
##               the SAME shared cell type so the pair is directly comparable
##   bar       = the ASD-SCZ difference
##   right column lists all cell types significant in both disorders
## 9 shared regulons from 7 RBPs. RBFOX1_CDS is flagged red -> Fig. 7f,g.
## Reads panelE_ASD_SCZ_shared_DATA.csv (same folder). Run:
##   Rscript panelE_ASD_SCZ_shared_dumbbell.R
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
})

## ------------------------- CONFIG (edit freely) -------------------------
here <- tryCatch(dirname(normalizePath(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) getwd())
if (length(here) == 0 || is.na(here)) here <- getwd()
data_csv <- file.path(here, "panelE_ASD_SCZ_shared_DATA.csv")
out_png  <- file.path(here, "figures", "panelE_ASD_SCZ_shared_R.png")
out_pdf  <- file.path(here, "figures", "panelE_ASD_SCZ_shared_R.pdf")

COL_ASD <- "#1B9E77"; COL_SCZ <- "#7570B3"
highlight_regulon <- "RBFOX1_CDS"; hi_col <- "#D7191C"
## ------------------------------------------------------------------------

d <- read_csv(data_csv, show_col_types = FALSE)
d <- d %>% mutate(mid = (ASD_z + SCZ_z) / 2) %>% arrange(mid)
d$regulon_lab <- factor(d$regulon_lab, levels = d$regulon_lab)   # weakest at bottom
d$ct_txt <- d$shared_ct

x_lo <- 2.5; x_hi <- 4.6
x_ct <- x_hi + 0.06

long <- d %>% select(regulon_lab, ASD_z, SCZ_z) %>%
  pivot_longer(c(ASD_z, SCZ_z), names_to = "disease", values_to = "z") %>%
  mutate(disease = sub("_z", "", disease))

ylab_cols <- ifelse(levels(d$regulon_lab) == highlight_regulon, hi_col, "#222222")
ylab_face <- ifelse(levels(d$regulon_lab) == highlight_regulon, "bold", "plain")

p <- ggplot() +
  geom_segment(data = d, aes(x = ASD_z, xend = SCZ_z, y = regulon_lab, yend = regulon_lab),
               color = "#BBBBBB", linewidth = 1.1) +
  geom_point(data = long, aes(x = z, y = regulon_lab, color = disease), size = 3.4) +
  geom_vline(xintercept = x_hi + 0.02, color = "#DDDDDD", linewidth = 0.3) +
  geom_text(data = d, aes(x = x_ct, y = regulon_lab, label = ct_txt),
            hjust = 0, size = 3.0, color = "#444444") +
  annotate("text", x = x_ct, y = nlevels(d$regulon_lab) + 0.7,
           label = "cell type", hjust = 0, size = 3.0, fontface = "italic", color = "#888888") +
  scale_color_manual(values = c(ASD = COL_ASD, SCZ = COL_SCZ), name = NULL) +
  scale_x_continuous(limits = c(x_lo, x_hi + 0.75), breaks = seq(2.5, 4.5, 0.5)) +
  coord_cartesian(clip = "off") +
  labs(title = "Shared ASD-SCZ rare-variant regulons are concentrated in ExN",
       subtitle = "9 regulons significant in both, shown at their shared cell type (7 RBPs)",
       x = "Rare-variant TRS", y = "shared regulon") +
  theme_classic(base_size = 11) +
  theme(axis.line = element_line(color = "#333333", linewidth = 0.4),
        axis.ticks = element_line(color = "#333333", linewidth = 0.4),
        panel.grid.major.x = element_line(color = "#F0F0F0", linewidth = 0.3),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 9, color = ylab_cols, face = ylab_face),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "#666666"),
        legend.position = "top", legend.justification = "right",
        legend.background = element_blank(), legend.text = element_text(size = 10),
        legend.margin = margin(0, 0, -4, 0),
        plot.margin = margin(6, 78, 6, 6))

ggsave(out_png, p, width = 6.6, height = 4.3, dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 6.6, height = 4.3, bg = "white")
cat("saved:", out_png, "\n"); cat("saved:", out_pdf, "\n")
