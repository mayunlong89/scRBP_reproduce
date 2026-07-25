#!/usr/bin/env Rscript
## Panel a — Proportion of ASD rare-variant-associated devRegulons.
## Of the 71 regulons significantly associated with ASD rare-variant burden,
## 45.1% (32/71) are developmentally dynamic regulons (devRegulons). This panel
## shows how those 32 devRegulons distribute across developmental trend clusters.
## Horizontal bars, one per trend, ordered by proportion (largest at top).
## Reads panelA_ASD_devRegulon_trend_DATA.csv (same folder). Run:
##   Rscript panelA_ASD_devRegulon_trend.R
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr)
})

## ------------------------- CONFIG (edit freely) -------------------------
here <- tryCatch(dirname(normalizePath(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) getwd())
if (length(here) == 0 || is.na(here)) here <- getwd()
data_csv <- file.path(here, "panelA_ASD_devRegulon_trend_DATA.csv")
out_png  <- file.path(here, "figures", "panelA_ASD_devRegulon_trend_R.png")
out_pdf  <- file.path(here, "figures", "panelA_ASD_devRegulon_trend_R.pdf")

BAR_COL <- "#C9A227"                              # gold
note    <- "45.1% (32/71) of the rare-variant\nregulons belonged to devRegulons"
## ------------------------------------------------------------------------

d <- read_csv(data_csv, show_col_types = FALSE)
d <- d %>% arrange(proportion)                    # ggplot draws bottom->top
d$trend <- factor(d$trend, levels = d$trend)      # largest ends up at top
d$lab   <- sprintf("%.2f%%", 100 * d$proportion)

p <- ggplot(d, aes(x = proportion, y = trend)) +
  geom_col(fill = BAR_COL, width = 0.72) +
  geom_text(aes(label = lab), hjust = -0.15, size = 3.3, color = "#333333") +
  annotate("text", x = max(d$proportion) * 0.95,
           y = 1.6, label = note, hjust = 1, vjust = 0,
           size = 3.0, color = "#555555") +
  scale_x_continuous(limits = c(0, max(d$proportion) * 1.18),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Proportion of ASD rare-variant-associated devRegulons",
       x = "Proportion", y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.line = element_line(color = "#333333", linewidth = 0.4),
        axis.ticks = element_line(color = "#333333", linewidth = 0.4),
        axis.text.y = element_text(size = 10, color = "#222222"),
        plot.title = element_text(face = "bold", size = 12))

ggsave(out_png, p, width = 5.2, height = 3.4, dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 5.2, height = 3.4, bg = "white")
cat("saved:", out_png, "\n"); cat("saved:", out_pdf, "\n")
