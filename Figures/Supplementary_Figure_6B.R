
#--2026-03-05
library(tidyverse)

# 1) read data
df <- readr::read_tsv("ALL_9aging_TRAITS.E_statistic_D.matrix_N_vs_ODC_window_size.tsv")

# 2) long format
long <- df %>%
  pivot_longer(
    cols = -window_size,
    names_to = "pheno_region",
    values_to = "Escore"
  ) %>%
  mutate(
    region = str_extract(pheno_region, "(?i)(5UTR|3UTR|CDS|INTRON[S]?)") %>%
      toupper() %>%
      str_replace("^INTRON[S]?$", "INTRON"),
    region = factor(region, levels = c("5UTR", "CDS", "INTRON", "3UTR")),
    window_size = factor(
      window_size,
      levels = c("0kb", "5kb", "10kb", "25kb", "50kb", "100kb")
    )
  )

# 3) single color
col_single <- "#C47BA6"

# 4) plot
p <- ggplot(long, aes(x = window_size, y = Escore)) +
  geom_boxplot(
    width = 0.62,
    alpha = 0.28,
    fill = col_single,
    color = "grey30",
    outlier.shape = NA,
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(shape = region),
    width = 0.13,
    size = 2.7,
    alpha = 0.9,
    color = col_single,
    stroke = 0.1
  ) +
  scale_shape_manual(
    values = c(
      "5UTR" = 16,   # circle
      "CDS" = 17,    # triangle
      "INTRON" = 15, # square
      "3UTR" = 18    # diamond
    ),
    labels = c("5’UTR", "CDS", "Intron", "3’UTR")
  ) +
  labs(
    title = "Assessment across 9 aging-related traits",
    x = "Window size (SNP-to-Gene)",
    y = "E-statistics"
  ) +
  guides(shape = guide_legend(title = NULL)) +
  theme_classic(base_size = 13) +
  theme(
    axis.ticks.length = unit(3, "mm"),
    axis.text.x = element_text(margin = margin(t = 6)),
    axis.text.y = element_text(margin = margin(r = 6))
  ) +
  coord_cartesian(
    ylim = c(min(long$Escore, na.rm = TRUE) - 0.01,
             max(long$Escore, na.rm = TRUE) + 0.015),
    clip = "off"
  )

print(p)


