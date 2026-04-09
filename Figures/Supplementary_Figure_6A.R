# --- 2025-10-29
# Final version — E-score evaluation of MAGMA window sizes

library(tidyverse)

# 1) Read data and reshape to long format
df <- readr::read_csv("01_MAGMA_window_sizes_assessment_Escore.csv")

long <- df %>%
  pivot_longer(cols = -window_size, names_to = "pheno_region", values_to = "Escore") %>%
  # Extract and standardize region names
  mutate(
    region = str_extract(pheno_region, "(?i)(5UTR|3UTR|CDS|INTRON[S]?)") %>%
      toupper() %>%
      str_replace("^INTRON[S]?$", "INTRONS"),
    region = factor(region, levels = c("5UTR", "CDS", "INTRONS", "3UTR")),
    window_size = factor(window_size, levels = c("0kb", "5kb", "10kb", "25kb", "50kb", "100kb"))
  )

# 2) Single color (consistent with Figure 1)
col_single <- "#CC79A7"

# 3) Plot: boxplot + jittered points (shape by region; fixed single color)
p <- ggplot(long, aes(x = window_size, y = Escore)) +
  geom_boxplot(width = 0.62, alpha = 0.28,
               fill = col_single, color = "grey20", outlier.shape = NA) +
  geom_jitter(aes(shape = region),
              width = 0.14, size = 2.2, alpha = 0.9,
              color = col_single, stroke = 0) +
  # Shape mapping for each mRNA region
  scale_shape_manual(values = c("5UTR"    = 16,   # circle
                                "CDS"     = 17,   # triangle
                                "INTRONS" = 15,   # square
                                "3UTR"    = 18)) + # diamond
  labs(x = "Window size", y = "E-score") +
  guides(shape = guide_legend(title = NULL)) +
  theme_classic(base_size = 13) +
  theme(
    plot.title.position = "plot",
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6))
  )

print(p)
