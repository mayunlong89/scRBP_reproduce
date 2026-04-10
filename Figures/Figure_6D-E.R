

suppressPackageStartupMessages({
  library(tidyverse)
  library(binom)
  library(patchwork)
  library(scales)
})

# -------------------------
# Input (replace with read.csv if you want)
# -------------------------
df <- tribble(
  ~Group,          ~Total, ~Sig,
  "devRegulons",     237,   117,
  "non-devRegulons", 339,   100,
  "Trend 1",          26,    18,
  "Trend 2",          21,    20,
  "Trend 3",          24,    23,
  "Trend 4",          11,     0,
  "Trend 5",          26,     8,
  "Trend 6",          44,     1,
  "Trend 7",          16,    12,
  "Trend 8",          15,     6,
  "Trend 9",          23,    19,
  "Trend 10",         31,    10
)

# -------------------------
# Precompute proportion + Wilson CI
# -------------------------
df <- df %>%
  mutate(nonsig = Total - Sig,
         prop = Sig / Total)

ci_all <- binom.confint(x = df$Sig, n = df$Total, methods = "wilson")
df <- df %>%
  mutate(ci_low = ci_all$lower,
         ci_high = ci_all$upper,
         label = sprintf("%d/%d (%.1f%%)", Sig, Total, 100*prop))

# -------------------------
# TOP: dev vs non-dev
# -------------------------
df_top <- df %>%
  filter(Group %in% c("devRegulons", "non-devRegulons")) %>%
  mutate(
    Group = factor(Group, levels = c("devRegulons", "non-devRegulons")),
    Type = Group
  )

mat <- matrix(c(df_top$Sig[df_top$Group=="devRegulons"],
                df_top$nonsig[df_top$Group=="devRegulons"],
                df_top$Sig[df_top$Group=="non-devRegulons"],
                df_top$nonsig[df_top$Group=="non-devRegulons"]),
              nrow = 2, byrow = TRUE)
ft <- fisher.test(mat, alternative = "two.sided")
or  <- unname(ft$estimate)
p   <- ft$p.value
lab_top <- sprintf("Fisher's exact test (two-sided): OR = %.2f, P = %.2e", or, p)

bg_prop <- df_top$prop[df_top$Group=="non-devRegulons"]

# a clean two-color palette (color-blind friendly-ish)
col_dev <- "#1B9E77"   # teal/green
col_non <- "#7570B3"   # purple

p1 <- ggplot(df_top, aes(x = prop, y = Group, color = Type)) +
  geom_vline(xintercept = bg_prop, linetype = "dashed", linewidth = 0.6, color = "grey65") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.18, linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 4.0) +
  geom_text(aes(label = label), x = 1.02, hjust = 0, size = 3.4, color = "grey20") +
  scale_color_manual(values = c("devRegulons" = col_dev, "non-devRegulons" = col_non), guide = "none") +
  scale_x_continuous(limits = c(0, 1.22), labels = percent_format(accuracy = 1),
                     breaks = seq(0, 1, by = 0.2)) +
  labs(
    x = NULL, y = NULL,
    title = "Enrichment of significant regulons",
    subtitle = lab_top
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.margin = margin(8, 60, 6, 8)  # extra right margin for labels
  )

# -------------------------
# BOTTOM: Trends forest plot
# option A: color by Up/Down/Trans (recommended)
#   Down: Trend1-4, Trans: Trend5-7, Up: Trend8-10 (based on your figure)
# -------------------------
df_trend <- df %>%
  filter(grepl("^Trend", Group)) %>%
  mutate(
    trend_id = as.integer(gsub("Trend\\s+", "", Group)),
    Class = case_when(
      trend_id %in% 1:4 ~ "Down",
      trend_id %in% 5:7 ~ "Trans",
      trend_id %in% 8:10 ~ "Up",
      TRUE ~ "Other"
    ),
    Group = fct_reorder(Group, prop, .desc = TRUE)
  )

col_trend <- c(
  "Down"  = "#1B9E77",  # green
  "Trans" = "#984EA3",  # purple
  "Up"    = "#E41A1C"   # red
)

p2 <- ggplot(df_trend, aes(x = prop, y = Group, color = Class)) +
  geom_vline(xintercept = bg_prop, linetype = "dashed", linewidth = 0.6, color = "grey65") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.18, linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 3.4) +
  geom_text(aes(label = label), x = 1.02, hjust = 0, size = 3.2, color = "grey20") +
  scale_color_manual(values = col_trend, name = NULL) +
  scale_x_continuous(limits = c(0, 1.22), labels = percent_format(accuracy = 1),
                     breaks = seq(0, 1, by = 0.2)) +
  labs(
    x = "Proportion of significant regulons (95% CI)",
    y = NULL,
    title = "Developmental trend clusters"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    plot.margin = margin(6, 60, 8, 8)
  )

# -------------------------
# Combine
# -------------------------
(p1 / p2) + plot_layout(heights = c(1, 2))




