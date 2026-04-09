#2026-03-18

# =========================================================
# Plot permutation null distribution of median Jaccard similarity
# =========================================================

# =========================================================
# ggplot2 version
# =========================================================

library(ggplot2)

# -------------------------
# 1. Read input file
# -------------------------
infile <- "66_sharedProteins_Jaccard_RBPregions_vs_TF__10kb_up_and_down_tss.minRankAgg_top1500.globalMedian_permutation_n1000.csv"
df <- read.csv(infile, header = TRUE, stringsAsFactors = FALSE)

# -------------------------
# 2. Extract observed and permutation values
# -------------------------
observed_value <- df$median_js[df$type == "observed"]
perm_df <- df[df$type == "permutation", ]

if (length(observed_value) != 1) {
  stop("Error: observed value should have exactly one row.")
}

if (nrow(perm_df) == 0) {
  stop("Error: no permutation values found.")
}

# permutation p-value
p_perm <- (sum(perm_df$median_js >= observed_value) + 1) / (nrow(perm_df) + 1)

# -------------------------
# 3. Plot
# -------------------------
p <- ggplot(perm_df, aes(x = median_js)) +
  geom_histogram(
    bins = 35,
    fill = "#CFE8CF",     # 浅绿色
    color = "gray50",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = observed_value,
    color = "red",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = observed_value,
    y = Inf,
    label = paste0("Observed = ", sprintf("%.3f", observed_value)),
    color = "red",
    vjust = 1.5,
    hjust = -0.05,
    size = 6
  ) +
  annotate(
    "text",
    x = min(perm_df$median_js),
    y = Inf,
    label = paste0("Permutation p = ", signif(p_perm, 3)),
    hjust = 0,
    vjust = 3.2,
    size = 5
  ) +
  labs(
    x = "Median Jaccard similarity",
    y = "Frequency"
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    plot.margin = margin(15, 25, 15, 15)
  )

print(p)


