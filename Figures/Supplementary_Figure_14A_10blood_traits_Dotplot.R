#2025-12-30
# =========================================================
# RBP-SNP -> Trait bubble dotplot
#  - x: Trait
#  - y: RBP:SNP_rs
#  - size: -log10(pValue)
#  - color: L2GScore
#  - outline: significant (optional)
# =========================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})
# =========================
# 0) Input
# =========================
infile <- "01_9_blood_cell_traits_SNPs_RBPs.csv"
df0 <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)
# Parse pValue: compatible with "3.2x10-5" / "3.2×10-5" / "3.2e-5"
parse_p <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- NA_character_
  
  x <- gsub("×", "x", x)
  x <- gsub("x10\\^?([+-]?)", "e\\1", x)
  x <- gsub("X10\\^?([+-]?)", "e\\1", x)
  x <- gsub("x10", "e", x, ignore.case = TRUE)
  
  suppressWarnings(as.numeric(x))
}
# Sample size column name (modify as needed)
sample_col <- "sampleSize"  # e.g. "N" / "n" / "sample_size"
if (!sample_col %in% colnames(df0)) {
  stop(sprintf("Cannot find sample size column '%s'. Current column names:\n%s",
               sample_col, paste(colnames(df0), collapse = ", ")))
}
# =========================
# 1) Clean
# =========================
df <- df0 %>%
  mutate(
    pValue_num = parse_p(pValue),
    L2GScore   = suppressWarnings(as.numeric(L2GScore)),
    N          = suppressWarnings(as.numeric(.data[[sample_col]]))
  ) %>%
  filter(
    !is.na(RBP), !is.na(SNP_rs), !is.na(Trait),
    is.finite(pValue_num), pValue_num > 0,
    is.finite(L2GScore),
    is.finite(N)
  ) %>%
  mutate(
    logp    = -log10(pValue_num),
    rbp_snp = paste0(RBP, " (", SNP_rs, ")")
  )
# =========================
# 2) For each Trait x RBP, keep the most significant SNP (smallest p first; ties broken by larger N)
# =========================
plot_df <- df %>%
  group_by(Trait, RBP) %>%
  arrange(pValue_num, desc(N)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    rbp_snp  = paste0(RBP, " (", SNP_rs, ")"),
    size_val = sqrt(logp)  # Compress extreme values to avoid oversized circles
  )
# =========================
# 3) Ordering: y = rbp_snp (sorted by overall most significant p, then reversed for display)
# =========================
y_levels <- plot_df %>%
  group_by(rbp_snp) %>%
  summarise(p_min = min(pValue_num), .groups = "drop") %>%
  arrange(p_min) %>%
  pull(rbp_snp)
plot_df <- plot_df %>%
  mutate(
    rbp_snp = factor(rbp_snp, levels = rev(y_levels)),  # <<< Reversed (as desired)
    Trait   = factor(Trait, levels = unique(Trait))
  )
# =========================
# 4) Plot
# =========================
p <- ggplot(plot_df, aes(x = Trait, y = rbp_snp)) +
  geom_point(
    aes(size = logp, fill = L2GScore),   # <<< Use size_val instead of logp
    shape = 21, colour = "grey25", stroke = 0.25
  ) +
  scale_size_continuous(
    name  = expression(-log[10](p)),
    range = c(1.2, 4.2)   # <<< Control max circle size (reduce 4.2 to 3 for smaller circles)
  ) +
  scale_fill_gradient(
    name   = "L2GScore",
    low    = "#F3F3F3",
    high   = "#F87C56",
    limits = c(0, 1),
    oob    = squish
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 6, angle = 60, hjust = 1, vjust = 1),  # <<< Reduce x label size
    axis.text.y = element_text(size = 4),                                    # <<< Adjust if y labels are too crowded
    plot.margin = margin(5.5, 12, 5.5, 5.5)
  )
p
# ggsave("RBP_SNP_byTrait_dotplot_x_is_Trait_y_is_RBP_SNP.pdf", p, width = 6.5, height = 12, useDingbats = FALSE)
