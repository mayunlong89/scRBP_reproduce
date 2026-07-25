
#2026-01-31

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(stringr)
})

# =========================
# 0) Input / Output
# =========================
infile  <- "03_devRegulon_disease_heatmap.csv"
out_pdf <- "SuppFig_trend_by_disease_proportion_plus_regulon_table_dedup.pdf"
out_csv <- "SuppFig_trend_by_disease_proportion_plus_regulon_table_dedup_long.csv"

TOP_N <- 12
ONLY_SIG <- TRUE
P_CUTOFF <- 0.05

disease_levels <- c("SCZ","BIP","ASD","ADHD","AN","TS","MDD","OCD")
trend_levels   <- paste0("Trend ", 1:10)

fill_low  <- "white"
fill_high <- "#B30000"
text_size <- 2.6

# 是否让 Proportion 也按“去重后的 regulon”计算？
# FALSE: Proportion按原始ID计数（与你主图一致）
# TRUE : Proportion按去重后的Regulon计数（更符合“regulon层面”）
PROP_DEDUP_BY_REGULON <- FALSE

# =========================
# 1) Read data
# =========================
df <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot("Trends" %in% colnames(df))

p_cols <- names(df)[grepl("_P$", names(df))]
if (length(p_cols) == 0) stop("No *_P columns found. Please check column names.")
diseases <- sub("_P$", "", p_cols)

disease_levels <- disease_levels[disease_levels %in% diseases]
if (length(disease_levels) == 0) disease_levels <- diseases

# 你要显示的 regulon 名称列：优先 Regulon，否则用 ID
name_col <- if ("Regulon" %in% names(df)) "Regulon" else if ("ID" %in% names(df)) "ID" else NA
if (is.na(name_col)) stop("Neither 'Regulon' nor 'ID' column found for labeling.")

# =========================
# 2) Long format: (ID, regulon_name, trend, disease, pval)
# =========================
df_long <- df %>%
  mutate(
    trend = as.character(Trends),
    ID = if ("ID" %in% names(.)) as.character(ID) else NA_character_,
    regulon_name = as.character(.data[[name_col]])
  ) %>%
  filter(trend %in% trend_levels) %>%
  pivot_longer(cols = all_of(p_cols), names_to = "disease_p", values_to = "pval") %>%
  mutate(
    disease = sub("_P$", "", disease_p),
    disease = factor(disease, levels = disease_levels),
    trend   = factor(trend, levels = trend_levels),
    sig     = ifelse(!is.na(pval) & pval < P_CUTOFF, 1L, 0L)
  ) %>%
  select(ID, regulon_name, trend, disease, pval, sig)

# =========================
# 3) Proportion (Sig/Total) per cell
# =========================
if (!PROP_DEDUP_BY_REGULON) {
  # (A) 和主图一致：按原始行（ID层面）计算
  sum_cell <- df_long %>%
    group_by(trend, disease) %>%
    summarise(
      Total = n(),
      Sig_N = sum(sig, na.rm = TRUE),
      Proportion = ifelse(Total > 0, Sig_N / Total, 0),
      .groups = "drop"
    )
} else {
  # (B) 按去重后的 regulon 计算：每个 (trend,disease,regulon) 取最小p作为代表
  sum_cell <- df_long %>%
    filter(!is.na(pval), !is.na(regulon_name), regulon_name != "") %>%
    group_by(trend, disease, regulon_name) %>%
    summarise(p_min = min(pval, na.rm = TRUE), .groups = "drop") %>%
    mutate(sig = as.integer(p_min < P_CUTOFF)) %>%
    group_by(trend, disease) %>%
    summarise(
      Total = n(),
      Sig_N = sum(sig, na.rm = TRUE),
      Proportion = ifelse(Total > 0, Sig_N / Total, 0),
      .groups = "drop"
    )
}

# =========================
# 4) Text per cell (去重后列出 regulons)
#    去重逻辑：同一格子内同一regulon取 p_min=min(pval)
# =========================
df_for_text <- df_long %>%
  filter(!is.na(pval), !is.na(regulon_name), regulon_name != "") %>%
  # 先去重：每个 (trend,disease,regulon) 取最小p
  group_by(trend, disease, regulon_name) %>%
  summarise(p_min = min(pval, na.rm = TRUE), .groups = "drop") %>%
  { if (ONLY_SIG) filter(., p_min < P_CUTOFF) else . }

text_cell <- df_for_text %>%
  group_by(trend, disease) %>%
  summarise(
    n_reg = n(),
    txt = {
      ord <- order(p_min, na.last = TRUE)
      regs <- regulon_name[ord]
      show_regs <- head(regs, TOP_N)
      more_n <- max(length(regs) - TOP_N, 0)
      txt0 <- paste(show_regs, collapse = "\n")
      if (more_n > 0) paste0(txt0, "\n", "+", more_n, " more") else txt0
    },
    .groups = "drop"
  )

# =========================
# 5) Complete grid & merge
# =========================
plot_df <- expand_grid(
  trend   = factor(trend_levels, levels = trend_levels),
  disease = factor(disease_levels, levels = disease_levels)
) %>%
  left_join(sum_cell,  by = c("trend","disease")) %>%
  left_join(text_cell, by = c("trend","disease")) %>%
  mutate(
    Total      = replace_na(Total, 0),
    Sig_N      = replace_na(Sig_N, 0),
    Proportion = replace_na(Proportion, 0),
    n_reg      = replace_na(n_reg, 0),
    txt        = replace_na(txt, "")
  )

write.csv(plot_df, out_csv, row.names = FALSE)

# =========================
# 6) Plot
# =========================
p_table <- ggplot(plot_df, aes(x = disease, y = trend)) +
  geom_tile(aes(fill = Proportion), color = "white", linewidth = 0.6) +
  geom_text(
    aes(label = txt),
    size = text_size,
    lineheight = 0.95,
    vjust = 0.5,
    hjust = 0.5
  ) +
  scale_fill_gradient(
    name   = "Proportion\n(sig/total)",
    low    = fill_low,
    high   = fill_high,
    limits = c(0, 1)
  ) +
  scale_y_discrete(limits = rev(trend_levels)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(8, 8, 8, 8)
  )

print(p_table)
ggsave(out_pdf, p_table, width = 12.5, height = 7.8, useDingbats = FALSE)
