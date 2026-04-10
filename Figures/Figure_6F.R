
#2026-01-26

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

setwd("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/03_scRBP_trs_fetal_brain_15times/01_Final_scRBP_TRS_fetal_brain_50times")
 
# =========================
# 0) Input / Output
# =========================
infile  <- "03_devRegulon_disease_heatmap.csv"
out_wide <- "03_devRegulon_disease_heatmap_trend_group_disease_sig_summary_wide.csv"
out_long <- "03_devRegulon_disease_heatmap_trend_group_disease_sig_summary_long.csv"
out_pdf  <- "03_devRegulon_disease_heatmap_trend_group_disease_heatmap_dot.pdf"

# =========================
# 1) Read data
# =========================
df <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)

# =========================
# 2) Identify disease P-value columns
#    Your columns look like: SCZ_P, ASD_P, BIP_P, TS_P, OCD_P, MDD_P, AN_P, ADHD_P
# =========================
p_cols <- names(df)[grepl("_P$", names(df))]
if (length(p_cols) == 0) stop("No *_P columns found. Please check column names.")

diseases <- sub("_P$", "", p_cols)

# Optional: set disease order (edit as you like)
disease_order <- c("SCZ","BIP","ASD","ADHD","AN","TS","MDD","OCD")
disease_order <- disease_order[disease_order %in% diseases]  # keep existing ones
if (length(disease_order) == 0) disease_order <- diseases

# =========================
# 3) Build the groups you want to summarize
#    - devRegulons / non-devRegulons (use df$Group)
#    - Trend 1–10 (use df$Trends)
#
#    Key idea:
#    We create TWO copies and then bind:
#      (A) group_type = "Group", grp = Group
#      (B) group_type = "Trend", grp = Trends
# =========================
stopifnot(all(c("Group","Trends") %in% names(df)))

df_group <- df %>%
  mutate(group_type = "Group",
         grp = Group)

df_trend <- df %>%
  mutate(group_type = "Trend",
         grp = Trends) %>%
  filter(!is.na(grp), grp != "")

df2 <- bind_rows(df_group, df_trend)

# Keep only the groups you care about (optional but recommended)
group_levels <- c("devRegulons","non-devRegulons", paste0("Trend ", 1:10))
df2 <- df2 %>% filter(grp %in% group_levels)

# =========================
# 4) Long format: one row per (regulon, disease)
# =========================
df_long <- df2 %>%
  pivot_longer(cols = all_of(p_cols), names_to = "disease_p", values_to = "pval") %>%
  mutate(
    disease = sub("_P$", "", disease_p),
    sig = ifelse(!is.na(pval) & pval < 0.05, 1L, 0L)
  )

# =========================
# 5) Summarize: Total / Sig_N / Proportion for each (grp, disease)
# =========================
sum_long <- df_long %>%
  group_by(grp, disease) %>%
  summarise(
    Total = n(),
    Sig_N = sum(sig, na.rm = TRUE),
    Proportion = Sig_N / Total,
    .groups = "drop"
  ) %>%
  mutate(
    grp = factor(grp, levels = group_levels),
    disease = factor(disease, levels = disease_order)
  )

# Save long summary (useful for plotting / debugging)
write.csv(sum_long, out_long, row.names = FALSE)

# =========================
# 6) Make the "Fig2-style" wide table:
#    columns: SCZ_N, SCZ_proportion, BIP_N, BIP_proportion, ...
#    rows: devRegulons, non-devRegulons, Trend 1..10
# =========================
sum_wide <- sum_long %>%
  select(grp, disease, Sig_N, Proportion) %>%
  pivot_wider(
    names_from = disease,
    values_from = c(Sig_N, Proportion),
    names_glue = "{disease}_{.value}"
  ) %>%
  arrange(grp)

# Reorder columns to: disease_N then disease_proportion (like your screenshot)
# Here we rename to *_N and *_proportion
names(sum_wide) <- gsub("_Sig_N$", "_N", names(sum_wide))
names(sum_wide) <- gsub("_Proportion$", "_proportion", names(sum_wide))

# Put grp as first column named "Group"
sum_wide <- sum_wide %>%
  rename(Group = grp)

write.csv(sum_wide, out_wide, row.names = FALSE)

# =========================
# 7) Plot heatmap (tile = proportion; dot size = Sig_N)
#    Similar to your dot heatmap style.
# =========================
p_heat <- ggplot(sum_long, aes(x = disease, y = grp)) +
  geom_tile(aes(fill = Proportion), color = "white", linewidth = 0.4) +
  geom_point(aes(size = Sig_N), color = "black") +
  scale_size_continuous(range = c(0.5, 6)) +
  scale_fill_gradient(low = "white", high = "#B30000", limits = c(0, 1)) +
  labs(
    x = NULL, y = NULL,
    fill = "Proportion\n(sig/total)",
    size = "Significant\nregulons (N)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  )

p_heat 
ggsave(out_pdf, p_heat, width = 10.5, height = 4.2, useDingbats = FALSE)
p_heat



 

# ---- 你想保留的行（只画 Trends） ----
trend_levels <- paste0("Trend ", 1:10)
disease_levels <- c("SCZ","BIP","ASD","ADHD","AN","TS","MDD","OCD")  # 按你的顺序
disease_levels <- disease_levels[disease_levels %in% unique(as.character(sum_long$disease))]

# 1) 只保留 Trend，并补齐所有 Trend × Disease 组合（缺失填 0）
sum_trend <- sum_long %>%
  mutate(
    grp = as.character(grp),
    disease = as.character(disease)
  ) %>%
  filter(grp %in% trend_levels) %>%
  right_join(
    expand_grid(grp = trend_levels, disease = disease_levels),
    by = c("grp","disease")
  ) %>%
  mutate(
    Total = replace_na(Total, 0),
    Sig_N = replace_na(Sig_N, 0),
    Proportion = replace_na(Proportion, 0),
    grp = factor(grp, levels = trend_levels),
    disease = factor(disease, levels = disease_levels)
  )

# 先定义你希望图例显示的 N 值
legend_N <- c(0, 10, 20, 50)

p_heat <- ggplot(sum_trend, aes(x = disease, y = grp)) +
  geom_tile(aes(fill = Proportion), color = "white", linewidth = 0.6) +
  geom_point(aes(size = Sig_N), color = "black") +
  scale_size_continuous(
    name   = "Significant\nregulons (N)",
    range  = c(0.8, 9),
    breaks = legend_N,
    limits = c(0, max(legend_N, max(sum_trend$Sig_N, na.rm = TRUE)))
    # 如果你用了 sqrt，就加上这一行，并且图例会自动一致
    # , trans = "sqrt"
  ) +
  scale_fill_gradient(
    name   = "Proportion\n(sig/total)",
    low    = "white",
    high   = "#B30000",
    limits = c(0, 1)
  ) +
  guides(
    size = guide_legend(override.aes = list(alpha = 1)), # 可选：确保图例点不透明
    fill = guide_colorbar()
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold")
  )

p_heat



