# --- 2025-07-07
# Gene-level analysis

# Load required packages
library(tidyverse)
library(dplyr)
library(ggplot2)

# Read data
df <- read.csv("02_TopN_sumary_AUC_thresholds_genes.csv", header = TRUE)

# Compute the incremental contribution of each rank tier
df_stack <- df %>%
  mutate(
    top1  = top1rank,
    top2  = top2rank  - top1rank,
    top3  = top3rank  - top2rank,
    top4  = top4rank  - top3rank,
    top5  = top5rank  - top4rank,
    top10 = top10rank - top5rank,
    top15 = top15rank - top10rank,
    top20 = top20rank - top15rank
  ) %>%
  dplyr::select(group, top1, top2, top3, top4, top5, top10, top15, top20)

# Convert to long format
df_long <- df_stack %>%
  pivot_longer(cols = -group, names_to = "rank_group", values_to = "value")

# Set stacking order from finest granularity to broadest
df_long$rank_group <- factor(df_long$rank_group,
                             levels = c("top1", "top2", "top3", "top4", "top5", "top10", "top15", "top20"))

# Preserve the original group order from the CSV
df_long$group <- factor(df_long$group, levels = df$group)

# Compute the label for the bar top (e.g., "X/223 RBP hits")
label_df <- df %>%
  mutate(total_detected = round(top20rank * 223)) %>%
  dplyr::select(group, total_detected)

# Define color palette for each rank tier
rank_colors <- c("top1"  = "#800000",
                 "top2"  = "#B22222",
                 "top3"  = "#FF4500",
                 "top4"  = "#FFA500",
                 "top5"  = "#FFD700",
                 "top10" = "#EEE8AA",
                 "top15" = "#FAFAD2",
                 "top20" = "#F9EFE5")

# Plot stacked bar chart
ggplot(df_long, aes(x = group, y = value, fill = rank_group)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  geom_text(data = label_df,
            aes(x = group, y = 1.02, label = paste0(total_detected, "/223")),
            inherit.aes = FALSE, size = 3.5) +
  scale_fill_manual(values = rank_colors, labels = names(rank_colors)) +
  labs(title = "Validation\nRank thresholds (%) for AUC calculation",
       y = "%RBP Recovery", x = "AUC threshold") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.title = element_blank())
