
res_df <- read.csv("66_sharedProteins_Jaccard_RBPregions_vs_TF__10kb_up_and_down_tss.minRankAgg_top1500.csv",
                   check.names = FALSE)

library(tidyverse)

df_long <- res_df %>%
  pivot_longer(
    cols = -Protein,
    names_to = "Region",
    values_to = "Jaccard"
  ) %>%
  mutate(
    Region = str_remove(Region, "RBP_"),
    Region = str_remove(Region, "_vs_TF"),
    Region = factor(Region, levels = c("UTR5","UTR3","CDS","Intron"))
  )

# 转宽格式方便 paired test
df_wide <- df_long %>%
  pivot_wider(names_from = Region, values_from = Jaccard)

# paired tests
p_5_vs_3  <- wilcox.test(df_wide$UTR5, df_wide$UTR3, paired = TRUE)$p.value
p_5_vs_C  <- wilcox.test(df_wide$UTR5, df_wide$CDS,  paired = TRUE)$p.value
p_5_vs_I  <- wilcox.test(df_wide$UTR5, df_wide$Intron, paired = TRUE)$p.value

pvals <- c(p_5_vs_3, p_5_vs_C, p_5_vs_I)
pvals_adj <- p.adjust(pvals, method = "BH")

names(pvals_adj) <- c("UTR5 vs UTR3","UTR5 vs CDS","UTR5 vs Intron")

print(pvals_adj)

library(ggplot2)

ggplot(df_long, aes(x = Region, y = Jaccard, fill = Region)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  scale_fill_manual(values = c(
    UTR5   = "#6F8F2F",
    UTR3   = "#A9C05B",
    CDS    = "#C9D48A",
    Intron = "#E5EBCF"
  )) +
  theme_classic(base_size = 14) +
  labs(
    y = "Jaccard similarity (TF vs RBP targets)",
    x = ""
  ) +
  theme(legend.position = "none")
