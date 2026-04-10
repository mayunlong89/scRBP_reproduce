

##--------Supplementary Figure S20h---proportion of significant RBP-regulons across diseases

library(tidyverse)
library(scales)


# sum_long <- read.csv("trend_group_disease_sig_summary_long.csv", stringsAsFactors = FALSE)

disease_order <- c("SCZ","BIP","ASD","ADHD","AN","TS","MDD","OCD")

df_bar <- sum_long %>%
  mutate(grp = as.character(grp),
         disease = as.character(disease)) %>%
  filter(grp %in% c("devRegulons", "non-devRegulons")) %>%
  mutate(
    disease = factor(disease, levels = disease_order),
    grp = factor(grp, levels = c("devRegulons", "non-devRegulons")),
    label_prop = sprintf("%.1f%%", 100 * Proportion),
    label_nt = sprintf("%d/%d", Sig_N, Total)
  )


p <- ggplot(df_bar, aes(x = Proportion, y = disease, fill = grp)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_text(
    aes(label = label_prop),
    position = position_dodge(width = 0.75),
    hjust = -0.15, size = 3.8
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, max(df_bar$Proportion, na.rm = TRUE) * 1.15)
  ) +
  scale_fill_manual(
    values = c("devRegulons" = "#E76F51", "non-devRegulons" = "#F4A261"),
    name = NULL,
    labels = c("devRegulons", "non-devRegulons")
  ) +
  labs(
    x = "Proportion of significant regulons",
    y = NULL,
    title = "Proportion of significant RBP-regulons across psychiatric disorders"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

p
ggsave("dev_vs_nondev_proportion_bar.pdf", p, width = 7.2, height = 3.6, useDingbats = FALSE)



ord <- df_bar %>%
  filter(grp == "devRegulons") %>%
  arrange(desc(Proportion)) %>%
  pull(disease) %>% as.character()

df_bar2 <- df_bar %>%
  mutate(disease = factor(as.character(disease), levels = ord))

p_sorted <- ggplot(df_bar2, aes(x = Proportion, y = disease, fill = grp)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_text(aes(label = label_prop),
            position = position_dodge(width = 0.75),
            hjust = -0.15, size = 3.8) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, max(df_bar2$Proportion, na.rm = TRUE) * 1.15)) +
  scale_fill_manual(values = c("devRegulons"="#E76F51","non-devRegulons"="#F4A261"), name=NULL) +
  labs(x="Proportion of significant regulons", y=NULL,
       title="Proportion of significant RBP-regulons across psychiatric disorders") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face="bold"))

p_sorted
