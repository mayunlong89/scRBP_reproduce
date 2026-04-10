

##--------proportion of significant RBP-regulons across diseases

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







#2026-02-11
setwd("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/03_scRBP_trs_fetal_brain_15times/05_Developmental_devRegulons")


suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ===== 1) 示例数据（替换成你的真实比例即可）=====
df_wide <- tribble(
  ~disease, ~devRegulons, ~non_devRegulons,
  "SCZ",    72.6,         77.0,
  "ADHD",   32.5,         25.0,
  "ASD",    28.2,         31.0,
  "BIP",    27.4,         37.0,
  "AN",     23.9,         20.0,
  "TS",     20.5,         14.0,
  "MDD",    12.8,          9.0,
  "OCD",     9.4,          6.0
)

disease_order <- c("SCZ","ADHD","ASD","BIP","AN","TS","MDD","OCD")

df_long <- df_wide %>%
  mutate(disease = factor(disease, levels = disease_order)) %>%
  pivot_longer(
    cols = c(devRegulons, non_devRegulons),
    names_to = "group",
    values_to = "prop_pct"
  ) %>%
  mutate(
    group = recode(group,
                   devRegulons = "devRegulons",
                   non_devRegulons = "non-devRegulons"),
    group = factor(group, levels = c("non-devRegulons", "devRegulons"))
  )

# ===== 2) 颜色（你可以换成自己 paper 的颜色）=====
pal_group <- c(
  "non-devRegulons" = "#F6C9C9",  # 浅粉
  "devRegulons"     = "#E07B7B"   # 深粉
)

# ===== 3) 作图：boxplot + 点(带黑边框) + 配对线 =====
p <- ggplot(df_long, aes(x = group, y = prop_pct)) +
  # 配对连线（灰色）
  geom_line(aes(group = disease), color = "grey60", linewidth = 0.8, alpha = 0.8) +
  
  # 箱线图（按组填充）
  geom_boxplot(
    aes(fill = group),
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.35,
    color = "black",
    linewidth = 0.9
  ) +
  
  # 点：shape=21 才能同时设置边框(color)和填充(fill)
  geom_point(
    aes(fill = group),
    shape = 21,
    size = 3.0,
    color = "black",
    stroke = 0.7
  ) +
  
  scale_fill_manual(values = pal_group) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, NA),
    expand = expansion(mult = c(0.05, 0.12))
  ) +
  labs(x = NULL, y = "Proportion of trait-relevant RBP regulons") +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13)
  )

p

ggsave("paired_boxplot_groupColor.pdf", p, width = 4.0, height = 4.2, useDingbats = FALSE)
ggsave("paired_boxplot_groupColor.tiff", p, width = 4.0, height = 4.2, dpi = 600, compression = "lzw")



##-----OR odds ratio for two groups
library(tidyverse)

df <- tribble(
  ~disease, ~dev, ~nondev,
  "SCZ", 72.6, 77.0,
  "ADHD", 32.5, 25.0,
  "ASD", 28.2, 31.0,
  "BIP", 27.4, 37.0,
  "AN", 23.9, 20.0,
  "TS", 20.5, 14.0,
  "MDD", 12.8, 9.0,
  "OCD", 9.4, 6.0
)

df_ratio <- df %>%
  mutate(
    PR = dev / nondev,
    logPR = log(PR)
  )

df_ratio

#calculate geometric mean proportion ratio
exp(mean(log(df_ratio$PR)))




