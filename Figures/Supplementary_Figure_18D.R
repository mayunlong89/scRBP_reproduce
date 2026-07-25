
#2026-02-23

setwd("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/05_scRBP_trs_Liver_nine_metabolic_traits")
df <- read.csv("01_heatmap_9metabolic_traits_Pvalue_Barplot.csv", check.names = FALSE)

library(dplyr)
library(tidyr)

# phenotype列
pheno_cols <- c("AKP_P","ALT_P","HDL_P","LDL_P",
                "SHBG_P","TBIL_P","TC_P","TG_P","TST_P")

# 转long format
df_long <- df %>%
  pivot_longer(cols = all_of(pheno_cols),
               names_to = "Phenotype",
               values_to = "Pvalue")

sig_total <- df_long %>%
  group_by(Celltype, Regulons) %>%
  summarise(sig_any = any(Pvalue < 0.05), .groups="drop") %>%
  group_by(Celltype) %>%
  summarise(
    n_significant_regulons = sum(sig_any),
    total_regulons = n(),
    proportion = n_significant_regulons / total_regulons
  ) %>%
  arrange(desc(n_significant_regulons))

sig_total




sig_by_trait <- df_long %>%
  filter(Pvalue < 0.05) %>%
  group_by(Celltype, Phenotype) %>%
  summarise(n_regulons = n_distinct(Regulons), .groups="drop") %>%
  arrange(Celltype, Phenotype)

sig_by_trait


heatmap_df <- sig_by_trait %>%
  pivot_wider(names_from = Phenotype,
              values_from = n_regulons,
              values_fill = 0)

heatmap_df


sig_ratio <- df_long %>%
  group_by(Celltype) %>%
  summarise(
    total_tests = n(),
    sig_tests = sum(Pvalue < 0.05),
    ratio = sig_tests / total_tests
  )

sig_ratio

write.csv(sig_ratio,file="01_heatmap_9metabolic_traits_Pvalue_Barplot_proportion.csv",quote=F, row.names = F)





library(dplyr)
library(ggplot2)

# 如果你的 sig_ratio 列叫 proportion / prop / ratio 之一，统一一下：
# sig_ratio <- sig_ratio %>% rename(ratio = proportion)

plot_df <- sig_ratio %>%
  mutate(Celltype = as.character(Celltype)) %>%
  arrange(desc(ratio)) %>%
  mutate(Celltype = factor(Celltype, levels = Celltype),
         Celltype_color = Celltype)

p_lollipop <- ggplot(plot_df, aes(x = Celltype, y = ratio)) +
  geom_segment(aes(x = Celltype, xend = Celltype, y = 0, yend = ratio),
               linewidth = 0.7, color = "grey70") +
  geom_point(aes(color = Celltype_color), size = 4) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Proportion",
    title = "Proportion of trait-relevant RBP-regulons"
  )

p_lollipop

ggsave("sig_ratio_lollipop.pdf", p_lollipop, width = 6.8, height = 3.2)
ggsave("sig_ratio_lollipop.tiff", p_lollipop, width = 6.8, height = 3.2, dpi = 300, compression = "lzw")





library(dplyr)
library(tidyr)
library(ggplot2)

alpha <- 0.05
pheno_cols <- c("AKP_P","ALT_P","HDL_P","LDL_P","SHBG_P","TBIL_P","TC_P","TG_P","TST_P")

# 1) long format
df_long <- df %>%
  pivot_longer(cols = all_of(pheno_cols),
               names_to = "Trait",
               values_to = "Pvalue")

# 2) 显著数目：每个 celltype 在所有 traits 中显著的 (regulon, trait) 对数
#    如果你想“同一个 regulon 在多个 traits 显著也重复计数”，用这个：
sig_count <- df_long %>%
  group_by(Celltype) %>%
  summarise(
    n_sig = sum(Pvalue < alpha),
    .groups = "drop"
  )

# ---- 如果你想“每个 regulon 只算一次（跨 traits 任意一个显著就算）”，用这个替代上面那段：
# sig_count <- df_long %>%
#   group_by(Celltype, Regulons) %>%
#   summarise(sig_any = any(Pvalue < alpha), .groups="drop") %>%
#   group_by(Celltype) %>%
#   summarise(n_sig = sum(sig_any), .groups="drop")

# 3) 排序（大到小），并让最大值在最上面
plot_df <- sig_count %>%
  arrange(desc(n_sig)) %>%
  mutate(Celltype = factor(Celltype, levels = rev(Celltype)))  # rev: 让最大在最上

# 4) 画竖起来的 lollipop（y=celltype, x=n_sig）
p_lollipop_h <- ggplot(plot_df, aes(x = n_sig, y = Celltype)) +
  geom_segment(aes(x = 0, xend = n_sig, y = Celltype, yend = Celltype),
               linewidth = 0.8, color = "grey70") +
  geom_point(size = 3.8) +
  theme_classic(base_size = 12) +
  labs(
    x = paste0("# significant tests (p < ", alpha, ")"),
    y = NULL,
    title = "Number of significant RBP-regulon associations by cell type"
  )

p_lollipop_h

ggsave("lollipop_sigCount_by_celltype.pdf", p_lollipop_h, width = 6.5, height = 4.8)
ggsave("lollipop_sigCount_by_celltype.tiff", p_lollipop_h, width = 6.5, height = 4.8, dpi = 300, compression = "lzw")







library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

# df <- read.csv("your_file.csv", check.names = FALSE, stringsAsFactors = FALSE)

pheno_cols <- c("AKP_P","ALT_P","HDL_P","LDL_P","SHBG_P","TBIL_P","TC_P","TG_P","TST_P")
alpha <- 0.05

df_long <- df %>%
  pivot_longer(cols = all_of(pheno_cols),
               names_to = "Trait",
               values_to = "Pvalue")


# 每个 (Celltype, Trait) 的显著 regulon 数
sig_by_trait <- df_long %>%
  group_by(Celltype, Trait) %>%
  summarise(
    n_regulons = n_distinct(Regulons),
    n_sig_regulons = n_distinct(Regulons[Pvalue < alpha]),
    prop_sig = n_sig_regulons / n_regulons,
    .groups = "drop"
  )

# （可选）按总显著数给 celltype 排序，让图更好看
cell_order <- sig_by_trait %>%
  group_by(Celltype) %>%
  summarise(total_sig = sum(n_sig_regulons), .groups="drop") %>%
  arrange(desc(total_sig)) %>%
  pull(Celltype)

sig_by_trait <- sig_by_trait %>%
  mutate(Celltype = factor(Celltype, levels = cell_order),
         Trait = factor(Trait, levels = pheno_cols))




p_bubble <- ggplot(sig_by_trait, aes(x = Trait, y = Celltype)) +
  geom_point(aes(size = n_sig_regulons, color = prop_sig), alpha = 0.85) +
  scale_size_continuous(name = paste0("# significant regulons\n(p < ", alpha, ")")) +
  scale_color_continuous(name = "Proportion significant") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(x = NULL, y = NULL)

p_bubble

# 画完之后直接加这一行
p_bubble_rev <- p_bubble +
  scale_y_discrete(limits = rev(levels(sig_by_trait$Celltype)))

p_bubble_rev

ggsave("bubble_trait_by_celltype.pdf", p_bubble, width = 8.5, height = 6.5)
ggsave("bubble_trait_by_celltype.tiff", p_bubble, width = 8.5, height = 6.5, dpi = 300, compression = "lzw")



