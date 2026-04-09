
### 1) R plot pheatmap

```R
##----FINAL VERISON
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

res_df <- read.csv(
  "66_sharedProteins_Jaccard_RBPregions_vs_TF__10kb_up_and_down_tss.minRankAgg_top1500.csv",
  check.names = FALSE
)

res <- as.matrix(res_df[,-1])
rownames(res) <- res_df[,1]

# =========================
# 强化版 olive palette
# =========================

olive_pal_strong <- colorRampPalette(c(
  "#F8FBEF",   # 更浅
  "#DCE5B5",
  "#B9C87A",
  "#8FA34E",
  "#4F6228"    # 更深
))(100)

# 关键：固定范围，增强对比
breaks_strong <- seq(0.02, 0.10, length.out = 101)

pdf("66_sharedProteins_Jaccard_RBPregions_vs_TF__10kb_up_and_down_tss.minRankAgg_top1500_olive_strong.pdf",
    width = 8,
    height = 10)

pheatmap(
  res,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,
  color = olive_pal_strong,
  breaks = breaks_strong,
  border_color = NA,
  main = "Jaccard similarity (min-rank aggregation, top 1500 per protein)"
)

dev.off()
```

### 2)R boxplot with Wilcoxlon test
```R

res_df <- read.csv("66_sharedProteins_Jaccard_RBPregions_vs_TF__10kb_up_and_down_tss.minRankAgg_top1500",
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
```

### 3) meme2images
```bash
#2026-03-02

#YBX1
meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF11117 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/YBX1_MTF11117.out_logo_eps


meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF10799 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/YBX1_MTF10799.out_logo_eps



meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF11628 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/YBX1_MTF11628.out_logo_eps

###############################################################

meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF10885 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/CPEB1_MTF10885.out_logo_eps


meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF11848 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/CPEB1_MTF11848.out_logo_eps


###############################################################

meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF25073 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/ZNF777_MTF25073.out_logo_eps


meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF25074 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/ZNF777_MTF25074.out_logo_eps



meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF25080 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/ZNF777_MTF25080.out_logo_eps




meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF25100 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/ZNF777_MTF25100.out_logo_eps



###############################################################


meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF23745 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/ZNF449_MTF23745.out_logo_eps


meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF23747 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/ZNF449_MTF23747.out_logo_eps


meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF23785 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/ZNF449_MTF23785.out_logo_eps




###############################################################

meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF15815 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/SP1_MTF15815.out_logo_eps


meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF15825 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/SP1_MTF15825.out_logo_eps


meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF15830 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/SP1_MTF15830.out_logo_eps





###############################################################

meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF02390 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/ZNF800_MTF02390.out_logo_eps


meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF07366 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/ZNF800_MTF07366.out_logo_eps


meme_dir=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2
output=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings/03_616RBP_vs_1390TFs_66common_proteins
meme2images -eps -motif MTF02401 $meme_dir/04_merged_all_25044motifs_all_databases.meme $output/ZNF800_MTF02401.out_logo_eps
###############################################################

```





