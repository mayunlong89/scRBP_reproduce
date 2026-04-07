#---2025-12-28
setwd("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/02_scRBP_trs_10blood_cell_traits_100times_v2/00-OpenTargets_monocyte_count")
##_--------dotplot for regulon-opentargets enrichment
##_--------dotplot for regulon-opentargets enrichment
##_--------dotplot for regulon-opentargets enrichment
##_--------dotplot for regulon-opentargets enrichment
##_-----In this version, we merge results from all 4 mRNA regions, then display enrichment and label significant results
########------------monocyte count all L2G genes
df <- read.csv("./01_monocyteCount_enrichment_GSEA_dotplot.csv")
df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)
head(df)
# Color palette: smaller padj = darker color
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)
# Labeling: only add labels for points with FDR < 0.05 (adjust threshold as needed)
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# Only label specified regulons
keep_regs <- c(
  "DDX21_3UTR",
  "DDX21_5UTR",
  "DDX3X_3UTR",
  "DDX3X_5UTR",
  "TRA2A_5UTR",
  "TRA2A_3UTR",
  "FAM120A_5UTR",
  "QKI_5UTR",
  "YTHDF3_5UTR",
  "CELF1_3UTR",
  "YBX3_Introns",
  "RBM47_Introns",
  "RYBP_3UTR",
  "AKAP1_5UTR"
)
# Optional: trim potential whitespace/invisible characters to avoid matching issues
df$Regulon <- trimws(df$Regulon)
df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)
# ---- Plot ----
library(ggplot2)
library(ggrepel)
p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # Reverse color scale: smaller padj = darker (use rev(pal))
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)
p
