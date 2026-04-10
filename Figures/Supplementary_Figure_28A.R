



##_--------dotplot for regulon-opentargets enrichment
##_--------dotplot for regulon-opentargets enrichment
##_--------dotplot for regulon-opentargets enrichment
##_--------dotplot for regulon-opentargets enrichment
# ---- 数据 ----
txt <- "Regulon	NES	logp	padj
QKI_5UTR	1.307591657	5.000004343	0.000559994
DDX3X_5UTR	1.238209431	3.366535887	0.01203988
DDX21_5UTR	1.274739439	3.148745994	0.013253201
NIPBL_5UTR	1.255780212	2.567035052	0.025386413
TRA2A_5UTR	1.229844926	2.565435439	0.025386413
CELF1_3UTR	1.426940745	2.086416332	0.17544036
RBM47_Introns	1.238423435	2.037161662	0.183598164
TRA2A_3UTR	1.298532184	1.866146218	0.17544036
TAF15_3UTR	1.294268227	1.835949051	0.17544036
IGF2BP1_5UTR	1.424708484	1.478071846	0.155215445
DDX3X_3UTR	1.18704609	1.474182391	0.18910922
FAM120A_5UTR	1.256236336	1.296738167	0.200415774
YBX3_Introns	1.210306842	1.258218304	0.5518
PUM2_5UTR	1.127683633	1.234559325	0.200415774
WTAP_5UTR	1.250467089	1.202892334	0.200415774
DDX21_3UTR	1.141993404	1.078682099	0.32819231
RPRD2_5UTR	1.267089543	0.923869202	0.290128899
EIF4A1_5UTR	1.091410405	0.877886807	0.296729833
RYBP_3UTR	1.094149572	0.620517529	0.51800091
HNRNPH2_5UTR	1.130069427	0.560225305	0.468543335
YTHDF3_5UTR	1.089278448	0.558924323	0.468543335
ZNF148_Introns	1.047549317	0.515320216	0.614464295
PABPC1_Introns	1.033660783	0.453683819	0.670126632
HNRNPC_Introns	0.992828913	0.271212429	0.715745274
ILF3_3UTR	0.8350837	0.110167746	0.946277254
XRN2_Introns	0.914503258	0.07780326	0.926769111"

#df <- read.csv(text = txt, header = TRUE, sep = "\t", check.names = FALSE)




########------------ADHD
df <- read.csv("./02_opentargts_ADHD_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon

labeled_regulons <- read.csv("02_opentargts_ADHD_enriched_dotplot_all_sig_label_regulon_list.csv")

keep_regs <- labeled_regulons$Regulons


# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)







##-----------AN
df <- read.csv("./02_opentargts_AN_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon
labeled_regulons <- read.csv("02_opentargts_AN_enriched_dotplot_all_sig_label_regulon_list.csv")

keep_regs <- labeled_regulons$Regulons

# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)




########------------SCZ
df <- read.csv("./02_opentargts_SCZ_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)

# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon
labeled_regulons <- read.csv("02_opentargts_SCZ_enriched_dotplot_all_sig_label_regulon_list.csv")

keep_regs <- labeled_regulons$Regulons


# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)





########------------ASD
df <- read.csv("./02_opentargts_ASD_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon
labeled_regulons <- read.csv("02_opentargts_ASD_enriched_dotplot_all_sig_label_regulon_list.csv")

keep_regs <- labeled_regulons$Regulons



# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)





########-----------BIP
df <- read.csv("./02_opentargts_BIP_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon
labeled_regulons <- read.csv("02_opentargts_BIP_enriched_dotplot_all_sig_label_regulon_list.csv")

keep_regs <- labeled_regulons$Regulons



# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)





########-----------MDD
df <- read.csv("./02_opentargts_MDD_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon
labeled_regulons <- read.csv("02_opentargts_MDD_enriched_dotplot_all_sig_label_regulon_list.csv")

keep_regs <- labeled_regulons$Regulons




# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)







########-----------OCD
df <- read.csv("./02_opentargts_OCD_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon
labeled_regulons <- read.csv("02_opentargts_OCD_enriched_dotplot_all_sig_label_regulon_list.csv")

keep_regs <- labeled_regulons$Regulons


# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)







########----------TS
df <- read.csv("./02_opentargts_TS_enriched_dotplot_all.csv")

df$logp <- -log10(df$pval)
df$NES <- abs(df$NES)

head(df)
# 调色板：padj 越小越深
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)



# 标注：只给 FDR<0.05 的点加标签（按需修改阈值）
#df$label <- ifelse(df$logp > 1.47418239, df$Regulon, NA)
# 只标注指定 regulon
labeled_regulons <- read.csv("02_opentargts_TS_enriched_dotplot_all_sig_label_regulon_list.csv")

keep_regs <- labeled_regulons$Regulons

# 可选：去掉潜在空格/不可见字符，避免匹配不上
df$Regulon <- trimws(df$Regulon)

df$label <- ifelse(df$Regulon %in% keep_regs, df$Regulon, NA)

# ---- 作图 ----
library(ggplot2)
library(ggrepel)

p <- ggplot(df, aes(x = logp, y = NES, color = NES)) +
  geom_point(size = 2) +
  # 反转调色：padj 越小越深（用 rev(pal) 即可）
  scale_color_gradientn(colors = pal, name = "NES") +
  geom_text_repel(aes(label = label),
                  max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25,
                  min.segment.length = 0, seed = 42, size = 2.5) +
  labs(x = expression(-log[10](italic(P))),
       y = "NES",
       title = "Regulon enrichment (OpenTargets)") +
  theme_classic(base_size = 10)

p
# ggsave("regulon_scatter.png", p, width = 5, height = 6, dpi = 300)



