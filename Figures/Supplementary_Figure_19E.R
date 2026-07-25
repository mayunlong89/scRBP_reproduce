
#-2026-2-25

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(grid)
})

# =========================
# 0) Input / Output
# =========================
infile  <- "01_Summary_6cardiovascular_traits_results_OpenTargets_SNP.csv"
out_pdf <- "l2g_fill_mlogp_by_disease_ALLSNP_ordered.pdf"

# 你想要的 facet 顺序（固定）
disease_want <- c("SBP","DBP","HT","PP","AF","CAD")

# sample size 图例 breaks（固定）
legend_N <- c(1e5, 2.5e5, 5e5, 1e6, 2e6)

# 颜色（和你原来一致）
pal_l2g <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)

# -------------------------
# helper: pick col
# -------------------------
pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# =========================
# 1) Read data
# =========================
df_raw <- read.csv(infile, check.names = FALSE, stringsAsFactors = FALSE)

p_col   <- pick_col(df_raw, c("p_value","pValue","Pvalue","PValue","pval","pval_nominal","P"))
n_col   <- pick_col(df_raw, c("sample_size","sampleSize","sampleSiz","N","n"))
rbp_col <- pick_col(df_raw, c("RBP","gene","Gene"))
snp_col <- pick_col(df_raw, c("SNP","rsid","rsID","variant","Variant"))
dis_col <- pick_col(df_raw, c("Diseases","Disease","reportedTrait","trait"))
l2g_col <- pick_col(df_raw, c("l2GScore","L2GScore","L2G","l2g","L2G_score","L2G.score"))

need <- c(p_col, n_col, rbp_col, snp_col, dis_col, l2g_col)
if (any(is.na(need))) {
  mapping <- c("p-value"=p_col, "sample_size"=n_col, "RBP"=rbp_col,
               "SNP"=snp_col, "Diseases"=dis_col, "L2GScore"=l2g_col)
  stop(
    "Missing required columns.\n\nFound columns:\n  ",
    paste(names(df_raw), collapse = ", "),
    "\n\nResolved mapping:\n",
    paste(names(mapping), "->", mapping, collapse = "\n")
  )
}

# =========================
# 2) Clean + compute -log10(P)
# =========================
df1 <- df_raw %>%
  transmute(
    RBP         = as.character(.data[[rbp_col]]),
    SNP         = as.character(.data[[snp_col]]),
    Diseases    = as.character(.data[[dis_col]]),
    sample_size = as.numeric(gsub(",", "", as.character(.data[[n_col]]))),
    p_value     = as.numeric(.data[[p_col]]),
    L2GScore    = as.numeric(.data[[l2g_col]])
  ) %>%
  mutate(
    p_value = pmax(p_value, .Machine$double.xmin),
    mlogp   = -log10(p_value),
    label   = paste0(RBP, " (", SNP, ")")
  ) %>%
  filter(
    is.finite(mlogp),
    is.finite(sample_size),
    !is.na(Diseases), Diseases != "",
    !is.na(label), label != ""
  )

# =========================
# 3) Force disease order (facet order)
# =========================
disease_levels <- disease_want[disease_want %in% unique(df1$Diseases)]
if (length(disease_levels) == 0) disease_levels <- unique(df1$Diseases)

df1 <- df1 %>%
  mutate(Diseases = factor(Diseases, levels = disease_levels))

# =========================
# 4) Sort y within each disease (descending -log10(P))
#    -> 你要的“反过来”效果：levels 用 rev(unique(label))
# =========================
df1 <- df1 %>%
  group_by(Diseases) %>%
  arrange(desc(mlogp), .by_group = TRUE) %>%
  mutate(y = factor(label, levels = rev(unique(label)))) %>%
  ungroup()

# =========================
# 5) Axis limits
# =========================
x_min <- max(0, floor(min(df1$mlogp, na.rm = TRUE)))
x_max <- ceiling(max(df1$mlogp, na.rm = TRUE) + 0.5)

# 参考线：Genome-wide significance p=5e-8
gws_line <- -log10(5e-8)   # 7.30103
# 如果你想 nominal p=0.05：nom_line <- -log10(0.05) # 1.30103

# =========================
# 6) Plot
# =========================
p <- ggplot(df1, aes(x = mlogp, y = y)) +
  facet_grid(Diseases ~ ., scales = "free_y", space = "free_y", switch = "y") +
  geom_vline(xintercept = gws_line, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_point(aes(fill = L2GScore, size = sample_size),
             shape = 21, color = "grey30", stroke = 0.3) +
  
  scale_fill_gradientn(
    colours = pal_l2g,
    name = "Open Targets\nL2G score",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    oob = scales::squish,
    guide = guide_colorbar(
      nbin = 256, raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4, "mm")
    )
  ) +
  
  scale_size_continuous(
    name   = "GWAS\n(sample sizes)",
    range  = c(1.8, 5.2),
    breaks = legend_N,
    labels = scales::label_number(
      scale_cut = scales::cut_short_scale()
    ),
    limits = c(min(legend_N, df1$sample_size, na.rm = TRUE),
               max(legend_N, df1$sample_size, na.rm = TRUE))
    # 如果跨度太大，可加：, trans = "sqrt"
  ) +
  
  scale_x_continuous(
    name = expression(-log[10](P)),
    limits = c(x_min, x_max),
    breaks = pretty(c(x_min, x_max), n = 6),
    expand = expansion(mult = c(0.01, 0.06))
  ) +
  
  labs(
    y = NULL,
    title = "Cardiovascular disease-associated variants in RBP master regulators are supported by Open Targets"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text.y.left = element_text(face = "bold", angle = 0),
    strip.placement = "outside",
    axis.text.y = element_text(size = 9, color = "grey20"),
    legend.position = "right",
    plot.margin = margin(5, 12, 5, 5)
  )

# =========================
# 7) Save
# =========================
n_rows <- nrow(df1)
n_dis  <- n_distinct(df1$Diseases)
h <- max(6, 0.22 * n_rows + 0.6 * n_dis)

ggsave(out_pdf, p, width = 9, height = h, useDingbats = FALSE)
p









