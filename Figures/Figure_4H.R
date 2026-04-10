
#2025-12-18



# =========================================================
# L2G plot by disease (ALL SNPs shown)
# x: GWAS Z-score (two-sided) = qnorm(1 - p/2)
# fill: Open Targets L2G score (0-1)
# size: sample size
# =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(grid)
})

df_raw <- read.csv(
  "01_Summary_8autoimmune_diseases_results_OpenTargets_SNP.csv",
  check.names = FALSE, stringsAsFactors = FALSE
)

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

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
    z = qnorm(1 - p_value/2),
    label = paste0(RBP, " (", SNP, ")")
  ) %>%
  filter(
    is.finite(z),
    is.finite(sample_size),
    !is.na(Diseases), Diseases != ""
  )

# facet 顺序：按 max(z) 从大到小（更显著的 disease 在上面）
disease_order <- df1 %>%
  group_by(Diseases) %>%
  summarise(max_z = max(z, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_z)) %>%
  pull(Diseases)

df1 <- df1 %>%
  mutate(Diseases = factor(Diseases, levels = disease_order)) %>%
  group_by(Diseases) %>%
  arrange(z, .by_group = TRUE) %>%
  mutate(y = factor(label, levels = unique(label))) %>%
  ungroup()

pal_l2g <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)

# 自动扩展上限，保证 rs11256593 这种极显著点也能显示
x_min <- max(2, floor(min(df1$z, na.rm = TRUE)))
x_max <- ceiling(max(df1$z, na.rm = TRUE) + 0.5)

p <- ggplot(df1, aes(x = z, y = y)) +
  facet_grid(Diseases ~ ., scales = "free_y", space = "free_y") +
  geom_vline(xintercept = 3, linetype = "dashed", color = "grey70", linewidth = 0.4) +
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
  scale_size_continuous(name = "Sample size", range = c(1.8, 5.0), labels = comma) +
  scale_x_continuous(name = "GWAS Z-score (two-sided)",
                     limits = c(x_min, x_max),
                     breaks = pretty(c(x_min, x_max), n = 6),
                     expand = expansion(mult = c(0.01, 0.06))) +
  labs(y = NULL,
       title = "Autoimmune-associated variants in RBP master regulators are supported by Open Targets") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 9, color = "grey20"),
    legend.position = "right",
    plot.margin = margin(5, 12, 5, 5)
  )

n_rows <- nrow(df1)
n_dis  <- n_distinct(df1$Diseases)
h <- max(6, 0.22 * n_rows + 0.6 * n_dis)

ggsave("l2g_fill_z_by_disease_ALLSNP.pdf", p, width = 9, height = h, useDingbats = FALSE)
p
