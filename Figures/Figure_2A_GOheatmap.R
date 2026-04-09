# --- 2025-11-09
# GO enrichment treemap for 616 RBPs using semantic similarity-based term reduction

# ============================================================
# Load required packages
# ============================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(rrvgo)
library(dplyr)
library(ggplot2)
library(treemapify)

# ============================================================
# Step 1: GO enrichment analysis (BP)
# ============================================================
# Read the list of 616 RBP gene symbols
rbp_symbols <- read.table("616_RBPs_list.csv", header = FALSE, stringsAsFactors = FALSE)$V1
n_rbp <- length(rbp_symbols)

# Run GO biological process enrichment
eg <- enrichGO(
  gene          = rbp_symbols,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# Extract significant terms and compute -log10(FDR) as the enrichment score
df <- as.data.frame(eg) %>%
  filter(p.adjust <= 0.05) %>%
  mutate(
    score = -log10(p.adjust),
    Count = as.integer(Count)
  )

# ============================================================
# Step 2: Compute GO term semantic similarity matrix
# ============================================================
semData <- godata("org.Hs.eg.db", ont = "BP", keytype = "SYMBOL")
simMat  <- calculateSimMatrix(
  df$ID,
  orgdb  = "org.Hs.eg.db",
  ont    = "BP",
  method = "Rel"
)

# ============================================================
# Step 3: Reduce redundant GO terms into representative clusters
# ============================================================
# Higher threshold = stricter merging; typical range: 0.5–0.8
reduced <- reduceSimMatrix(
  simMat,
  scores    = setNames(df$score, df$ID),
  threshold = 0.7,
  orgdb     = "org.Hs.eg.db"
)

# ============================================================
# Step 4: Merge enrichment statistics back into reduced terms
# ============================================================
# Detect column names (compatible across rrvgo versions)
id_col      <- intersect(c("term_ID", "ID", "term.id", "go", "term"), names(reduced))[1]
name_col    <- intersect(c("term_name", "name", "Description"), names(reduced))[1]
cluster_col <- intersect(c("cluster", "parent"), names(reduced))[1]
stopifnot(!is.na(id_col), !is.na(name_col), !is.na(cluster_col))

# Join Count and p.adjust from enrichGO results
plot_df <- reduced %>%
  rename(term_id = !!id_col, term_name = !!name_col, cluster = !!cluster_col) %>%
  left_join(
    df %>% select(ID, Count, p.adjust) %>% rename(term_id = ID),
    by = "term_id"
  ) %>%
  mutate(
    ratio = ifelse(!is.na(Count), Count / n_rbp, NA_real_),
    label_show = if_else(
      !is.na(Count),
      sprintf("%s\nn=%d | %.1f%%", term_name, Count, 100 * ratio),
      term_name
    )
  ) %>%
  arrange(desc(score))

# ============================================================
# Step 5: Treemap visualization
# ============================================================
# Area: enrichment score; Fill: enrichment score; Label: term + gene count + ratio
p <- ggplot(
  plot_df,
  aes(area = score, fill = score, label = label_show, subgroup = cluster)
) +
  geom_treemap(show.legend = TRUE, start = "topleft") +
  geom_treemap_subgroup_border(color = "white", size = 0.8) +
  geom_treemap_text(
    colour   = "white",
    place    = "centre",
    grow     = TRUE,
    reflow   = TRUE,
    family   = "sans",
    min.size = 2
  ) +
  scale_fill_gradient(
    name = expression(-log[10] ~ "FDR"),
    low  = "#E7F2A4",
    high = "#7F815A"
  ) +
  guides(fill = guide_colorbar(barheight = unit(4, "cm"))) +
  theme_void()

print(p)
# ggsave("RBP_GO_treemap.svg", p, width = 10, height = 7)
