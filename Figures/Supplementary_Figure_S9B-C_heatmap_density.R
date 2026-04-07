#2026-03-04

setwd("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/09_scRBP_infer_grn_adultbrain_fetal_brain_Lung")

##---------Fetal brain Nature Neuroscience paper---#################################################################################################
##---------Fetal brain Nature Neuroscience paper---#################################################################################################
##---------Fetal brain Nature Neuroscience paper---#################################################################################################
##---------Fetal brain Nature Neuroscience paper---#################################################################################################

#-------Part I---- Convert files into the long format we need
# ---- packages ----
library(tidyverse)

# ---- Parameters: fill in your four gmt file paths here ----
gmt_files <- c(
  "scRBP_prune_matched_results_min10genes_5UTR_fetal.gmt",
  "scRBP_prune_matched_results_min10genes_CDS_fetal.gmt",
  "scRBP_prune_matched_results_min10genes_Introns_fetal.gmt",
  "scRBP_prune_matched_results_min10genes_3UTR_fetal.gmt"
)

# Extract region from filename
extract_region <- function(path) {
  m <- stringr::str_match(basename(path), "(5UTR|CDS|Introns|3UTR)")
  ifelse(is.na(m[,2]), NA_character_, m[,2])
}

# Read one gmt -> long table: RBP, region, gene
read_gmt_long <- function(path) {
  region <- extract_region(path)
  lines  <- readr::read_lines(path, progress = FALSE)
  if (length(lines) == 0) return(tibble(RBP=character(), region=character(), gene=character()))

  parse_line <- function(line) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(f) < 3) return(NULL)
    name1 <- f[1]
    desc  <- f[2]
    # Get RBP name: prefer desc with "_regulon" stripped; fall back to name1 if missing
    rbp   <- if (!is.na(desc)) gsub("_regulon$", "", desc) else name1
    rbp   <- if (is.na(rbp) || rbp == "") name1 else rbp
    genes <- f[-c(1,2)] |> trimws()
    genes <- genes[genes != "" & !is.na(genes)]
    if (length(genes) == 0) return(NULL)
    tibble(RBP = rbp, region = region, gene = genes)
  }
  purrr::map_dfr(lines, parse_line)
}

# Merge four files
df_long <- purrr::map_dfr(gmt_files, read_gmt_long) |>
  mutate(region = factor(region, levels = c("5UTR","CDS","Introns","3UTR"))) |>
  distinct()  # deduplicate (when the same RBP-region-gene appears multiple times)

# Save
readr::write_tsv(df_long, "rbp_targets_long_fetal.tsv")

# Quick check
message("RBPs: ", n_distinct(df_long$RBP),
        " | rows: ", nrow(df_long),
        " | per-region counts:\n",
        paste(capture.output(df_long |> count(region)), collapse = "\n"))



#--------Part II---- heatmap plotting----

# ---- pkgs ----
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ComplexHeatmap)
library(circlize)

# ---- Input: long table ----
df <- read_tsv("rbp_targets_long_fetal.tsv", show_col_types = FALSE)
regions <- c("5UTR","CDS","Introns","3UTR")
df <- df %>% filter(region %in% regions) %>%
  mutate(region = factor(region, levels = regions))

# ---- Gene set for each (RBP, region) ----
sets <- df %>% group_by(RBP, region) %>%
  summarise(genes = list(unique(gene)), .groups = "drop")

rbps <- sets %>% pull(RBP) %>% unique() %>% sort()

# =======================
# Heatmap 1: counts (359 x 4)
# =======================
count_df <- sets %>%
  mutate(n = lengths(genes)) %>%
  dplyr::select(RBP, region, n) %>%
  tidyr::complete(RBP = rbps, region = factor(regions, levels = regions), fill = list(n = 0)) %>%
  dplyr::arrange(factor(RBP, levels = rbps))

# Convert to matrix (explicitly set rownames to avoid issues with as.matrix)
count_mat <- count_df %>%
  tidyr::pivot_wider(names_from = region, values_from = n, values_fill = 0) %>%
  as.data.frame()
rownames(count_mat) <- count_mat$RBP
count_mat <- as.matrix(count_mat[, regions, drop = FALSE])


library(pheatmap)

## ---------- Heatmap 1: counts (359 x 4) ----------
# Assumes count_mat already exists (rows = RBP, cols = c("5UTR","CDS","Introns","3UTR"))
pal_cnt <- colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(101)



# Ensure numeric matrix
count_num <- as.matrix(count_mat)
storage.mode(count_num) <- "double"

# --- Discrete breaks (use 0.5 offset to align integer intervals) ---
maxv <- max(count_num, na.rm = TRUE)

# 0, 5-20, 21-50, 51-100, 101-200, >=201
bk <- c(-0.5, 0.5, 20.5, 50.5, 100.5, 200.5, maxv + 0.5)
#bk <- c(-0.5, 0.5, 30.5, 50.5, 150.5, 250.5, maxv + 0.5)

# 6 discrete colors (replace with your preferred palette)
cols <- c("white",           # 0
          "#fee5d9",         # 5-20  (light)
          "#fcae91",         # 21-50
          "#fb6a4a",         # 51-100
          "#de2d26",         # 101-200
          "#a50f15")         # >=201 (dark)

# Legend breaks and labels (corresponding to the intervals above)
leg_brks <- c(0, 12, 35, 75, 150, min(250, maxv))
leg_labs <- c("0", "5-20", "21-50", "51-100", "101-200",
              ifelse(maxv >= 201, "\u2265201", paste0("\u2264", maxv)))

pheatmap(
  count_num,
  color = cols,
  breaks = bk,                 # key: discrete bin mapping
  na_col = "#eeeeee",
  border_color = "grey90",
  cluster_rows = TRUE,
  cluster_cols = F,        # column order fixed as 5UTR/CDS/Introns/3UTR
  show_rownames = TRUE,        # turn on to show row names
  fontsize_row = 5,
  cellwidth = 8, cellheight = 1.1, fontsize = 6,
  legend_breaks = leg_brks,
  legend_labels = leg_labs,
  main = "Regulon size (discrete bins)"
)




# ==============================
# Heatmap 2: pairwise Jaccard (359x6)
# ==============================
# Ensure sets is correct (each row = 1 RBP x region, genes is a list column)
# sets <- df %>% group_by(RBP, region) %>% summarise(genes = list(unique(gene)), .groups = "drop")

# --- Safe grouping approach: split into list by RBP ---
sets$RBP   <- as.character(sets$RBP)     # prevent issues caused by factor/encoding
sets$region<- factor(sets$region, levels = c("5UTR","CDS","Introns","3UTR"))
set_map <- split(sets[, c("region","genes")], sets$RBP)
length(set_map)          # should be approximately 359
head(names(set_map), 3)  # check the first few RBP names

# --- Jaccard helper function ---
jaccard <- function(a, b){
  a <- unique(a); b <- unique(b)
  u <- union(a,b); if (length(u)==0) return(NA_real_)
  length(intersect(a,b)) / length(u)
}

regions <- levels(sets$region)
pairs   <- combn(regions, 2, simplify = FALSE)
pair_names <- sapply(pairs, function(x) paste(x, collapse="|"))

# --- Compute Jaccard for 6 pairs per RBP ---
calc_one <- function(tbl){     # tbl is the data.frame(region, genes) for one RBP
  getg <- function(r){
    idx <- which(tbl$region == r)
    if (length(idx)==0) character(0) else tbl$genes[[idx[1]]]
  }
  sapply(pairs, function(p) jaccard(getg(p[1]), getg(p[2])))
}

jac_mat <- t(sapply(set_map, calc_one))
colnames(jac_mat) <- pair_names


library(pheatmap)

# Color scale and range: 0 ~ 0.30
MAX_COL <- 0.5
#MAX_COL <- 1
pal  <- colorRampPalette(c("#f7fbff", "#74add1", "#08306b"))(121)   # 120-level color scale
#pal  <- colorRampPalette(c("#F7F5CE", "#E77A77", "firebrick3"))(121)   # 120-level color scale
brks <- seq(0, MAX_COL, length.out = 122)                           # one more than number of colors

# (optional) clip display values to 0.30 before plotting
# jac_disp <- pmin(jac_num, MAX_COL)
jac_disp <- pmin(jac_mat, MAX_COL)  # clip values to MAX_COL to prevent exceeding color range
hr <- hclust(dist(jac_mat))   # hierarchical clustering of RBPs (rows) by Euclidean distance
hc <- hclust(dist(t(jac_mat)))  # hierarchical clustering of columns (pairwise Jaccard)

pheatmap(
  jac_disp,
  color = pal, breaks = brks, na_col = "#eeeeee",
  border_color = NA,cellwidth =8, cellheight =1.1,fontsize=6,
  cluster_rows = hr, cluster_cols = hc,          # use the clustering trees computed above
  show_rownames = T,
  frontsize_row = 1,
  legend_breaks = seq(0, MAX_COL, by = 0.05),    # legend tick marks (adjustable)
  legend_labels = sprintf("%.2f", seq(0, MAX_COL, by = 0.05)),
  main = "Pairwise Jaccard across regions (0-0.30)"
)



####--------4 domains---> region sizes distributions

library(readr)
library(dplyr)
library(ggplot2)

# 1) Read in and compute regulon size
regions <- c("5UTR","CDS","Introns","3UTR")

df <- read_tsv("rbp_targets_long_fetal.tsv", show_col_types = FALSE) |>
  distinct(RBP, region, gene) |>                               # deduplicate to prevent duplicates
  filter(region %in% regions) |>
  mutate(region = factor(region, levels = regions))

sizes <- df |>
  group_by(region, RBP) |>
  summarise(size = n_distinct(gene), .groups = "drop")         # number of genes per RBP x region

# 2) Basic statistics
summary_tab <- sizes |>
  group_by(region) |>
  summarise(
    n_regulons = n(),
    median = median(size),
    mean   = mean(size),
    IQR    = IQR(size),
    min    = min(size),
    max    = max(size),
    .groups = "drop"
  )
print(summary_tab)

# 3) Plotting (choose one)



###_------ Important ----density plot
library(dplyr)
library(ggplot2)

breaks <- log1p(c(5,10,15,20,30,50,100,200,500,1000))
# Compute median per region (original scale) and convert to log1p scale
meds <- sizes %>%
  group_by(region) %>%
  summarise(med_size = median(size, na.rm = TRUE), .groups = "drop") %>%
  mutate(med_log1p = log1p(med_size))

ggplot(sizes, aes(x = log1p(size), color = region, fill = region)) +
  geom_density(alpha = 0.25, adjust = 1.2, linewidth = 1) +
  geom_vline(
    data = meds,
    aes(xintercept = med_log1p, color = region),
    linetype = "dashed", linewidth = 0.8, alpha = 0.9
  ) +
  scale_x_continuous(
    name = "Regulon size (# genes)",
    breaks = breaks,
    labels = ~ scales::comma(round(expm1(.)))
  ) +
  labs(y = "Density") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic()


## 3a. Violin + boxplot (log scale, suitable for long-tailed distributions)
ggplot(sizes, aes(region, size, fill = region)) +
  geom_violin(trim = FALSE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.size = 0.6, alpha = 0.7, color = "grey20") +
  scale_y_continuous(trans = "log1p") +
  labs(x = "Region", y = "Regulon size (# genes)", title = "Regulon size distribution by region") +
  theme_bw() + theme(legend.position = "none")


ggplot(sizes, aes(region, size, fill = region)) +
  geom_violin(trim = FALSE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.size = 0.6, alpha = 0.7, color = "grey20") +
  scale_y_continuous(trans = "log1p") +
  labs(x = "Region", y = "Regulon size (# genes)", title = "Regulon size distribution by region") +
  theme_classic() + theme(legend.position = "none")

ticks <- c(5, 10, 25, 50, 100, 300, 600, 900)
ggplot(sizes, aes(region, size, fill = region)) +
  geom_violin(trim = FALSE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.size = 0.6, alpha = 0.7, color = "grey20") +
  scale_y_continuous(
    trans = "log1p",
    breaks = ticks, labels = ticks,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = "Region", y = "Regulon size (# genes)",
       title = "Regulon size distribution by region") +
  theme_classic() + theme(legend.position = "none")


## 3b. Faceted histogram (one panel per region)
ggplot(sizes, aes(size)) +
  geom_histogram(bins = 30, color = "white") +
  facet_wrap(~ region, scales = "free_y") +
  scale_x_continuous(trans = "log1p") +
  labs(x = "Regulon size (# genes)", y = "Count") +
  theme_bw()

## 3c. ECDF (cumulative distribution curve, useful for comparing overall size differences)
ggplot(sizes, aes(size, color = region)) +
  stat_ecdf(size = 1) +
  scale_x_continuous(trans = "log1p") +
  labs(x = "Regulon size (# genes)", y = "ECDF", title = "ECDF of regulon sizes") +
  theme_classic()

# Compute quantiles per region
sizes %>%
  group_by(region) %>%
  summarise(n = n(),
            q25 = quantile(size, .25),
            median = median(size),
            q75 = quantile(size, .75),
            p90 = quantile(size, .90),
            max = max(size))

# Add median dashed lines to ECDF plot

qs <- sizes %>% group_by(region) %>% summarise(med = median(size))
ggplot(sizes, aes(size, color = region)) +
  stat_ecdf(size = 1) +
  geom_vline(data = qs, aes(xintercept = med, color = region),
             linetype = "dashed", alpha = .7) +
  scale_x_continuous(trans = "log1p") +
  labs(x = "Regulon size (# genes)", y = "ECDF",
       title = "ECDF of regulon sizes (dashed = median)") +
  theme_bw()


qs <- sizes %>% group_by(region) %>% summarise(med = median(size))
ticks <- c(5, 10, 25, 50, 100,200, 300, 500, 1000,1500)
ggplot(sizes, aes(size, color = region)) +
  stat_ecdf(size = 1) +
  geom_vline(data = qs, aes(xintercept = med, color = region),
             linetype = "dashed", alpha = .6) +
  scale_x_continuous(
    trans = "log1p",
    breaks = ticks, labels = ticks,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = "Regulon size (# genes)", y = "ECDF",
       title = "ECDF of regulon sizes (dashed = median)") +
  theme_classic()

# 4) (Optional) Non-parametric test: are regulon sizes different across regions?
kruskal.test(size ~ region, data = sizes)
pairwise.wilcox.test(sizes$size, sizes$region, p.adjust.method = "BH")



##---------Palmer et al. PNAS paper---#################################################################################################
##---------Palmer et al. PNAS paper---#################################################################################################
##---------Palmer et al. PNAS paper---#################################################################################################
##---------Palmer et al. PNAS paper---#################################################################################################

#-------Part I---- Convert files into the long format we need
# ---- packages ----
library(tidyverse)

# ---- Parameters: fill in your four gmt file paths here ---- #100K cells based on 50 times
gmt_files <- c(
  "scRBP_prune_matched_results_min10genes.symbol_default_palmer_adult_100Kcells_5UTR.gmt",
  "scRBP_prune_matched_results_min10genes.symbol_default_palmer_adult_100Kcells_CDS.gmt",
  "scRBP_prune_matched_results_min10genes.symbol_default_palmer_adult_100Kcells_Introns.gmt",
  "scRBP_prune_matched_results_min10genes.symbol_default_palmer_adult_100Kcells_3UTR.gmt"
)

# Extract region from filename
extract_region <- function(path) {
  m <- stringr::str_match(basename(path), "(5UTR|CDS|Introns|3UTR)")
  ifelse(is.na(m[,2]), NA_character_, m[,2])
}

# Read one gmt -> long table: RBP, region, gene
read_gmt_long <- function(path) {
  region <- extract_region(path)
  lines  <- readr::read_lines(path, progress = FALSE)
  if (length(lines) == 0) return(tibble(RBP=character(), region=character(), gene=character()))

  parse_line <- function(line) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(f) < 3) return(NULL)
    name1 <- f[1]
    desc  <- f[2]
    # Get RBP name: prefer desc with "_regulon" stripped; fall back to name1 if missing
    rbp   <- if (!is.na(desc)) gsub("_regulon$", "", desc) else name1
    rbp   <- if (is.na(rbp) || rbp == "") name1 else rbp
    genes <- f[-c(1,2)] |> trimws()
    genes <- genes[genes != "" & !is.na(genes)]
    if (length(genes) == 0) return(NULL)
    tibble(RBP = rbp, region = region, gene = genes)
  }
  purrr::map_dfr(lines, parse_line)
}

# Merge four files
df_long <- purrr::map_dfr(gmt_files, read_gmt_long) |>
  mutate(region = factor(region, levels = c("5UTR","CDS","Introns","3UTR"))) |>
  distinct()  # deduplicate (when the same RBP-region-gene appears multiple times)

# Save
readr::write_tsv(df_long, "rbp_targets_long_palmer_100K.tsv")

# Quick check
message("RBPs: ", n_distinct(df_long$RBP),
        " | rows: ", nrow(df_long),
        " | per-region counts:\n",
        paste(capture.output(df_long |> count(region)), collapse = "\n"))



#--------Part II---- heatmap plotting----

# ---- pkgs ----
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ComplexHeatmap)
library(circlize)

# ---- Input: long table ----
df <- read_tsv("rbp_targets_long_palmer_100K.tsv", show_col_types = FALSE)
regions <- c("5UTR","CDS","Introns","3UTR")
df <- df %>% filter(region %in% regions) %>%
  mutate(region = factor(region, levels = regions))

# ---- Gene set for each (RBP, region) ----
sets <- df %>% group_by(RBP, region) %>%
  summarise(genes = list(unique(gene)), .groups = "drop")

rbps <- sets %>% pull(RBP) %>% unique() %>% sort()

# =======================
# Heatmap 1: counts (359 x 4)
# =======================
count_df <- sets %>%
  mutate(n = lengths(genes)) %>%
  dplyr::select(RBP, region, n) %>%
  tidyr::complete(RBP = rbps, region = factor(regions, levels = regions), fill = list(n = 0)) %>%
  dplyr::arrange(factor(RBP, levels = rbps))

# Convert to matrix (explicitly set rownames to avoid issues with as.matrix)
count_mat <- count_df %>%
  tidyr::pivot_wider(names_from = region, values_from = n, values_fill = 0) %>%
  as.data.frame()
rownames(count_mat) <- count_mat$RBP
count_mat <- as.matrix(count_mat[, regions, drop = FALSE])


library(pheatmap)

## ---------- Heatmap 1: counts (359 x 4) ----------
# Assumes count_mat already exists (rows = RBP, cols = c("5UTR","CDS","Introns","3UTR"))
pal_cnt <- colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(101)



# Ensure numeric matrix
count_num <- as.matrix(count_mat)
storage.mode(count_num) <- "double"

# --- Discrete breaks (use 0.5 offset to align integer intervals) ---
maxv <- max(count_num, na.rm = TRUE)

# 0, 5-20, 21-50, 51-100, 101-200, >=201
bk <- c(-0.5, 0.5, 20.5, 50.5, 100.5, 200.5, maxv + 0.5)
#bk <- c(-0.5, 0.5, 30.5, 50.5, 150.5, 250.5, maxv + 0.5)

# 6 discrete colors (replace with your preferred palette)
cols <- c("white",           # 0
          "#fee5d9",         # 5-20  (light)
          "#fcae91",         # 21-50
          "#fb6a4a",         # 51-100
          "#de2d26",         # 101-200
          "#a50f15")         # >=201 (dark)

# Legend breaks and labels (corresponding to the intervals above)
leg_brks <- c(0, 12, 35, 75, 150, min(250, maxv))
leg_labs <- c("0", "5-20", "21-50", "51-100", "101-200",
              ifelse(maxv >= 201, "\u2265201", paste0("\u2264", maxv)))

pheatmap(
  count_num,
  color = cols,
  breaks = bk,                 # key: discrete bin mapping
  na_col = "#eeeeee",
  border_color = "grey90",
  cluster_rows = T,
  cluster_cols = F,        # column order fixed as 5UTR/CDS/Introns/3UTR
  show_rownames = TRUE,        # turn on to show row names
  fontsize_row = 5,
  cellwidth = 8, cellheight = 1.1, fontsize = 6,
  legend_breaks = leg_brks,
  legend_labels = leg_labs,
  main = "Regulon size (discrete bins)"
)




# ==============================
# Heatmap 2: pairwise Jaccard (359x6)
# ==============================
# Ensure sets is correct (each row = 1 RBP x region, genes is a list column)
# sets <- df %>% group_by(RBP, region) %>% summarise(genes = list(unique(gene)), .groups = "drop")

# --- Safe grouping approach: split into list by RBP ---
sets$RBP   <- as.character(sets$RBP)     # prevent issues caused by factor/encoding
sets$region<- factor(sets$region, levels = c("5UTR","CDS","Introns","3UTR"))
set_map <- split(sets[, c("region","genes")], sets$RBP)
length(set_map)          # should be approximately 359
head(names(set_map), 3)  # check the first few RBP names

# --- Jaccard helper function ---
jaccard <- function(a, b){
  a <- unique(a); b <- unique(b)
  u <- union(a,b); if (length(u)==0) return(NA_real_)
  length(intersect(a,b)) / length(u)
}

regions <- levels(sets$region)
pairs   <- combn(regions, 2, simplify = FALSE)
pair_names <- sapply(pairs, function(x) paste(x, collapse="|"))

# --- Compute Jaccard for 6 pairs per RBP ---
calc_one <- function(tbl){     # tbl is the data.frame(region, genes) for one RBP
  getg <- function(r){
    idx <- which(tbl$region == r)
    if (length(idx)==0) character(0) else tbl$genes[[idx[1]]]
  }
  sapply(pairs, function(p) jaccard(getg(p[1]), getg(p[2])))
}

jac_mat <- t(sapply(set_map, calc_one))
colnames(jac_mat) <- pair_names


library(pheatmap)

# Color scale and range: 0 ~ 0.30
MAX_COL <- 0.5
#MAX_COL <- 1
pal  <- colorRampPalette(c("#f7fbff", "#74add1", "#08306b"))(121)   # 120-level color scale
#pal  <- colorRampPalette(c("#F7F5CE", "#E77A77", "firebrick3"))(121)   # 120-level color scale
brks <- seq(0, MAX_COL, length.out = 122)                           # one more than number of colors

# (optional) clip display values to 0.30 before plotting
# jac_disp <- pmin(jac_num, MAX_COL)
jac_disp <- pmin(jac_mat, MAX_COL)  # clip values to MAX_COL to prevent exceeding color range
hr <- hclust(dist(jac_mat))   # hierarchical clustering of RBPs (rows) by Euclidean distance
hc <- hclust(dist(t(jac_mat)))  # hierarchical clustering of columns (pairwise Jaccard)

pheatmap(
  jac_disp,
  color = pal, breaks = brks, na_col = "#eeeeee",
  border_color = NA,cellwidth =8, cellheight =1.1,fontsize=6,
  cluster_rows = hr, cluster_cols = hc,          # use the clustering trees computed above
  show_rownames = T,
  frontsize_row = 1,
  legend_breaks = seq(0, MAX_COL, by = 0.05),    # legend tick marks (adjustable)
  legend_labels = sprintf("%.2f", seq(0, MAX_COL, by = 0.05)),
  main = "Pairwise Jaccard across regions (0-0.30)"
)



####--------4 domains---> region sizes distributions

library(readr)
library(dplyr)
library(ggplot2)

# 1) Read in and compute regulon size
regions <- c("5UTR","CDS","Introns","3UTR")

df <- read_tsv("rbp_targets_long_palmer_100K.tsv", show_col_types = FALSE) |>
  distinct(RBP, region, gene) |>                               # deduplicate to prevent duplicates
  filter(region %in% regions) |>
  mutate(region = factor(region, levels = regions))

sizes <- df |>
  group_by(region, RBP) |>
  summarise(size = n_distinct(gene), .groups = "drop")         # number of genes per RBP x region

# 2) Basic statistics
summary_tab <- sizes |>
  group_by(region) |>
  summarise(
    n_regulons = n(),
    median = median(size),
    mean   = mean(size),
    IQR    = IQR(size),
    min    = min(size),
    max    = max(size),
    .groups = "drop"
  )
print(summary_tab)

# 3) Plotting (choose one)


###_------ Important ----density plot
library(dplyr)
library(ggplot2)

breaks <- log1p(c(5,10,15,20,30,50,100,200,500,1000))
# Compute median per region (original scale) and convert to log1p scale
meds <- sizes %>%
  group_by(region) %>%
  summarise(med_size = median(size, na.rm = TRUE), .groups = "drop") %>%
  mutate(med_log1p = log1p(med_size))

ggplot(sizes, aes(x = log1p(size), color = region, fill = region)) +
  geom_density(alpha = 0.25, adjust = 1.2, linewidth = 1) +
  geom_vline(
    data = meds,
    aes(xintercept = med_log1p, color = region),
    linetype = "dashed", linewidth = 0.8, alpha = 0.9
  ) +
  scale_x_continuous(
    name = "Regulon size (# genes)",
    breaks = breaks,
    labels = ~ scales::comma(round(expm1(.)))
  ) +
  labs(y = "Density") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic()


## 3a. Violin + boxplot (log scale, suitable for long-tailed distributions)
ggplot(sizes, aes(region, size, fill = region)) +
  geom_violin(trim = FALSE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.size = 0.6, alpha = 0.7, color = "grey20") +
  scale_y_continuous(trans = "log1p") +
  labs(x = "Region", y = "Regulon size (# genes)", title = "Regulon size distribution by region") +
  theme_bw() + theme(legend.position = "none")


ggplot(sizes, aes(region, size, fill = region)) +
  geom=violin(trim = FALSE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.size = 0.6, alpha = 0.7, color = "grey20") +
  scale_y_continuous(trans = "log1p") +
  labs(x = "Region", y = "Regulon size (# genes)", title = "Regulon size distribution by region") +
  theme_classic() + theme(legend.position = "none")

ticks <- c(5, 10, 25, 50, 100, 300, 600, 900)
ggplot(sizes, aes(region, size, fill = region)) +
  geom_violin(trim = FALSE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.size = 0.6, alpha = 0.7, color = "grey20") +
  scale_y_continuous(
    trans = "log1p",
    breaks = ticks, labels = ticks,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = "Region", y = "Regulon size (# genes)",
       title = "Regulon size distribution by region") +
  theme_classic() + theme(legend.position = "none")


## 3b. Faceted histogram (one panel per region)
ggplot(sizes, aes(size)) +
  geom_histogram(bins = 30, color = "white") +
  facet_wrap(~ region, scales = "free_y") +
  scale_x_continuous(trans = "log1p") +
  labs(x = "Regulon size (# genes)", y = "Count") +
  theme_bw()

## 3c. ECDF (cumulative distribution curve, useful for comparing overall size differences)
ggplot(sizes, aes(size, color = region)) +
  stat_ecdf(size = 1) +
  scale_x_continuous(trans = "log1p") +
  labs(x = "Regulon size (# genes)", y = "ECDF", title = "ECDF of regulon sizes") +
  theme_classic()

# Compute quantiles per region
sizes %>%
  group_by(region) %>%
  summarise(n = n(),
            q25 = quantile(size, .25),
            median = median(size),
            q75 = quantile(size, .75),
            p90 = quantile(size, .90),
            max = max(size))

# Add median dashed lines to ECDF plot

qs <- sizes %>% group_by(region) %>% summarise(med = median(size))
ggplot(sizes, aes(size, color = region)) +
  stat_ecdf(size = 1) +
  geom_vline(data = qs, aes(xintercept = med, color = region),
             linetype = "dashed", alpha = .7) +
  scale_x_continuous(trans = "log1p") +
  labs(x = "Regulon size (# genes)", y = "ECDF",
       title = "ECDF of regulon sizes (dashed = median)") +
  theme_bw()


qs <- sizes %>% group_by(region) %>% summarise(med = median(size))
ticks <- c(5, 10, 25, 50, 100,200, 300, 500, 1000,1500, 2000, 2500)
ggplot(sizes, aes(size, color = region)) +
  stat_ecdf(size = 1) +
  geom_vline(data = qs, aes(xintercept = med, color = region),
             linetype = "dashed", alpha = .6) +
  scale_x_continuous(
    trans = "log1p",
    breaks = ticks, labels = ticks,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = "Regulon size (# genes)", y = "ECDF",
       title = "ECDF of regulon sizes (dashed = median)") +
  theme_classic()

# 4) (Optional) Non-parametric test: are regulon sizes different across regions?
kruskal.test(size ~ region, data = sizes)
pairwise.wilcox.test(sizes$size, sizes$region, p.adjust.method = "BH")






##---------Lung Science Advance paper---#################################################################################################
##---------Lung Science Advance paper---#################################################################################################
##---------Lung Science Advance paper---#################################################################################################
##---------Lung Science Advance paper---#################################################################################################

#-------Part I---- Convert files into the long format we need
# ---- packages ----
library(tidyverse)

# ---- Parameters: fill in your four gmt file paths here ----
gmt_files <- c(
  "scRBP_prune_matched_results_min10genes_5UTR_lung.gmt",
  "scRBP_prune_matched_results_min10genes_CDS_lung.gmt",
  "scRBP_prune_matched_results_min10genes_Introns_lung.gmt",
  "scRBP_prune_matched_results_min10genes_3UTR_lung.gmt"
)

# Extract region from filename
extract_region <- function(path) {
  m <- stringr::str_match(basename(path), "(5UTR|CDS|Introns|3UTR)")
  ifelse(is.na(m[,2]), NA_character_, m[,2])
}

# Read one gmt -> long table: RBP, region, gene
read_gmt_long <- function(path) {
  region <- extract_region(path)
  lines  <- readr::read_lines(path, progress = FALSE)
  if (length(lines) == 0) return(tibble(RBP=character(), region=character(), gene=character()))
  
  parse_line <- function(line) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(f) < 3) return(NULL)
    name1 <- f[1]
    desc  <- f[2]
    # Get RBP name: prefer desc with "_regulon" stripped; fall back to name1 if missing
    rbp   <- if (!is.na(desc)) gsub("_regulon$", "", desc) else name1
    rbp   <- if (is.na(rbp) || rbp == "") name1 else rbp
    genes <- f[-c(1,2)] |> trimws()
    genes <- genes[genes != "" & !is.na(genes)]
    if (length(genes) == 0) return(NULL)
    tibble(RBP = rbp, region = region, gene = genes)
  }
  purrr::map_dfr(lines, parse_line)
}

# Merge four files
df_long <- purrr::map_dfr(gmt_files, read_gmt_long) |>
  mutate(region = factor(region, levels = c("5UTR","CDS","Introns","3UTR"))) |>
  distinct()  # deduplicate (when the same RBP-region-gene appears multiple times)

# Save
readr::write_tsv(df_long, "rbp_targets_long_lung.tsv")

# Quick check
message("RBPs: ", n_distinct(df_long$RBP),
        " | rows: ", nrow(df_long),
        " | per-region counts:\n",
        paste(capture.output(df_long |> count(region)), collapse = "\n"))



#--------Part II---- heatmap plotting----

# ---- pkgs ----
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ComplexHeatmap)
library(circlize)

# ---- Input: long table ----
df <- read_tsv("rbp_targets_long_lung.tsv", show_col_types = FALSE)
regions <- c("5UTR","CDS","Introns","3UTR")
df <- df %>% filter(region %in% regions) %>%
  mutate(region = factor(region, levels = regions))

# ---- Gene set for each (RBP, region) ----
sets <- df %>% group_by(RBP, region) %>%
  summarise(genes = list(unique(gene)), .groups = "drop")

rbps <- sets %>% pull(RBP) %>% unique() %>% sort()

# =======================
# Heatmap 1: counts (359 x 4)
# =======================
count_df <- sets %>%
  mutate(n = lengths(genes)) %>%
  dplyr::select(RBP, region, n) %>%
  tidyr::complete(RBP = rbps, region = factor(regions, levels = regions), fill = list(n = 0)) %>%
  dplyr::arrange(factor(RBP, levels = rbps))

# Convert to matrix (explicitly set rownames to avoid issues with as.matrix)
count_mat <- count_df %>%
  tidyr::pivot_wider(names_from = region, values_from = n, values_fill = 0) %>%
  as.data.frame()
rownames(count_mat) <- count_mat$RBP
count_mat <- as.matrix(count_mat[, regions, drop = FALSE])


library(pheatmap)

## ---------- Heatmap 1: counts (359 x 4) ----------
# Assumes count_mat already exists (rows = RBP, cols = c("5UTR","CDS","Introns","3UTR"))
pal_cnt <- colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(101)



# Ensure numeric matrix
count_num <- as.matrix(count_mat)
storage.mode(count_num) <- "double"

# --- Discrete breaks (use 0.5 offset to align integer intervals) ---
maxv <- max(count_num, na.rm = TRUE)

# 0, 5-20, 21-50, 51-100, 101-200, >=201
bk <- c(-0.5, 0.5, 20.5, 50.5, 100.5, 200.5, maxv + 0.5)
#bk <- c(-0.5, 0.5, 30.5, 50.5, 150.5, 250.5, maxv + 0.5)

# 6 discrete colors (replace with your preferred palette)
cols <- c("white",           # 0
          "#fee5d9",         # 5-20  (light)
          "#fcae91",         # 21-50
          "#fb6a4a",         # 51-100
          "#de2d26",         # 101-200
          "#a50f15")         # >=201 (dark)

# Legend breaks and labels (corresponding to the intervals above)
leg_brks <- c(0, 12, 35, 75, 150, min(250, maxv))
leg_labs <- c("0", "5-20", "21-50", "51-100", "101-200",
              ifelse(maxv >= 201, "\u2265201", paste0("\u2264", maxv)))

pheatmap(
  count_num,
  color = cols,
  breaks = bk,                 # key: discrete bin mapping
  na_col = "#eeeeee",
  border_color = "grey90",
  cluster_rows = TRUE,
  cluster_cols = F,        # column order fixed as 5UTR/CDS/Introns/3UTR
  show_rownames = TRUE,        # turn on to show row names
  fontsize_row = 5,
  cellwidth = 8, cellheight = 1.1, fontsize = 6,
  legend_breaks = leg_brks,
  legend_labels = leg_labs,
  main = "Regulon size (discrete bins)"
)




# ==============================
# Heatmap 2: pairwise Jaccard (359x6)
# ==============================
# Ensure sets is correct (each row = 1 RBP x region, genes is a list column)
# sets <- df %>% group_by(RBP, region) %>% summarise(genes = list(unique(gene)), .groups = "drop")

# --- Safe grouping approach: split into list by RBP ---
sets$RBP   <- as.character(sets$RBP)     # prevent issues caused by factor/encoding
sets$region<- factor(sets$region, levels = c("5UTR","CDS","Introns","3UTR"))
set_map <- split(sets[, c("region","genes")], sets$RBP)
length(set_map)          # should be approximately 359
head(names(set_map), 3)  # check the first few RBP names

# --- Jaccard helper function ---
jaccard <- function(a, b){
  a <- unique(a); b <- unique(b)
  u <- union(a,b); if (length(u)==0) return(NA_real_)
  length(intersect(a,b)) / length(u)
}

regions <- levels(sets$region)
pairs   <- combn(regions, 2, simplify = FALSE)
pair_names <- sapply(pairs, function(x) paste(x, collapse="|"))

# --- Compute Jaccard for 6 pairs per RBP ---
calc_one <- function(tbl){     # tbl is the data.frame(region, genes) for one RBP
  getg <- function(r){
    idx <- which(tbl$region == r)
    if (length(idx)==0) character(0) else tbl$genes[[idx[1]]]
  }
  sapply(pairs, function(p) jaccard(getg(p[1]), getg(p[2])))
}

jac_mat <- t(sapply(set_map, calc_one))
colnames(jac_mat) <- pair_names


library(pheatmap)

# Color scale and range: 0 ~ 0.30
MAX_COL <- 0.5
#MAX_COL <- 1
pal  <- colorRampPalette(c("#f7fbff", "#74add1", "#08306b"))(121)   # 120-level color scale
#pal  <- colorRampPalette(c("#F7F5CE", "#E77A77", "firebrick3"))(121)   # 120-level color scale
brks <- seq(0, MAX_COL, length.out = 122)                           # one more than number of colors

# (optional) clip display values to 0.30 before plotting
# jac_disp <- pmin(jac_num, MAX_COL)
jac_disp <- pmin(jac_mat, MAX_COL)  # clip values to MAX_COL to prevent exceeding color range
hr <- hclust(dist(jac_mat))   # hierarchical clustering of RBPs (rows) by Euclidean distance
hc <- hclust(dist(t(jac_mat)))  # hierarchical clustering of columns (pairwise Jaccard)

pheatmap(
  jac_disp,
  color = pal, breaks = brks, na_col = "#eeeeee",
  border_color = NA,cellwidth =8, cellheight =1.1,fontsize=6,
  cluster_rows = hr, cluster_cols = hc,          # use the clustering trees computed above
  show_rownames = T,
  frontsize_row = 1,
  legend_breaks = seq(0, MAX_COL, by = 0.05),    # legend tick marks (adjustable)
  legend_labels = sprintf("%.2f", seq(0, MAX_COL, by = 0.05)),
  main = "Pairwise Jaccard across regions (0-0.30)"
)



####--------4 domains---> region sizes distributions

library(readr)
library(dplyr)
library(ggplot2)

# 1) Read in and compute regulon size
regions <- c("5UTR","CDS","Introns","3UTR")

df <- read_tsv("rbp_targets_long_lung.tsv", show_col_types = FALSE) |>
  distinct(RBP, region, gene) |>                               # deduplicate to prevent duplicates
  filter(region %in% regions) |>
  mutate(region = factor(region, levels = regions))

sizes <- df |>
  group_by(region, RBP) |>
  summarise(size = n_distinct(gene), .groups = "drop")         # number of genes per RBP x region

# 2) Basic statistics
summary_tab <- sizes |>
  group_by(region) |>
  summarise(
    n_regulons = n(),
    median = median(size),
    mean   = mean(size),
    IQR    = IQR(size),
    min    = min(size),
    max    = max(size),
    .groups = "drop"
  )
print(summary_tab)

# 3) Plotting (choose one)



###_------ Important ----density plot
library(dplyr)
library(ggplot2)

breaks <- log1p(c(5,10,15,20,30,50,100,200,500,1000))
# Compute median per region (original scale) and convert to log1p scale
meds <- sizes %>%
  group_by(region) %>%
  summarise(med_size = median(size, na.rm = TRUE), .groups = "drop") %>%
  mutate(med_log1p = log1p(med_size))

ggplot(sizes, aes(x = log1p(size), color = region, fill = region)) +
  geom_density(alpha = 0.25, adjust = 1.2, linewidth = 1) +
  geom_vline(
    data = meds,
    aes(xintercept = med_log1p, color = region),
    linetype = "dashed", linewidth = 0.8, alpha = 0.9
  ) +
  scale_x_continuous(
    name = "Regulon size (# genes)",
    breaks = breaks,
    labels = ~ scales::comma(round(expm1(.)))
  ) +
  labs(y = "Density") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic()


## 3a. Violin + boxplot (log scale, suitable for long-tailed distributions)
ggplot(sizes, aes(region, size, fill = region)) +
  geom_violin(trim = FALSE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.size = 0.6, alpha = 0.7, color = "grey20") +
  scale_y_continuous(trans = "log1p") +
  labs(x = "Region", y = "Regulon size (# genes)", title = "Regulon size distribution by region") +
  theme_bw() + theme(legend.position = "none")


ggplot(sizes, aes(region, size, fill = region)) +
  geom_violin(trim = FALSE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.size = 0.6, alpha = 0.7, color = "grey20") +
  scale_y_continuous(trans = "log1p") +
  labs(x = "Region", y = "Regulon size (# genes)", title = "Regulon size distribution by region") +
  theme_classic() + theme(legend.position = "none")

ticks <- c(5, 10, 25, 50, 100, 300, 600, 900)
ggplot(sizes, aes(region, size, fill = region)) +
  geom_violin(trim = FALSE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.size = 0.6, alpha = 0.7, color = "grey20") +
  scale_y_continuous(
    trans = "log1p",
    breaks = ticks, labels = ticks,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = "Region", y = "Regulon size (# genes)",
       title = "Regulon size distribution by region") +
  theme_classic() + theme(legend.position = "none")


## 3b. Faceted histogram (one panel per region)
ggplot(sizes, aes(size)) +
  geom_histogram(bins = 30, color = "white") +
  facet_wrap(~ region, scales = "free_y") +
  scale_x_continuous(trans = "log1p") +
  labs(x = "Regulon size (# genes)", y = "Count") +
  theme_bw()

## 3c. ECDF (cumulative distribution curve, useful for comparing overall size differences)
ggplot(sizes, aes(size, color = region)) +
  stat_ecdf(size = 1) +
  scale_x_continuous(trans = "log1p") +
  labs(x = "Regulon size (# genes)", y = "ECDF", title = "ECDF of regulon sizes") +
  theme_classic()

# Compute quantiles per region
sizes %>%
  group_by(region) %>%
  summarise(n = n(),
            q25 = quantile(size, .25),
            median = median(size),
            q75 = quantile(size, .75),
            p90 = quantile(size, .90),
            max = max(size))

# Add median dashed lines to ECDF plot

qs <- sizes %>% group_by(region) %>% summarise(med = median(size))
ggplot(sizes, aes(size, color = region)) +
  stat_ecdf(size = 1) +
  geom_vline(data = qs, aes(xintercept = med, color = region),
             linetype = "dashed", alpha = .7) +
  scale_x_continuous(trans = "log1p") +
  labs(x = "Regulon size (# genes)", y = "ECDF",
       title = "ECDF of regulon sizes (dashed = median)") +
  theme_bw()


qs <- sizes %>% group_by(region) %>% summarise(med = median(size))
ticks <- c(5, 10, 25, 50, 100,200, 300, 500, 1000,1500)
ggplot(sizes, aes(size, color = region)) +
  stat_ecdf(size = 1) +
  geom_vline(data = qs, aes(xintercept = med, color = region),
             linetype = "dashed", alpha = .6) +
  scale_x_continuous(
    trans = "log1p",
    breaks = ticks, labels = ticks,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = "Regulon size (# genes)", y = "ECDF",
       title = "ECDF of regulon sizes (dashed = median)") +
  theme_classic()

# 4) (Optional) Non-parametric test: are regulon sizes different across regions?
kruskal.test(size ~ region, data = sizes)
pairwise.wilcox.test(sizes$size, sizes$region, p.adjust.method = "BH")
