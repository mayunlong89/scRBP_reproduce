
# --- 2026-01-11
# Fetal brain atlas: re-clustering, UMAP embedding, and visualization
# 599K cells, integrated assay, multiple UMAP parameter explorations


# "Type.v2"---> For main Figure 5A.
# "Subtype.v2" ---> For Sup. Figure 19F.
# "GW_group" ---> For Sup. Figure 19A.
# "Dataset" ---> For Sup. Figure 19B.



# ============================================================
# Load required packages
# ============================================================
library(Seurat)
library(ggplot2)
library(future)
library(parallelly)

# ============================================================
# Configure parallelization
# ============================================================
avail   <- parallelly::availableCores()
workers <- max(1, min(avail, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", avail))))
message(sprintf("Using %d worker(s) [availableCores=%d]", workers, avail))

if (workers <= 1) {
  plan(sequential)
} else {
  plan(multisession, workers = workers)
}
options(future.globals.maxSize = 30 * 1024^3)

# ============================================================
# Part 1: Re-cluster and compute new UMAP embeddings
# ============================================================

# Load the integrated Seurat object
seu <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells.rds")
DefaultAssay(seu) <- "integrated"
set.seed(12345)

# --- PCA (100 components for downstream flexibility) ---
seu <- RunPCA(seu, npcs = 100, verbose = FALSE, reduction.name = "pca_v2")
dims_use <- 1:50

# --- Neighbor graph and multi-resolution Leiden clustering ---
seu <- FindNeighbors(
  seu, reduction = "pca_v2", dims = dims_use,
  k.param = 30, nn.method = "annoy", annoy.metric = "euclidean",
  graph.name = "integrated_snn_v2", verbose = FALSE
)

seu <- FindClusters(
  seu, graph.name = "integrated_snn_v2",
  resolution = c(0.2, 0.4, 0.8, 1.2, 2.0),
  algorithm = 4, verbose = FALSE
)

Idents(seu) <- "integrated_snn_v2_res.1.2"

# --- Standard UMAP (cosine distance) ---
seu <- RunUMAP(
  seu, reduction = "pca_v2", dims = dims_use,
  n.neighbors = 40, min.dist = 0.25, metric = "cosine",
  seed.use = 12345, verbose = FALSE,
  reduction.name = "umap_v2"
)

# ============================================================
# Part 2: Continuum UMAP (trajectory-like layout)
# ============================================================
# Higher n.neighbors + spectral init produces a smoother developmental arc
seu <- RunUMAP(
  seu, reduction = "pca_v2", dims = dims_use,
  n.neighbors = 120, min.dist = 0.6, metric = "cosine",
  init = "spectral", seed.use = 9, verbose = FALSE,
  reduction.name = "umap_continuum"
)

# ============================================================
# Part 3: Define gestational week groups (9 developmental stages)
# ============================================================
gw_raw <- trimws(as.character(seu$Gestational_week))
gw_num <- suppressWarnings(as.numeric(gw_raw))

seu$GW_group <- dplyr::case_when(
  !is.na(gw_num) & gw_num >=  6 & gw_num < 10 ~ "1: GW 6-10",
  !is.na(gw_num) & gw_num >= 10 & gw_num < 12 ~ "2: GW 10-12",
  !is.na(gw_num) & gw_num >= 12 & gw_num < 15 ~ "3: GW 12-15",
  !is.na(gw_num) & gw_num >= 15 & gw_num < 18 ~ "4: GW 15-18",
  !is.na(gw_num) & gw_num >= 18 & gw_num < 21 ~ "5: GW 18-21",
  !is.na(gw_num) & gw_num >= 21 & gw_num < 26 ~ "6: GW 21-26",
  !is.na(gw_num) & gw_num >= 26 & gw_num <= 40 ~ "7: GW 26-40",
  grepl("postnatal", gw_raw, ignore.case = TRUE) ~ "8: 8 mo postnatal",
  TRUE ~ NA_character_
)

seu$GW_group <- factor(seu$GW_group, levels = c(
  "1: GW 6-10", "2: GW 10-12", "3: GW 12-15", "4: GW 15-18",
  "5: GW 18-21", "6: GW 21-26", "7: GW 26-40", "8: 8 mo postnatal"
), ordered = TRUE)

# ============================================================
# Part 4: Cell type color palette
# ============================================================
seu$Type.v2 <- factor(seu$Type.v2)
lev <- levels(seu$Type.v2)

my_cols <- c(
  "Excitatory Neuron" = "#8FBF88",
  "Inhibitory Neuron" = "#3A7D44",
  "Newborn Neuron"    = "#A5D6A7",
  "IPC"               = "#66C2A5",
  "RG"                = "#C79ECF",
  "Div"               = "#F4A259",
  "Astrocyte"         = "#F28E2B",
  "Microglia"         = "#76B7B2",
  "OPC"               = "#4E79A7",
  "OPC.div"           = "#A0CBE8",
  "Oligodendrocyte"   = "#59A14F",
  "Endothelial"       = "#E15759",
  "Pericyte"          = "#B07AA1",
  "Red blood cells"   = "#D37295"
)

# Auto-fill missing cell types with grey
missing_lv <- setdiff(lev, names(my_cols))
if (length(missing_lv) > 0) {
  my_cols <- c(my_cols, setNames(rep("#BDBDBD", length(missing_lv)), missing_lv))
}
my_cols <- my_cols[lev]

# ============================================================
# Part 5: Generate all UMAP plots
# ============================================================

# --- 5a. Cell type UMAP (unlabeled) ---
p_ct <- DimPlot(seu, reduction = "umap_continuum", group.by = "Type.v2",
                raster = TRUE, shuffle = TRUE) +
  scale_color_manual(values = my_cols)

ggsave("umap_continuum_celltype.pdf",  p_ct, width = 9, height = 7, device = cairo_pdf)
ggsave("umap_continuum_celltype.tiff", p_ct, width = 9, height = 7,
       units = "in", dpi = 600, compression = "lzw")

# --- 5b. Cell type UMAP (labeled, no legend) ---
p_ct_lab <- DimPlot(seu, reduction = "umap_continuum", group.by = "Type.v2",
                    raster = TRUE, shuffle = TRUE,
                    label = TRUE, repel = TRUE, label.size = 4) +
  scale_color_manual(values = my_cols) + NoLegend()

ggsave("umap_continuum_celltype_labeled.tiff", p_ct_lab, width = 9, height = 7,
       units = "in", dpi = 600, compression = "lzw")

# --- 5c. Gestational week group UMAP (blue gradient) ---
cols_stage_blue <- c(
  "1: GW 6-10"        = "#DCEBFA",
  "2: GW 10-12"       = "#BBD7F5",
  "3: GW 12-15"       = "#8CC6EA",
  "4: GW 15-18"       = "#5BB7D8",
  "5: GW 18-21"       = "#2BA7C7",
  "6: GW 21-26"       = "#1D86B8",
  "7: GW 26-40"       = "#2E5FA8",
  "8: 8 mo postnatal" = "#7A3E9D"
)

p_gw_blue <- DimPlot(seu, reduction = "umap_continuum", group.by = "GW_group",
                     raster = TRUE) +
  scale_color_manual(values = cols_stage_blue, na.value = "grey80") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave("umap_continuum_GW_group_blue.tiff", p_gw_blue, width = 9, height = 7,
       units = "in", dpi = 300, compression = "lzw")

# --- 5d. Gestational week group UMAP (green gradient) ---
cols_stage_green <- c(
  "1: GW 6-10"        = "#F4F8F4",
  "2: GW 10-12"       = "#E6F2E6",
  "3: GW 12-15"       = "#D0E7D0",
  "4: GW 15-18"       = "#B8DAB8",
  "5: GW 18-21"       = "#9CCB9C",
  "6: GW 21-26"       = "#79B679",
  "7: GW 26-40"       = "#4C9A4C",
  "8: 8 mo postnatal" = "#2F7D2F"
)

p_gw_green <- DimPlot(seu, reduction = "umap_continuum", group.by = "GW_group",
                      raster = TRUE) +
  scale_color_manual(values = cols_stage_green) +
  theme(legend.title = element_blank())

ggsave("umap_continuum_GW_group_green.tiff", p_gw_green, width = 9, height = 7,
       units = "in", dpi = 300, compression = "lzw")

# --- 5e. Dataset of origin UMAP ---
seu$Dataset <- factor(seu$Dataset, levels = c(
  "Bhaduri2021", "Fan2018", "Fan2020",
  "Nowakowski2017", "Polioudakis2019",
  "Smith2021", "Zhong2018"
))

dataset_cols <- c(
  "Bhaduri2021"      = "#F0C0C0",
  "Fan2018"          = "#66A61E",
  "Fan2020"          = "#1F9E89",
  "Nowakowski2017"   = "#7570B3",
  "Polioudakis2019"  = "#1B7837",
  "Smith2021"        = "#D95F02",
  "Zhong2018"        = "#CB5BF0"
)

p_dataset <- DimPlot(seu, reduction = "umap_continuum", group.by = "Dataset",
                     raster = TRUE) +
  scale_color_manual(values = dataset_cols, drop = FALSE) +
  theme(legend.title = element_blank())

ggsave("umap_continuum_Dataset.tiff", p_dataset, width = 9, height = 7,
       units = "in", dpi = 300, compression = "lzw")

# ============================================================
# Part 6: Alternative auto-generated palettes (optional)
# ============================================================

# Option A: Glasbey palette (high-contrast discrete colors for many classes)
library(pals)
pal_gb <- setNames(pals::glasbey(length(lev)), lev)

p_gb <- DimPlot(seu, reduction = "umap_continuum", group.by = "Type.v2",
                raster = TRUE, shuffle = TRUE) +
  scale_color_manual(values = pal_gb) + NoLegend()
ggsave("umap_continuum_glasbey.pdf", p_gb, width = 9, height = 7, device = cairo_pdf)

# Option B: NPG (Nature Publishing Group) palette
library(ggsci)
pal_npg <- setNames(rep(ggsci::pal_npg("nrc")(10), length.out = length(lev)), lev)

p_npg <- DimPlot(seu, reduction = "umap_continuum", group.by = "Type.v2",
                 raster = TRUE, shuffle = TRUE) +
  scale_color_manual(values = pal_npg) + NoLegend()
ggsave("umap_continuum_npg.pdf", p_npg, width = 9, height = 7, device = cairo_pdf)

# ============================================================
# Save final object
# ============================================================
saveRDS(seu, file = "clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

