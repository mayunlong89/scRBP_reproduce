# --- 2026-02-10
# Fetal brain UMAP visualization (599K cells, continuum embedding)

# ============================================================
# Load required packages
# ============================================================
library(Seurat)
library(ggplot2)

# ============================================================
# Step 1: Load Seurat object
# ============================================================
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

# ============================================================
# Step 2: UMAP using the continuum embedding (final version)
# ============================================================

# Unlabeled UMAP
p7 <- DimPlot(
  seurat_obj,
  reduction = "umap_continuum",
  group.by  = "Type.v2",
  raster    = TRUE,
  shuffle   = TRUE
)

ggsave("umap_v2_Type_v2_umap_continuum_nonlabeled.pdf",  p7, width = 9, height = 7, dpi = 600)
ggsave("umap_v2_Type_v2_umap_continuum_nonlabeled.png",  p7, width = 9, height = 7, dpi = 600)
ggsave("umap_v2_Type_v2_umap_continuum_nonlabeled.tiff", p7, width = 9, height = 7,
       units = "in", dpi = 600, compression = "lzw")

# Labeled UMAP (with repelled text annotations, no legend)
p8 <- DimPlot(
  seurat_obj,
  reduction  = "umap_continuum",
  group.by   = "Type.v2",
  raster     = TRUE,
  shuffle    = TRUE,
  label      = TRUE,
  repel      = TRUE,
  label.size = 4
) + NoLegend()

ggsave("umap_v2_Type_v2_umap_continuum_labeled.tiff", p8, width = 9, height = 7,
       units = "in", dpi = 600, compression = "lzw")
ggsave("umap_v2_Type_v2_umap_continuum_labeled.png",  p8, width = 9, height = 7, dpi = 600)

# ============================================================
# Step 3: Custom color palette (manually curated, journal-style)
# ============================================================
seurat_obj$Type.v2 <- factor(seurat_obj$Type.v2)
lev <- levels(seurat_obj$Type.v2)

my_cols <- c(
  "Excitatory Neuron" = "#3A7D44",
  "Inhibitory Neuron" = "#8FBF88",
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

# Auto-fill any cell types present in data but missing from the palette
missing_lv <- setdiff(lev, names(my_cols))
if (length(missing_lv) > 0) {
  my_cols <- c(my_cols, setNames(rep("#BDBDBD", length(missing_lv)), missing_lv))
}

# Reorder palette to match factor levels (prevents color misalignment)
my_cols <- my_cols[lev]

p_named <- DimPlot(seurat_obj, reduction = "umap_continuum",
                   group.by = "Type.v2", raster = TRUE, shuffle = TRUE) +
  scale_color_manual(values = my_cols) +
  NoLegend()

ggsave("umap_continuum_named.pdf", p_named, width = 9, height = 7, device = cairo_pdf)
ggsave("umap_continuum_named.png", p_named, width = 9, height = 7, dpi = 300)

# ============================================================
# Step 4: Alternative palettes for many cell types (>20 classes)
# ============================================================

# Option A: Glasbey palette (high-contrast discrete colors)
library(pals)
pal_gb <- pals::glasbey(length(lev))
names(pal_gb) <- lev

p_gb <- DimPlot(seurat_obj, reduction = "umap_continuum",
                group.by = "Type.v2", raster = TRUE, shuffle = TRUE) +
  scale_color_manual(values = pal_gb) + NoLegend()

ggsave("umap_continuum_glasbey.pdf", p_gb, width = 9, height = 7, device = cairo_pdf)

# Option B: NPG (Nature Publishing Group) palette via ggsci
library(ggsci)
pal_npg <- rep(ggsci::pal_npg("nrc")(10), length.out = length(lev))
names(pal_npg) <- lev

p_npg <- DimPlot(seurat_obj, reduction = "umap_continuum",
                 group.by = "Type.v2", raster = TRUE, shuffle = TRUE) +
  scale_color_manual(values = pal_npg) + NoLegend()

ggsave("umap_continuum_npg.pdf", p_npg, width = 9, height = 7, device = cairo_pdf)
