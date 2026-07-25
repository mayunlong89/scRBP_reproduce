#-----2025-09-04
library(Seurat)
library(ggplot2)
library(grid)  # for unit()
obj <- readRDS("10X_PBMC_subset_6.5K.rds")         # your Seurat object
                 # RBP column name to plot
# 1) Read TRS matrix (rows = cells, columns = RBPs); keep cells present in the object and align to Seurat order
trs <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_6.5Kcell_100runs_results/01_monocyte_count_trs/10X_PBMC_subset_6.5K_100runs_TRS_sc_monocyte_count_5UTR.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)
trs <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_6.5Kcell_100runs_results/01_monocyte_count_trs/10X_PBMC_subset_6.5K_100runs_TRS_sc_monocyte_count_Introns.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)
trs <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_6.5Kcell_100runs_results/01_monocyte_count_trs/10X_PBMC_subset_6.5K_100runs_TRS_sc_monocyte_count_3UTR.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)
common <- intersect(Cells(obj), rownames(trs))
obj    <- subset(obj, cells = common)
rbp <- "QKI"
rbp <- "YTHDF3"
rbp <- "RBM47"  
rbp <- "CELF1"
rbp <- "YBX3"
rbp <- "RYBP"
rbp <- "RPRD2"
rbp <- "TAF15"
rbp <- "HNRNPC"
v      <- trs[common, rbp, drop = TRUE]
# 2) (Optional) Quantile clipping to prevent extreme values from distorting the color scale
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
# 3) Store in meta.data and plot (using umap.rna as coordinates)
obj[[paste0(rbp, "_TRS")]] <- v
# Finer grey -> white -> orange-red gradient (you can also change to your preferred two/three endpoints)
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)
p <- FeaturePlot(
  obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap.rna",
  pt.size   = 1,
  order     = TRUE
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,        # smoother gradient in the legend
      raster = TRUE,     # use raster rendering to avoid banding artifacts
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP (umap.rna)"))
ggsave("./subset6.5K_figures_plots/HNRNPC_monocyte_count.umap_trs.pdf", p, width = 7, height = 5)
ggsave("./subset6.5K_figures_plots/TAF15_monocyte_count.umap_trs.pdf", p, width = 7, height = 5)
ggsave("QKI_monocyte_count.umap_trs.pdf", p, width = 7, height = 5)
ggsave("YTHDF3_monocyte_count.umap_trs.pdf", p, width = 7, height = 5)
ggsave("RBM47_monocyte_count.umap_trs.pdf", p, width = 7, height = 5)
ggsave("CELF1_monocyte_count.umap_trs.pdf", p, width = 7, height = 5)
ggsave("./subset6.5K_figures_plots/YBX3_monocyte_count.umap_trs.pdf", p, width = 7, height = 5)
ggsave("./subset6.5K_figures_plots/RYBP_monocyte_count.umap_trs.pdf", p, width = 7, height = 5)
ggsave("./subset6.5K_figures_plots/RPRD2_monocyte_count.umap_trs.pdf", p, width = 7, height = 5)
