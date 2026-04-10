

# ######################



# FeaturePlot---Sig RBP regulon -UMAP---590K---SCZ



# ######################

# #-----SCZ-3UTR---FeaturePlot

```R
#-----2026-2-10


library(Seurat)
library(ggplot2)
library(grid)  # 用于 unit()
library(SeuratDisk)

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")



#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ")

#----------SCZ 3UTR
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_SCZ_3UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ/z_SCZ_TRS/SCZ_processed_magma_results_590Kcells_3UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


common <- intersect(Cells(seurat_obj), rownames(trs_SCZ_3UTR))
seurat_obj  <- subset(seurat_obj, cells = common)

#SCZ 3UTR
rbp <- "ELAVL4"
rbp <- "CELF2"
rbp <- "CPEB4"
rbp <- "EIF3B"

rbp <- "RBFOX2"


rbp <- "HNRNPR"
rbp <- "HNRNPU"

rbp <- "SRRM4"

rbp <- "SRSF6"

rbp <- "YBX3"






v <- trs_SCZ_3UTR[common, rbp, drop = TRUE]

# 2) （可选）分位数截断，避免极端值把色标拉爆
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])

# 3) 写进 meta.data 并画图（用 umap.rna 做坐标）
seurat_obj[[paste0(rbp, "_TRS")]] <- v


# 更细腻的灰→白→橙红渐变（你也可以改成你喜欢的两个/三个端点）
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)

p <- FeaturePlot(
  seurat_obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap_continuum",
  pt.size   = 0.3,        # 你可以直接 4~6 试
  order     = TRUE,
  raster    = FALSE     # 关键：关闭栅格化
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,
      raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP"))



 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
ggsave(
  "./ELAVL4_3UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  "./CELF2_3UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./CPEB4_3UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  "./EIF3B_3UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./HNRNPR_3UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./HNRNPU_3UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./RBFOX2_3UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  "./SRRM4_3UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)




ggsave(
  "./SRSF6_3UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./YBX3_3UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)

```



## SCZ-3UTR boxplot

```bash
#---2026-2-10

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ")


#----------07_SCZ 3UTR
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_SCZ_3UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ/z_SCZ_TRS/SCZ_processed_magma_results_590Kcells_3UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


#SCZ 3UTR
rbp <- "ELAVL4"
rbp <- "CELF2"
rbp <- "CPEB4"
rbp <- "EIF3B"

rbp <- "RBFOX2"


rbp <- "HNRNPR"
rbp <- "HNRNPU"

rbp <- "SRRM4"

rbp <- "SRSF6"

rbp <- "YBX3"




###-Step 1 function: plot_trs_box_by_celltype()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

plot_trs_box_by_celltype <- function(seurat_obj,
                                     trs_vec,
                                     celltype_col,
                                     rbp_label = "RBP",
                                     pal_cols = c("#FCEFEA", "#F7C4B2", "#F26B4F"),
                                     ylab = "scRBP TRS",
                                     sort_by = c("median","mean"),
                                     show_outliers = FALSE,
                                     box_width = 0.65) {

  sort_by <- match.arg(sort_by)

  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  df <- data.frame(
    cell      = names(trs_vec),
    TRS       = as.numeric(trs_vec),
    cell_type = seurat_obj@meta.data[names(trs_vec), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      n = dplyr::n(),
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # 颜色渐进：按排序后的 cell type 映射到渐变调色板
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp_label, " TRS across cell types"),
      x = NULL, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}


###-----Step 1: Processing data
# v: 保证顺序=Seurat cells，并且带names
v <- trs_SCZ_3UTR[Cells(seurat_obj), rbp, drop = TRUE]
names(v) <- Cells(seurat_obj)

# 截断（不丢names的写法）
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
names(v) <- Cells(seurat_obj)   # 再保险补一次


###-----Step 2: run function
celltype_col <- "Type.v2"   # 你实际用的列名
 
col_m     = c("#FCEFEA", "#F7C4B2", "#F26B4F")
p_box <- plot_trs_box_by_celltype(
  seurat_obj   = seurat_obj,
  trs_vec      = v,
  celltype_col = celltype_col,
  rbp_label    = rbp,
  sort_by = c("median"),
  pal_cols     = rev(col_m)  # 图2那种浅->深
)


###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
ggsave(
  filename = "./ELAVL4_3UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)



ggsave(
  filename = "./CELF2_3UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename = "./CPEB4_3UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename = "./EIF3B_3UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename = "./RBFOX2_3UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)




ggsave(
  filename = "./HNRNPR_3UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "./HNRNPU_3UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename = "./SRRM4_3UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename = "./SRSF6_3UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename = "./YBX3_3UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)




```













# #-----SCZ-5UTR---FeaturePlot

```R
#-----2026-2-10


library(Seurat)
library(ggplot2)
library(grid)  # 用于 unit()
library(SeuratDisk)

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")



#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ")

#----------SCZ 5UTR
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_SCZ_5UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ/z_SCZ_TRS/SCZ_processed_magma_results_590Kcells_5UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


common <- intersect(Cells(seurat_obj), rownames(trs_SCZ_5UTR))
seurat_obj  <- subset(seurat_obj, cells = common)

#SCZ 5UTR
rbp <- "CELF4"
rbp <- "EIF4A1"
rbp <- "EIF4G2"

rbp <- "RBFOX2"


rbp <- "RBM17"


rbp <- "HNRNPA1"




rbp <- "HNRNPH1"



 




v <- trs_SCZ_5UTR[common, rbp, drop = TRUE]

# 2) （可选）分位数截断，避免极端值把色标拉爆
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])

# 3) 写进 meta.data 并画图（用 umap.rna 做坐标）
seurat_obj[[paste0(rbp, "_TRS")]] <- v


# 更细腻的灰→白→橙红渐变（你也可以改成你喜欢的两个/三个端点）
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)

p <- FeaturePlot(
  seurat_obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap_continuum",
  pt.size   = 0.3,        # 你可以直接 4~6 试
  order     = TRUE,
  raster    = FALSE     # 关键：关闭栅格化
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,
      raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP"))



 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
ggsave(
  "./CELF4_5UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  "./EIF4A1_5UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  "./EIF4G2_5UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)


ggsave(
  "./RBFOX2_5UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  "./RBM17_5UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)



ggsave(
  "./HNRNPA1_5UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)



ggsave(
  "./HNRNPH1_5UTR_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)



```



## SCZ-5UTR boxplot

```bash
#---2026-2-10

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ")


#----------07_SCZ 5UTR
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_SCZ_5UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ/z_SCZ_TRS/SCZ_processed_magma_results_590Kcells_5UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


#SCZ 5UTR

rbp <- "CELF4"


rbp <- "EIF4A1"

rbp <- "EIF4G2"

rbp <- "RBFOX2"


rbp <- "RBM17"


rbp <- "HNRNPA1"




rbp <- "HNRNPH1"





###-Step 1 function: plot_trs_box_by_celltype()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

plot_trs_box_by_celltype <- function(seurat_obj,
                                     trs_vec,
                                     celltype_col,
                                     rbp_label = "RBP",
                                     pal_cols = c("#FCEFEA", "#F7C4B2", "#F26B4F"),
                                     ylab = "scRBP TRS",
                                     sort_by = c("median","mean"),
                                     show_outliers = FALSE,
                                     box_width = 0.65) {

  sort_by <- match.arg(sort_by)

  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  df <- data.frame(
    cell      = names(trs_vec),
    TRS       = as.numeric(trs_vec),
    cell_type = seurat_obj@meta.data[names(trs_vec), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      n = dplyr::n(),
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # 颜色渐进：按排序后的 cell type 映射到渐变调色板
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp_label, " TRS across cell types"),
      x = NULL, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}


###-----Step 1: Processing data
# v: 保证顺序=Seurat cells，并且带names
v <- trs_SCZ_5UTR[Cells(seurat_obj), rbp, drop = TRUE]
names(v) <- Cells(seurat_obj)

# 截断（不丢names的写法）
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
names(v) <- Cells(seurat_obj)   # 再保险补一次


###-----Step 2: run function
celltype_col <- "Type.v2"   # 你实际用的列名
 
col_m     = c("#FCEFEA", "#F7C4B2", "#F26B4F")
p_box <- plot_trs_box_by_celltype(
  seurat_obj   = seurat_obj,
  trs_vec      = v,
  celltype_col = celltype_col,
  rbp_label    = rbp,
  sort_by = c("median"),
  pal_cols     = rev(col_m)  # 图2那种浅->深
)


###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
ggsave(
  filename = "./CELF4_5UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)



ggsave(
  filename = "./EIF4A1_5UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename = "./EIF4G2_5UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)



ggsave(
  filename = "./RBFOX2_5UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)



ggsave(
  filename = "./RBM17_5UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename = "./HNRNPA1_5UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename = "./HNRNPH1_5UTR_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)




```















# #-----SCZ-CDS---FeaturePlot

```R
#-----2026-2-10


library(Seurat)
library(ggplot2)
library(grid)  # 用于 unit()
library(SeuratDisk)

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")



#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ")

#----------SCZ CDS
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_SCZ_CDS <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ/z_SCZ_TRS/SCZ_processed_magma_results_590Kcells_CDS_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


common <- intersect(Cells(seurat_obj), rownames(trs_SCZ_CDS))
seurat_obj  <- subset(seurat_obj, cells = common)

#SCZ CDS
rbp <- "HNRNPH1"

rbp <- "NOVA1"

rbp <- "PABPC1"


rbp <- "RBFOX1"
 
rbp <- "RPS3"
 
rbp <- "SRRM4"
 
rbp <- "SRSF9"
 
 rbp <- "YWHAE"
 
 rbp <- "YWHAH"
 


v <- trs_SCZ_CDS[common, rbp, drop = TRUE]

# 2) （可选）分位数截断，避免极端值把色标拉爆
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])

# 3) 写进 meta.data 并画图（用 umap.rna 做坐标）
seurat_obj[[paste0(rbp, "_TRS")]] <- v


# 更细腻的灰→白→橙红渐变（你也可以改成你喜欢的两个/三个端点）
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)

p <- FeaturePlot(
  seurat_obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap_continuum",
  pt.size   = 0.3,        # 你可以直接 4~6 试
  order     = TRUE,
  raster    = FALSE     # 关键：关闭栅格化
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,
      raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP"))






 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
ggsave(
  "./HNRNPH1_CDS_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)





ggsave(
  "./NOVA1_CDS_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  "./PABPC1_CDS_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)


ggsave(
  "./RBFOX1_CDS_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./RPS3_CDS_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./SRRM4_CDS_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)



ggsave(
  "./SRSF9_CDS_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)




ggsave(
  "./YWHAE_CDS_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)




ggsave(
  "./YWHAH_CDS_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)









```



## SCZ-CDS boxplot

```bash
#---2026-2-10

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ")


#----------07_SCZ CDS
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_SCZ_CDS <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ/z_SCZ_TRS/SCZ_processed_magma_results_590Kcells_CDS_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


#SCZ CDS

rbp <- "HNRNPH1"

rbp <- "NOVA1"

rbp <- "PABPC1"


rbp <- "RBFOX1"
 
rbp <- "RPS3"
 
rbp <- "SRRM4"
 
rbp <- "SRSF9"
 
 rbp <- "YWHAE"
 
 rbp <- "YWHAH"
 





###-Step 1 function: plot_trs_box_by_celltype()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

plot_trs_box_by_celltype <- function(seurat_obj,
                                     trs_vec,
                                     celltype_col,
                                     rbp_label = "RBP",
                                     pal_cols = c("#FCEFEA", "#F7C4B2", "#F26B4F"),
                                     ylab = "scRBP TRS",
                                     sort_by = c("median","mean"),
                                     show_outliers = FALSE,
                                     box_width = 0.65) {

  sort_by <- match.arg(sort_by)

  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  df <- data.frame(
    cell      = names(trs_vec),
    TRS       = as.numeric(trs_vec),
    cell_type = seurat_obj@meta.data[names(trs_vec), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      n = dplyr::n(),
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # 颜色渐进：按排序后的 cell type 映射到渐变调色板
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp_label, " TRS across cell types"),
      x = NULL, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}



###-----Step 1: Processing data
# v: 保证顺序=Seurat cells，并且带names
v <- trs_SCZ_CDS[Cells(seurat_obj), rbp, drop = TRUE]
names(v) <- Cells(seurat_obj)

# 截断（不丢names的写法）
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
names(v) <- Cells(seurat_obj)   # 再保险补一次


###-----Step 2: run function
celltype_col <- "Type.v2"   # 你实际用的列名
 
col_m     = c("#FCEFEA", "#F7C4B2", "#F26B4F")
p_box <- plot_trs_box_by_celltype(
  seurat_obj   = seurat_obj,
  trs_vec      = v,
  celltype_col = celltype_col,
  rbp_label    = rbp,
  sort_by = c("median"),
  pal_cols     = rev(col_m)  # 图2那种浅->深
)


###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
ggsave(
  filename = "./HNRNPH1_CDS_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)



ggsave(
  filename = "./NOVA1_CDS_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)



ggsave(
  filename = "./PABPC1_CDS_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)




ggsave(
  filename = "./RBFOX1_CDS_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)



ggsave(
  filename = "./RPS3_CDS_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename = "./SRRM4_CDS_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)




ggsave(
  filename = "./SRSF9_CDS_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename = "./YWHAE_CDS_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename = "./YWHAH_CDS_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

```















# #-----SCZ-Introns---FeaturePlot

```R
#-----2026-2-10


library(Seurat)
library(ggplot2)
library(grid)  # 用于 unit()
library(SeuratDisk)

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")



#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ")

#----------SCZ Introns
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_SCZ_Introns <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ/z_SCZ_TRS/SCZ_processed_magma_results_590Kcells_Introns_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


common <- intersect(Cells(seurat_obj), rownames(trs_SCZ_Introns))
seurat_obj  <- subset(seurat_obj, cells = common)

#SCZ Introns
rbp <- "CELF2"

rbp <- "DDX5"


rbp <- "KHSRP"


rbp <- "RBFOX2"
 
rbp <- "RBFOX3"
 
rbp <- "TRA2A"
 


v <- trs_SCZ_Introns[common, rbp, drop = TRUE]

# 2) （可选）分位数截断，避免极端值把色标拉爆
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])

# 3) 写进 meta.data 并画图（用 umap.rna 做坐标）
seurat_obj[[paste0(rbp, "_TRS")]] <- v


# 更细腻的灰→白→橙红渐变（你也可以改成你喜欢的两个/三个端点）
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)

p <- FeaturePlot(
  seurat_obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap_continuum",
  pt.size   = 0.3,        # 你可以直接 4~6 试
  order     = TRUE,
  raster    = FALSE     # 关键：关闭栅格化
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,
      raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP"))






 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
ggsave(
  "./CELF2_Introns_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)




ggsave(
  "./DDX5_Introns_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)



ggsave(
  "./KHSRP_Introns_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  "./RBFOX2_Introns_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  "./RBFOX3_Introns_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)



ggsave(
  "./TRA2A_Introns_SCZ.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)

```



## SCZ-Introns boxplot

```bash
#---2026-2-10

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ")


#----------07_SCZ Introns
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_SCZ_CDS <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/07_SCZ/z_SCZ_TRS/SCZ_processed_magma_results_590Kcells_Introns_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


#SCZ Introns

rbp <- "CELF2"

rbp <- "DDX5"

rbp <- "KHSRP"


rbp <- "RBFOX2"
 
rbp <- "RBFOX3"
 
rbp <- "TRA2A"




###-Step 1 function: plot_trs_box_by_celltype()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

plot_trs_box_by_celltype <- function(seurat_obj,
                                     trs_vec,
                                     celltype_col,
                                     rbp_label = "RBP",
                                     pal_cols = c("#FCEFEA", "#F7C4B2", "#F26B4F"),
                                     ylab = "scRBP TRS",
                                     sort_by = c("median","mean"),
                                     show_outliers = FALSE,
                                     box_width = 0.65) {

  sort_by <- match.arg(sort_by)

  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  df <- data.frame(
    cell      = names(trs_vec),
    TRS       = as.numeric(trs_vec),
    cell_type = seurat_obj@meta.data[names(trs_vec), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      n = dplyr::n(),
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # 颜色渐进：按排序后的 cell type 映射到渐变调色板
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp_label, " TRS across cell types"),
      x = NULL, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}



###-----Step 1: Processing data
# v: 保证顺序=Seurat cells，并且带names
v <- trs_SCZ_Introns[Cells(seurat_obj), rbp, drop = TRUE]
names(v) <- Cells(seurat_obj)

# 截断（不丢names的写法）
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
names(v) <- Cells(seurat_obj)   # 再保险补一次


###-----Step 2: run function
celltype_col <- "Type.v2"   # 你实际用的列名
 
col_m     = c("#FCEFEA", "#F7C4B2", "#F26B4F")
p_box <- plot_trs_box_by_celltype(
  seurat_obj   = seurat_obj,
  trs_vec      = v,
  celltype_col = celltype_col,
  rbp_label    = rbp,
  sort_by = c("median"),
  pal_cols     = rev(col_m)  # 图2那种浅->深
)


###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
ggsave(
  filename = "./CELF2_Introns_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename = "./DDX5_Introns_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename = "./KHSRP_Introns_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename = "./RBFOX2_Introns_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)



ggsave(
  filename = "./RBFOX3_Introns_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)



ggsave(
  filename = "./TRA2A_Introns_SCZ.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)








```









# ######################



# FeaturePlot---展示显著的RBP-UMAP---590K---BIP



# ######################

# #-----BIP-3UTR---FeaturePlot

```R
#-----2026-2-10


library(Seurat)
library(ggplot2)
library(grid)  # 用于 unit()
library(SeuratDisk)

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")



#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/04_BIP")

#----------BIP 3UTR
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_BIP_3UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/04_BIP/z_BIP_TRS/BIP_processed_magma_results_590Kcells_3UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


common <- intersect(Cells(seurat_obj), rownames(trs_BIP_3UTR))
seurat_obj  <- subset(seurat_obj, cells = common)

#BIP 3UTR

rbp <- "CELF2"
 


rbp <- "SRRM4"

rbp <- "KHSRP"

 






v <- trs_BIP_3UTR[common, rbp, drop = TRUE]

# 2) （可选）分位数截断，避免极端值把色标拉爆
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])

# 3) 写进 meta.data 并画图（用 umap.rna 做坐标）
seurat_obj[[paste0(rbp, "_TRS")]] <- v


# 更细腻的灰→白→橙红渐变（你也可以改成你喜欢的两个/三个端点）
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)

p <- FeaturePlot(
  seurat_obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap_continuum",
  pt.size   = 0.3,        # 你可以直接 4~6 试
  order     = TRUE,
  raster    = FALSE     # 关键：关闭栅格化
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,
      raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP"))



 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）

ggsave(
  "./CELF2_3UTR_BIP.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./SRRM4_3UTR_BIP.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)




ggsave(
  "./KHSRP_3UTR_BIP.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)





```



## BIP-3UTR boxplot

```bash
#---2026-2-10

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/04_BIP")


#----------BIP 3UTR
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_BIP_3UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/04_BIP/z_BIP_TRS/BIP_processed_magma_results_590Kcells_3UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


#BIP 3UTR
rbp <- "CELF2"
 


rbp <- "SRRM4"

rbp <- "KHSRP"

 

###-Step 1 function: plot_trs_box_by_celltype()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

plot_trs_box_by_celltype <- function(seurat_obj,
                                     trs_vec,
                                     celltype_col,
                                     rbp_label = "RBP",
                                     pal_cols = c("#FCEFEA", "#F7C4B2", "#F26B4F"),
                                     ylab = "scRBP TRS",
                                     sort_by = c("median","mean"),
                                     show_outliers = FALSE,
                                     box_width = 0.65) {

  sort_by <- match.arg(sort_by)

  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  df <- data.frame(
    cell      = names(trs_vec),
    TRS       = as.numeric(trs_vec),
    cell_type = seurat_obj@meta.data[names(trs_vec), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      n = dplyr::n(),
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # 颜色渐进：按排序后的 cell type 映射到渐变调色板
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp_label, " TRS across cell types"),
      x = NULL, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}


###-----Step 1: Processing data
# v: 保证顺序=Seurat cells，并且带names
v <- trs_BIP_3UTR[Cells(seurat_obj), rbp, drop = TRUE]
names(v) <- Cells(seurat_obj)

# 截断（不丢names的写法）
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
names(v) <- Cells(seurat_obj)   # 再保险补一次


###-----Step 2: run function
celltype_col <- "Type.v2"   # 你实际用的列名
 
col_m     = c("#FCEFEA", "#F7C4B2", "#F26B4F")
p_box <- plot_trs_box_by_celltype(
  seurat_obj   = seurat_obj,
  trs_vec      = v,
  celltype_col = celltype_col,
  rbp_label    = rbp,
  sort_by = c("median"),
  pal_cols     = rev(col_m)  # 图2那种浅->深
)


###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format


ggsave(
  filename = "./CELF2_3UTR_BIP.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "./SRRM4_3UTR_BIP.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

 

ggsave(
  filename = "./KHSRP_3UTR_BIP.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)




```















# ######################



# FeaturePlot---展示显著的RBP-UMAP---590K---ASD



# ######################

# #-----ASD-3UTR---FeaturePlot

```R
#-----2026-2-10


library(Seurat)
library(ggplot2)
library(grid)  # 用于 unit()
library(SeuratDisk)

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")



#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/01_ASD")

#----------ASD 3UTR
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_ASD_3UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/01_ASD/z_ASD_TRS/ASD_processed_magma_results_590Kcells_3UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


common <- intersect(Cells(seurat_obj), rownames(trs_ASD_3UTR))
seurat_obj  <- subset(seurat_obj, cells = common)

#ASD 3UTR

rbp <- "CELF2"
 


rbp <- "SRRM4"








v <- trs_ASD_3UTR[common, rbp, drop = TRUE]

# 2) （可选）分位数截断，避免极端值把色标拉爆
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])

# 3) 写进 meta.data 并画图（用 umap.rna 做坐标）
seurat_obj[[paste0(rbp, "_TRS")]] <- v


# 更细腻的灰→白→橙红渐变（你也可以改成你喜欢的两个/三个端点）
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)

p <- FeaturePlot(
  seurat_obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap_continuum",
  pt.size   = 0.3,        # 你可以直接 4~6 试
  order     = TRUE,
  raster    = FALSE     # 关键：关闭栅格化
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,
      raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP"))



 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）

ggsave(
  "./CELF2_3UTR_ASD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./SRRM4_3UTR_ASD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)








```



## ASD-3UTR boxplot

```bash
#---2026-2-10

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/01_ASD")


#---------ASD 3UTR
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_ASD_3UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/01_ASD/z_ASD_TRS/ASD_processed_magma_results_590Kcells_3UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


#ASD 3UTR
rbp <- "CELF2"
 


rbp <- "SRRM4"

 
 

###-Step 1 function: plot_trs_box_by_celltype()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

plot_trs_box_by_celltype <- function(seurat_obj,
                                     trs_vec,
                                     celltype_col,
                                     rbp_label = "RBP",
                                     pal_cols = c("#FCEFEA", "#F7C4B2", "#F26B4F"),
                                     ylab = "scRBP TRS",
                                     sort_by = c("median","mean"),
                                     show_outliers = FALSE,
                                     box_width = 0.65) {

  sort_by <- match.arg(sort_by)

  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  df <- data.frame(
    cell      = names(trs_vec),
    TRS       = as.numeric(trs_vec),
    cell_type = seurat_obj@meta.data[names(trs_vec), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      n = dplyr::n(),
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # 颜色渐进：按排序后的 cell type 映射到渐变调色板
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp_label, " TRS across cell types"),
      x = NULL, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}


###-----Step 1: Processing data
# v: 保证顺序=Seurat cells，并且带names
v <- trs_ASD_3UTR[Cells(seurat_obj), rbp, drop = TRUE]
names(v) <- Cells(seurat_obj)

# 截断（不丢names的写法）
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
names(v) <- Cells(seurat_obj)   # 再保险补一次


###-----Step 2: run function
celltype_col <- "Type.v2"   # 你实际用的列名
 
col_m     = c("#FCEFEA", "#F7C4B2", "#F26B4F")
p_box <- plot_trs_box_by_celltype(
  seurat_obj   = seurat_obj,
  trs_vec      = v,
  celltype_col = celltype_col,
  rbp_label    = rbp,
  sort_by = c("median"),
  pal_cols     = rev(col_m)  # 图2那种浅->深
)


###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format


ggsave(
  filename = "./CELF2_3UTR_ASD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "./SRRM4_3UTR_ASD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

 
 


```

















# ######################

# #-----ASD-5UTR---FeaturePlot

```R
#-----2026-2-10


library(Seurat)
library(ggplot2)
library(grid)  # 用于 unit()
library(SeuratDisk)

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")



#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/01_ASD")

#----------ASD 5UTR
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_ASD_5UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/01_ASD/z_ASD_TRS/ASD_processed_magma_results_590Kcells_5UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


common <- intersect(Cells(seurat_obj), rownames(trs_ASD_5UTR))
seurat_obj  <- subset(seurat_obj, cells = common)

#ASD 5UTR

rbp <- "CELF4"
 

rbp <- "DDX24"

rbp <- "FTO"

rbp <- "ILF3" 

rbp <- "RBMX"







v <- trs_ASD_5UTR[common, rbp, drop = TRUE]

# 2) （可选）分位数截断，避免极端值把色标拉爆
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])

# 3) 写进 meta.data 并画图（用 umap.rna 做坐标）
seurat_obj[[paste0(rbp, "_TRS")]] <- v


# 更细腻的灰→白→橙红渐变（你也可以改成你喜欢的两个/三个端点）
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)

p <- FeaturePlot(
  seurat_obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap_continuum",
  pt.size   = 0.3,        # 你可以直接 4~6 试
  order     = TRUE,
  raster    = FALSE     # 关键：关闭栅格化
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,
      raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP"))



 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）

ggsave(
  "./CELF4_5UTR_ASD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./DDX24_5UTR_ASD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)


ggsave(
  "./FTO_5UTR_ASD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  "./ILF3_5UTR_ASD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  "./RBMX_5UTR_ASD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)





```



## ASD-5UTR boxplot

```bash
#---2026-2-10

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/01_ASD")


#---------ASD 5UTR
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_ASD_5UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/01_ASD/z_ASD_TRS/ASD_processed_magma_results_590Kcells_5UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


#ASD 5UTR
rbp <- "CELF4"
 


rbp <- "DDX24"

rbp <- "FTO"

 
rbp <- "ILF3"
rbp <- "RBMX"

 
 
 
 

###-Step 1 function: plot_trs_box_by_celltype()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

plot_trs_box_by_celltype <- function(seurat_obj,
                                     trs_vec,
                                     celltype_col,
                                     rbp_label = "RBP",
                                     pal_cols = c("#FCEFEA", "#F7C4B2", "#F26B4F"),
                                     ylab = "scRBP TRS",
                                     sort_by = c("median","mean"),
                                     show_outliers = FALSE,
                                     box_width = 0.65) {

  sort_by <- match.arg(sort_by)

  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  df <- data.frame(
    cell      = names(trs_vec),
    TRS       = as.numeric(trs_vec),
    cell_type = seurat_obj@meta.data[names(trs_vec), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      n = dplyr::n(),
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # 颜色渐进：按排序后的 cell type 映射到渐变调色板
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp_label, " TRS across cell types"),
      x = NULL, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}


###-----Step 1: Processing data
# v: 保证顺序=Seurat cells，并且带names
v <- trs_ASD_5UTR[Cells(seurat_obj), rbp, drop = TRUE]
names(v) <- Cells(seurat_obj)

# 截断（不丢names的写法）
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
names(v) <- Cells(seurat_obj)   # 再保险补一次


###-----Step 2: run function
celltype_col <- "Type.v2"   # 你实际用的列名
 
col_m     = c("#FCEFEA", "#F7C4B2", "#F26B4F")
p_box <- plot_trs_box_by_celltype(
  seurat_obj   = seurat_obj,
  trs_vec      = v,
  celltype_col = celltype_col,
  rbp_label    = rbp,
  sort_by = c("median"),
  pal_cols     = rev(col_m)  # 图2那种浅->深
)


###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format


ggsave(
  filename = "./CELF4_5UTR_ASD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "./DDX24_5UTR_ASD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

 
 ggsave(
  filename = "./FTO_5UTR_ASD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

  ggsave(
  filename = "./ILF3_5UTR_ASD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

 
  ggsave(
  filename = "./RBMX_5UTR_ASD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


```













# ######################

# #-----ASD-CDS---FeaturePlot

```R
#-----2026-2-10


library(Seurat)
library(ggplot2)
library(grid)  # 用于 unit()
library(SeuratDisk)

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")



#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/01_ASD")

#----------ASD CDS
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_ASD_CDS <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/01_ASD/z_ASD_TRS/ASD_processed_magma_results_590Kcells_CDS_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


common <- intersect(Cells(seurat_obj), rownames(trs_ASD_CDS))
seurat_obj  <- subset(seurat_obj, cells = common)

#ASD CDS

rbp <- "DHX9"
 

rbp <- "ELAVL3"


rbp <- "HNRNPM"







v <- trs_ASD_CDS[common, rbp, drop = TRUE]

# 2) （可选）分位数截断，避免极端值把色标拉爆
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])

# 3) 写进 meta.data 并画图（用 umap.rna 做坐标）
seurat_obj[[paste0(rbp, "_TRS")]] <- v


# 更细腻的灰→白→橙红渐变（你也可以改成你喜欢的两个/三个端点）
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)

p <- FeaturePlot(
  seurat_obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap_continuum",
  pt.size   = 0.3,        # 你可以直接 4~6 试
  order     = TRUE,
  raster    = FALSE     # 关键：关闭栅格化
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,
      raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP"))



 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）

ggsave(
  "./DHX9_CDS_ASD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)





ggsave(
  "./ELAVL3_CDS_ASD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)



ggsave(
  "./HNRNPM_CDS_ASD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)



 



```



## ASD-CDS boxplot

```bash
#---2026-2-10

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/01_ASD")


#---------ASD CDS
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_ASD_CDS <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/01_ASD/z_ASD_TRS/ASD_processed_magma_results_590Kcells_CDS_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


#ASD CDS

rbp <- "DHX9"
 

rbp <- "ELAVL3"


rbp <- "HNRNPM"


 
 

###-Step 1 function: plot_trs_box_by_celltype()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

plot_trs_box_by_celltype <- function(seurat_obj,
                                     trs_vec,
                                     celltype_col,
                                     rbp_label = "RBP",
                                     pal_cols = c("#FCEFEA", "#F7C4B2", "#F26B4F"),
                                     ylab = "scRBP TRS",
                                     sort_by = c("median","mean"),
                                     show_outliers = FALSE,
                                     box_width = 0.65) {

  sort_by <- match.arg(sort_by)

  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  df <- data.frame(
    cell      = names(trs_vec),
    TRS       = as.numeric(trs_vec),
    cell_type = seurat_obj@meta.data[names(trs_vec), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      n = dplyr::n(),
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # 颜色渐进：按排序后的 cell type 映射到渐变调色板
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp_label, " TRS across cell types"),
      x = NULL, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}


###-----Step 1: Processing data
# v: 保证顺序=Seurat cells，并且带names
v <- trs_ASD_CDS[Cells(seurat_obj), rbp, drop = TRUE]
names(v) <- Cells(seurat_obj)

# 截断（不丢names的写法）
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
names(v) <- Cells(seurat_obj)   # 再保险补一次


###-----Step 2: run function
celltype_col <- "Type.v2"   # 你实际用的列名
 
col_m     = c("#FCEFEA", "#F7C4B2", "#F26B4F")
p_box <- plot_trs_box_by_celltype(
  seurat_obj   = seurat_obj,
  trs_vec      = v,
  celltype_col = celltype_col,
  rbp_label    = rbp,
  sort_by = c("median"),
  pal_cols     = rev(col_m)  # 图2那种浅->深
)


###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format


ggsave(
  filename = "./DHX9_CDS_ASD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)



ggsave(
  filename = "./ELAVL3_CDS_ASD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)



  ggsave(
  filename = "./HNRNPM_CDS_ASD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


```

















# ######################



# FeaturePlot---展示显著的RBP-UMAP---590K---ADHD



# ######################

# #-----ADHD-3UTR---FeaturePlot

```R
#-----2026-2-10


library(Seurat)
library(ggplot2)
library(grid)  # 用于 unit()
library(SeuratDisk)

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")



#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/02_ADHD")

#----------ADHD 3UTR
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_ADHD_3UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/02_ADHD/z_ADHD_TRS/ADHD_processed_magma_results_590Kcells_3UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


common <- intersect(Cells(seurat_obj), rownames(trs_ADHD_3UTR))
seurat_obj  <- subset(seurat_obj, cells = common)

#ADHD 3UTR

rbp <- "CELF2"
 


rbp <- "SRRM4"

rbp <- "RBFOX2"







v <- trs_ADHD_3UTR[common, rbp, drop = TRUE]

# 2) （可选）分位数截断，避免极端值把色标拉爆
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])

# 3) 写进 meta.data 并画图（用 umap.rna 做坐标）
seurat_obj[[paste0(rbp, "_TRS")]] <- v


# 更细腻的灰→白→橙红渐变（你也可以改成你喜欢的两个/三个端点）
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)

p <- FeaturePlot(
  seurat_obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap_continuum",
  pt.size   = 0.3,        # 你可以直接 4~6 试
  order     = TRUE,
  raster    = FALSE     # 关键：关闭栅格化
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,
      raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP"))



 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）

ggsave(
  "./CELF2_3UTR_ADHD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./SRRM4_3UTR_ADHD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)




ggsave(
  "./RBFOX2_3UTR_ADHD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)






```



## ADHD-3UTR boxplot

```bash
#---2026-2-10

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/02_ADHD")


#---------ADHD 3UTR
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_ADHD_3UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/02_ADHD/z_ADHD_TRS/ADHD_processed_magma_results_590Kcells_3UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


#ADHD 3UTR
rbp <- "CELF2"
 


rbp <- "SRRM4"

rbp <- "RBFOX2"


 
 

###-Step 1 function: plot_trs_box_by_celltype()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

plot_trs_box_by_celltype <- function(seurat_obj,
                                     trs_vec,
                                     celltype_col,
                                     rbp_label = "RBP",
                                     pal_cols = c("#FCEFEA", "#F7C4B2", "#F26B4F"),
                                     ylab = "scRBP TRS",
                                     sort_by = c("median","mean"),
                                     show_outliers = FALSE,
                                     box_width = 0.65) {

  sort_by <- match.arg(sort_by)

  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  df <- data.frame(
    cell      = names(trs_vec),
    TRS       = as.numeric(trs_vec),
    cell_type = seurat_obj@meta.data[names(trs_vec), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      n = dplyr::n(),
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # 颜色渐进：按排序后的 cell type 映射到渐变调色板
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp_label, " TRS across cell types"),
      x = NULL, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}


###-----Step 1: Processing data
# v: 保证顺序=Seurat cells，并且带names
v <- trs_ADHD_3UTR[Cells(seurat_obj), rbp, drop = TRUE]
names(v) <- Cells(seurat_obj)

# 截断（不丢names的写法）
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
names(v) <- Cells(seurat_obj)   # 再保险补一次


###-----Step 2: run function
celltype_col <- "Type.v2"   # 你实际用的列名
 
col_m     = c("#FCEFEA", "#F7C4B2", "#F26B4F")
p_box <- plot_trs_box_by_celltype(
  seurat_obj   = seurat_obj,
  trs_vec      = v,
  celltype_col = celltype_col,
  rbp_label    = rbp,
  sort_by = c("median"),
  pal_cols     = rev(col_m)  # 图2那种浅->深
)


###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format


ggsave(
  filename = "./CELF2_3UTR_ADHD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "./SRRM4_3UTR_ADHD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

 
 ggsave(
  filename = "./RBFOX2_3UTR_ADHD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)



```









# #-----ADHD-5UTR---FeaturePlot

```R
#-----2026-2-10


library(Seurat)
library(ggplot2)
library(grid)  # 用于 unit()
library(SeuratDisk)

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")



#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/02_ADHD")

#----------ADHD 5UTR
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_ADHD_5UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/02_ADHD/z_ADHD_TRS/ADHD_processed_magma_results_590Kcells_5UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


common <- intersect(Cells(seurat_obj), rownames(trs_ADHD_5UTR))
seurat_obj  <- subset(seurat_obj, cells = common)

 
rbp <- "CELF4"
 


rbp <- "HNRNPH1"





rbp <- "RPS19"

rbp <- "EIF4G2"

rbp <- "RBFOX2"





v <- trs_ADHD_5UTR[common, rbp, drop = TRUE]

# 2) （可选）分位数截断，避免极端值把色标拉爆
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])

# 3) 写进 meta.data 并画图（用 umap.rna 做坐标）
seurat_obj[[paste0(rbp, "_TRS")]] <- v


# 更细腻的灰→白→橙红渐变（你也可以改成你喜欢的两个/三个端点）
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)

p <- FeaturePlot(
  seurat_obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap_continuum",
  pt.size   = 0.3,        # 你可以直接 4~6 试
  order     = TRUE,
  raster    = FALSE     # 关键：关闭栅格化
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,
      raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP"))



 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）

ggsave(
  "./CELF4_5UTR_ADHD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./HNRNPH1_5UTR_ADHD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)




ggsave(
  "./RPS19_5UTR_ADHD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)



ggsave(
  "./EIF4G2_5UTR_ADHD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)


ggsave(
  "./RBFOX2_5UTR_ADHD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)


```



## ADHD-5UTR boxplot

```bash
#---2026-2-10

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/02_ADHD")


#---------ADHD 
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_ADHD_5UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/02_ADHD/z_ADHD_TRS/ADHD_processed_magma_results_590Kcells_5UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


#ADHD  
rbp <- "CELF4"
 


rbp <- "HNRNPH1"





rbp <- "RPS19"

rbp <- "EIF4G2"

rbp <- "RBFOX2"




 
 

###-Step 1 function: plot_trs_box_by_celltype()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

plot_trs_box_by_celltype <- function(seurat_obj,
                                     trs_vec,
                                     celltype_col,
                                     rbp_label = "RBP",
                                     pal_cols = c("#FCEFEA", "#F7C4B2", "#F26B4F"),
                                     ylab = "scRBP TRS",
                                     sort_by = c("median","mean"),
                                     show_outliers = FALSE,
                                     box_width = 0.65) {

  sort_by <- match.arg(sort_by)

  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  df <- data.frame(
    cell      = names(trs_vec),
    TRS       = as.numeric(trs_vec),
    cell_type = seurat_obj@meta.data[names(trs_vec), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      n = dplyr::n(),
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # 颜色渐进：按排序后的 cell type 映射到渐变调色板
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp_label, " TRS across cell types"),
      x = NULL, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}


###-----Step 1: Processing data
# v: 保证顺序=Seurat cells，并且带names
v <- trs_ADHD_5UTR[Cells(seurat_obj), rbp, drop = TRUE]
names(v) <- Cells(seurat_obj)

# 截断（不丢names的写法）
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
names(v) <- Cells(seurat_obj)   # 再保险补一次


###-----Step 2: run function
celltype_col <- "Type.v2"   # 你实际用的列名
 
col_m     = c("#FCEFEA", "#F7C4B2", "#F26B4F")
p_box <- plot_trs_box_by_celltype(
  seurat_obj   = seurat_obj,
  trs_vec      = v,
  celltype_col = celltype_col,
  rbp_label    = rbp,
  sort_by = c("median"),
  pal_cols     = rev(col_m)  # 图2那种浅->深
)


###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format


ggsave(
  filename = "./CELF4_5UTR_ADHD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "./HNRNPH1_5UTR_ADHD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

 
 ggsave(
  filename = "./RPS19_5UTR_ADHD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


  
 ggsave(
  filename = "./EIF4G2_5UTR_ADHD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)



 ggsave(
  filename = "./RBFOX2_5UTR_ADHD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

```















# #-----ADHD-Introns---FeaturePlot

```R
#-----2026-2-10


library(Seurat)
library(ggplot2)
library(grid)  # 用于 unit()
library(SeuratDisk)

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")



#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/02_ADHD")

#----------ADHD Introns
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_ADHD_Introns <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/02_ADHD/z_ADHD_TRS/ADHD_processed_magma_results_590Kcells_Introns_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


common <- intersect(Cells(seurat_obj), rownames(trs_ADHD_Introns))
seurat_obj  <- subset(seurat_obj, cells = common)

 
rbp <- "CELF2"
 


rbp <- "RPS27A"







v <- trs_ADHD_Introns[common, rbp, drop = TRUE]

# 2) （可选）分位数截断，避免极端值把色标拉爆
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])

# 3) 写进 meta.data 并画图（用 umap.rna 做坐标）
seurat_obj[[paste0(rbp, "_TRS")]] <- v


# 更细腻的灰→白→橙红渐变（你也可以改成你喜欢的两个/三个端点）
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)

p <- FeaturePlot(
  seurat_obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap_continuum",
  pt.size   = 0.3,        # 你可以直接 4~6 试
  order     = TRUE,
  raster    = FALSE     # 关键：关闭栅格化
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,
      raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP"))



 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）

ggsave(
  "./CELF2_Introns_ADHD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./RPS27A_Introns_ADHD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)
 





```









## ADHD-Introns boxplot

```bash
#---2026-2-10

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/02_ADHD")


#---------ADHD 
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_ADHD_Introns <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/02_ADHD/z_ADHD_TRS/ADHD_processed_magma_results_590Kcells_Introns_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


#ADHD  
rbp <- "CELF2"
 


rbp <- "RPS27A"




 
 

###-Step 1 function: plot_trs_box_by_celltype()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

plot_trs_box_by_celltype <- function(seurat_obj,
                                     trs_vec,
                                     celltype_col,
                                     rbp_label = "RBP",
                                     pal_cols = c("#FCEFEA", "#F7C4B2", "#F26B4F"),
                                     ylab = "scRBP TRS",
                                     sort_by = c("median","mean"),
                                     show_outliers = FALSE,
                                     box_width = 0.65) {

  sort_by <- match.arg(sort_by)

  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  df <- data.frame(
    cell      = names(trs_vec),
    TRS       = as.numeric(trs_vec),
    cell_type = seurat_obj@meta.data[names(trs_vec), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      n = dplyr::n(),
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # 颜色渐进：按排序后的 cell type 映射到渐变调色板
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp_label, " TRS across cell types"),
      x = NULL, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}


###-----Step 1: Processing data
# v: 保证顺序=Seurat cells，并且带names
v <- trs_ADHD_Introns[Cells(seurat_obj), rbp, drop = TRUE]
names(v) <- Cells(seurat_obj)

# 截断（不丢names的写法）
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
names(v) <- Cells(seurat_obj)   # 再保险补一次


###-----Step 2: run function
celltype_col <- "Type.v2"   # 你实际用的列名
 
col_m     = c("#FCEFEA", "#F7C4B2", "#F26B4F")
p_box <- plot_trs_box_by_celltype(
  seurat_obj   = seurat_obj,
  trs_vec      = v,
  celltype_col = celltype_col,
  rbp_label    = rbp,
  sort_by = c("median"),
  pal_cols     = rev(col_m)  # 图2那种浅->深
)


###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format


ggsave(
  filename = "./CELF2_Introns_ADHD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "./RPS27A_Introns_ADHD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

 


```













# ######################



# FeaturePlot---展示显著的RBP-UMAP---590K---TS



# ######################

# #-----TS-3UTR---FeaturePlot

```R
#-----2026-2-10


library(Seurat)
library(ggplot2)
library(grid)  # 用于 unit()
library(SeuratDisk)

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")



#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/08_TS")

#----------TS
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_TS_3UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/08_TS/z_TS_TRS/TS_processed_magma_results_590Kcells_3UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


common <- intersect(Cells(seurat_obj), rownames(trs_TS_3UTR))
seurat_obj  <- subset(seurat_obj, cells = common)

#TS

rbp <- "ELAVL3"
 


rbp <- "ELAVL4"

rbp <- "CPEB4"




v <- trs_TS_3UTR[common, rbp, drop = TRUE]

# 2) （可选）分位数截断，避免极端值把色标拉爆
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])

# 3) 写进 meta.data 并画图（用 umap.rna 做坐标）
seurat_obj[[paste0(rbp, "_TRS")]] <- v


# 更细腻的灰→白→橙红渐变（你也可以改成你喜欢的两个/三个端点）
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)

p <- FeaturePlot(
  seurat_obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap_continuum",
  pt.size   = 0.3,        # 你可以直接 4~6 试
  order     = TRUE,
  raster    = FALSE     # 关键：关闭栅格化
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,
      raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP"))



 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）

ggsave(
  "./ELAVL3_3UTR_TS.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./ELAVL4_3UTR_TS.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)




ggsave(
  "./CPEB4_3UTR_TS.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)






```



## TS-3UTR boxplot

```bash
#---2026-2-10

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/08_TS")


#---------TS 3UTR
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_TS_3UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/08_TS/z_TS_TRS/TS_processed_magma_results_590Kcells_3UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


#TS
rbp <- "ELAVL3"
 


rbp <- "ELAVL4"

rbp <- "CPEB4"


 
 

###-Step 1 function: plot_trs_box_by_celltype()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

plot_trs_box_by_celltype <- function(seurat_obj,
                                     trs_vec,
                                     celltype_col,
                                     rbp_label = "RBP",
                                     pal_cols = c("#FCEFEA", "#F7C4B2", "#F26B4F"),
                                     ylab = "scRBP TRS",
                                     sort_by = c("median","mean"),
                                     show_outliers = FALSE,
                                     box_width = 0.65) {

  sort_by <- match.arg(sort_by)

  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  df <- data.frame(
    cell      = names(trs_vec),
    TRS       = as.numeric(trs_vec),
    cell_type = seurat_obj@meta.data[names(trs_vec), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      n = dplyr::n(),
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # 颜色渐进：按排序后的 cell type 映射到渐变调色板
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp_label, " TRS across cell types"),
      x = NULL, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}


###-----Step 1: Processing data
# v: 保证顺序=Seurat cells，并且带names
v <- trs_TS_3UTR[Cells(seurat_obj), rbp, drop = TRUE]
names(v) <- Cells(seurat_obj)

# 截断（不丢names的写法）
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
names(v) <- Cells(seurat_obj)   # 再保险补一次


###-----Step 2: run function
celltype_col <- "Type.v2"   # 你实际用的列名
 
col_m     = c("#FCEFEA", "#F7C4B2", "#F26B4F")
p_box <- plot_trs_box_by_celltype(
  seurat_obj   = seurat_obj,
  trs_vec      = v,
  celltype_col = celltype_col,
  rbp_label    = rbp,
  sort_by = c("median"),
  pal_cols     = rev(col_m)  # 图2那种浅->深
)


###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format


ggsave(
  filename = "./ELAVL3_3UTR_TS.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "./ELAVL4_3UTR_TS.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


 ggsave(
  filename = "./CPEB4_3UTR_TS.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)



```















# ######################



# FeaturePlot---展示显著的RBP-UMAP---590K---AN



# ######################

# #-----AN-5UTR---FeaturePlot

```R
#-----2026-2-10


library(Seurat)
library(ggplot2)
library(grid)  # 用于 unit()
library(SeuratDisk)

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")



#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/03_AN")

#----------AN
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_AN_5UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/03_AN/z_AN_TRS/AN_processed_magma_results_590Kcells_5UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


common <- intersect(Cells(seurat_obj), rownames(trs_AN_5UTR))
seurat_obj  <- subset(seurat_obj, cells = common)

#AN

rbp <- "HNRNPA1"
 


rbp <- "HNRNPH1"




v <- trs_AN_5UTR[common, rbp, drop = TRUE]

# 2) （可选）分位数截断，避免极端值把色标拉爆
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])

# 3) 写进 meta.data 并画图（用 umap.rna 做坐标）
seurat_obj[[paste0(rbp, "_TRS")]] <- v


# 更细腻的灰→白→橙红渐变（你也可以改成你喜欢的两个/三个端点）
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)

p <- FeaturePlot(
  seurat_obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap_continuum",
  pt.size   = 0.3,        # 你可以直接 4~6 试
  order     = TRUE,
  raster    = FALSE     # 关键：关闭栅格化
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,
      raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP"))



 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）

ggsave(
  "./HNRNPA1_5UTR_AN.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./HNRNPH1_5UTR_AN.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)



```



## AN-5UTR boxplot

```bash
#---2026-2-10

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/08_TS")


#---------AN 5UTR
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_AN_5UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/03_AN/z_AN_TRS/AN_processed_magma_results_590Kcells_5UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


#TS
rbp <- "HNRNPA1"
 


rbp <- "HNRNPH1"



 
 

###-Step 1 function: plot_trs_box_by_celltype()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

plot_trs_box_by_celltype <- function(seurat_obj,
                                     trs_vec,
                                     celltype_col,
                                     rbp_label = "RBP",
                                     pal_cols = c("#FCEFEA", "#F7C4B2", "#F26B4F"),
                                     ylab = "scRBP TRS",
                                     sort_by = c("median","mean"),
                                     show_outliers = FALSE,
                                     box_width = 0.65) {

  sort_by <- match.arg(sort_by)

  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  df <- data.frame(
    cell      = names(trs_vec),
    TRS       = as.numeric(trs_vec),
    cell_type = seurat_obj@meta.data[names(trs_vec), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      n = dplyr::n(),
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # 颜色渐进：按排序后的 cell type 映射到渐变调色板
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp_label, " TRS across cell types"),
      x = NULL, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}


###-----Step 1: Processing data
# v: 保证顺序=Seurat cells，并且带names
v <- trs_AN_5UTR[Cells(seurat_obj), rbp, drop = TRUE]
names(v) <- Cells(seurat_obj)

# 截断（不丢names的写法）
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
names(v) <- Cells(seurat_obj)   # 再保险补一次


###-----Step 2: run function
celltype_col <- "Type.v2"   # 你实际用的列名
 
col_m     = c("#FCEFEA", "#F7C4B2", "#F26B4F")
p_box <- plot_trs_box_by_celltype(
  seurat_obj   = seurat_obj,
  trs_vec      = v,
  celltype_col = celltype_col,
  rbp_label    = rbp,
  sort_by = c("median"),
  pal_cols     = rev(col_m)  # 图2那种浅->深
)


###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format

 
ggsave(
  filename = "./HNRNPA1_5UTR_AN.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


 ggsave(
  filename = "./HNRNPH1_5UTR_AN.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

 

```











# ######################



# FeaturePlot---展示显著的RBP-UMAP---590K---MDD



# ######################

# #-----MDD-5UTR---FeaturePlot

```R
#-----2026-2-10


library(Seurat)
library(ggplot2)
library(grid)  # 用于 unit()
library(SeuratDisk)

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")



#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/05_MDD")

#----------MDD
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_MDD_5UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/05_MDD/z_MDD_TRS/MDD_processed_magma_results_590Kcells_5UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


common <- intersect(Cells(seurat_obj), rownames(trs_MDD_5UTR))
seurat_obj  <- subset(seurat_obj, cells = common)

#MDD

rbp <- "ERI2"
 


rbp <- "CELF4"




v <- trs_MDD_5UTR[common, rbp, drop = TRUE]

# 2) （可选）分位数截断，避免极端值把色标拉爆
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])

# 3) 写进 meta.data 并画图（用 umap.rna 做坐标）
seurat_obj[[paste0(rbp, "_TRS")]] <- v


# 更细腻的灰→白→橙红渐变（你也可以改成你喜欢的两个/三个端点）
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)

p <- FeaturePlot(
  seurat_obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap_continuum",
  pt.size   = 0.3,        # 你可以直接 4~6 试
  order     = TRUE,
  raster    = FALSE     # 关键：关闭栅格化
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,
      raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP"))



 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）

ggsave(
  "./ERI2_5UTR_MDD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./CELF4_5UTR_MDD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)



```



## MDD-5UTR boxplot

```bash
#---2026-2-10

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/05_MDD")


#---------MDD
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_MDD_5UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/05_MDD/z_MDD_TRS/MDD_processed_magma_results_590Kcells_5UTR_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


#MDD
rbp <- "ERI2"
 


rbp <- "CELF4"



 
 

###-Step 1 function: plot_trs_box_by_celltype()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

plot_trs_box_by_celltype <- function(seurat_obj,
                                     trs_vec,
                                     celltype_col,
                                     rbp_label = "RBP",
                                     pal_cols = c("#FCEFEA", "#F7C4B2", "#F26B4F"),
                                     ylab = "scRBP TRS",
                                     sort_by = c("median","mean"),
                                     show_outliers = FALSE,
                                     box_width = 0.65) {

  sort_by <- match.arg(sort_by)

  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  df <- data.frame(
    cell      = names(trs_vec),
    TRS       = as.numeric(trs_vec),
    cell_type = seurat_obj@meta.data[names(trs_vec), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      n = dplyr::n(),
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # 颜色渐进：按排序后的 cell type 映射到渐变调色板
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp_label, " TRS across cell types"),
      x = NULL, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}


###-----Step 1: Processing data
# v: 保证顺序=Seurat cells，并且带names
v <- trs_MDD_5UTR[Cells(seurat_obj), rbp, drop = TRUE]
names(v) <- Cells(seurat_obj)

# 截断（不丢names的写法）
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
names(v) <- Cells(seurat_obj)   # 再保险补一次


###-----Step 2: run function
celltype_col <- "Type.v2"   # 你实际用的列名
 
col_m     = c("#FCEFEA", "#F7C4B2", "#F26B4F")
p_box <- plot_trs_box_by_celltype(
  seurat_obj   = seurat_obj,
  trs_vec      = v,
  celltype_col = celltype_col,
  rbp_label    = rbp,
  sort_by = c("median"),
  pal_cols     = rev(col_m)  # 图2那种浅->深
)


###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format

 
ggsave(
  filename = "./ERI2_5UTR_MDD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)


 ggsave(
  filename = "./CELF4_5UTR_MDD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

 

```











# ######################



# FeaturePlot---展示显著的RBP-UMAP---590K---OCD



# ######################

# #-----OCD-Introns--FeaturePlot

```R
#-----2026-2-10


library(Seurat)
library(ggplot2)
library(grid)  # 用于 unit()
library(SeuratDisk)

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")



#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/06_OCD")

#----------OCD
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_OCD_Introns <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/06_OCD/z_OCD_TRS/OCD_processed_magma_results_590Kcells_Introns_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


common <- intersect(Cells(seurat_obj), rownames(trs_OCD_Introns))
seurat_obj  <- subset(seurat_obj, cells = common)

#Introns

rbp <- "LSM11"
rbp <- "RBFOX2"
rbp <- "RBFOX3"




v <- trs_OCD_Introns[common, rbp, drop = TRUE]

# 2) （可选）分位数截断，避免极端值把色标拉爆
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])

# 3) 写进 meta.data 并画图（用 umap.rna 做坐标）
seurat_obj[[paste0(rbp, "_TRS")]] <- v


# 更细腻的灰→白→橙红渐变（你也可以改成你喜欢的两个/三个端点）
pal <- colorRampPalette(c("grey85", "white", "#F87C56"))(256)

p <- FeaturePlot(
  seurat_obj,
  features  = paste0(rbp, "_TRS"),
  reduction = "umap_continuum",
  pt.size   = 0.3,        # 你可以直接 4~6 试
  order     = TRUE,
  raster    = FALSE     # 关键：关闭栅格化
) +
  scale_color_gradientn(
    colours = pal,
    name    = "TRS",
    oob     = scales::squish,
    guide   = guide_colorbar(
      nbin = 256,
      raster = TRUE,
      barheight = unit(35, "mm"),
      barwidth  = unit(4,  "mm")
    )
  ) +
  theme_void() +
  ggtitle(paste(rbp, "TRS on UMAP"))



 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）
 # 或者 TIFF（期刊更喜欢）

ggsave(
  "./LSM11_Introns_OCD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  "./RBFOX2_Introns_OCD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  "./RBFOX3_Introns_OCD.umap_trs3.tiff",
  p,
  width = 7, height = 5, units = "in",
  dpi = 300,
  compression = "lzw"
)

```



## OCD-Introns boxplot

```bash
#---2026-2-10

# 读取你的 Seurat 对象
seurat_obj <- readRDS("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/clean.Dev_metaatlas_integrated_65580genesx599211cells_recluster_umap_v2.rds")

#location
setwd("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/06_OCD")


#---------OCD
# 1) 读 TRS 矩阵（行=cell，列=RBP）；保留对象里的细胞并按 Seurat 顺序排列
trs_MDD_5UTR <- read.csv("/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/replication_2_developing_brain_cells/z_seed50times_100K_scRBP/zzz_summary/scRBP_trs_8brain_traits_single_cellModes/06_OCD/z_OCD_TRS/OCD_processed_magma_results_590Kcells_Introns_REAL_only.sc.TRS_matrix.csv", row.names = 1, check.names = FALSE)


#OCD
rbp <- "LSM11"
rbp <- "RBFOX2"
rbp <- "RBFOX3"

 
 

###-Step 1 function: plot_trs_box_by_celltype()
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

plot_trs_box_by_celltype <- function(seurat_obj,
                                     trs_vec,
                                     celltype_col,
                                     rbp_label = "RBP",
                                     pal_cols = c("#FCEFEA", "#F7C4B2", "#F26B4F"),
                                     ylab = "scRBP TRS",
                                     sort_by = c("median","mean"),
                                     show_outliers = FALSE,
                                     box_width = 0.65) {

  sort_by <- match.arg(sort_by)

  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))

  df <- data.frame(
    cell      = names(trs_vec),
    TRS       = as.numeric(trs_vec),
    cell_type = seurat_obj@meta.data[names(trs_vec), celltype_col, drop = TRUE],
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(TRS), !is.na(cell_type))

  stat_df <- df |>
    group_by(cell_type) |>
    summarise(
      n = dplyr::n(),
      center = if (sort_by == "median") median(TRS) else mean(TRS),
      .groups = "drop"
    ) |>
    arrange(desc(center))

  df$cell_type <- factor(df$cell_type, levels = stat_df$cell_type)

  # 颜色渐进：按排序后的 cell type 映射到渐变调色板
  pal_n <- colorRampPalette(pal_cols)(nlevels(df$cell_type))
  names(pal_n) <- levels(df$cell_type)

  ggplot(df, aes(x = cell_type, y = TRS, fill = cell_type)) +
    geom_boxplot(
      width = box_width,
      outlier.shape = if (show_outliers) 16 else NA,
      size = 0.25
    ) +
    scale_fill_manual(values = pal_n, guide = "none") +
    labs(
      title = paste0(rbp_label, " TRS across cell types"),
      x = NULL, y = ylab
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      plot.title  = element_text(hjust = 0.5)
    )
}


###-----Step 1: Processing data
# v: 保证顺序=Seurat cells，并且带names
v <- trs_OCD_Introns[Cells(seurat_obj), rbp, drop = TRUE]
names(v) <- Cells(seurat_obj)

# 截断（不丢names的写法）
q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
v <- pmin(pmax(v, q[1]), q[2])
names(v) <- Cells(seurat_obj)   # 再保险补一次


###-----Step 2: run function
celltype_col <- "Type.v2"   # 你实际用的列名
 
col_m     = c("#FCEFEA", "#F7C4B2", "#F26B4F")
p_box <- plot_trs_box_by_celltype(
  seurat_obj   = seurat_obj,
  trs_vec      = v,
  celltype_col = celltype_col,
  rbp_label    = rbp,
  sort_by = c("median"),
  pal_cols     = rev(col_m)  # 图2那种浅->深
)


###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format
###------Step 3: save .tiff format

 
ggsave(
  filename = "./LSM11_Introns_OCD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)



 ggsave(
  filename = "./RBFOX2_Introns_OCD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)

 
 ggsave(
  filename = "./RBFOX3_Introns_OCD.TRS_box_by_celltype.tiff",
  plot = p_box,
  width = 10, height = 4.2, units = "in",
  dpi = 600,
  compression = "lzw"
)
```



