
##--R Seurat---> plot UMAP
library(zellkonverter)
library(SingleCellExperiment)
library(Seurat)
library(Matrix)

# 1) 读 h5ad（不依赖 reticulate）
sce <- readH5AD("PBMC_healthy_subset.h5ad", use_hdf5 = TRUE)

# 2) 如果只有 "X" 没有 "counts"，复制一份作为 counts
if (!"counts" %in% assayNames(sce)) {
  assay(sce, "counts") <- assay(sce, "X")
}

# 3) 确保是稀疏 dgCMatrix（内存友好；Seurat更稳）
assay(sce, "counts") <- as(assay(sce, "counts"), "dgCMatrix")
assay(sce, "X")      <- as(assay(sce, "X"),      "dgCMatrix")

# 4) 行名/列名唯一
rownames(sce) <- make.unique(rownames(sce))
colnames(sce) <- make.unique(colnames(sce))

# 5) 转成 Seurat：显式指定 counts 和 data
obj <- as.Seurat(sce, counts = "counts", data = "X")
DefaultAssay(obj) <- "originalexp"
obj <- RenameAssays(obj, originalexp = "RNA")
DefaultAssay(obj) <- "RNA"

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)

# 1) 可变基因
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
obj <- ScaleData(obj)
# 2) PCA（v5 支持分层，这里默认用 data 层）
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

# 3) 邻接图 & UMAP
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj, resolution = 1.0)
obj <- RunUMAP(obj, dims = 1:50)


# 4) 画图 + 保存 PDF
p <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.3)
ggplot2::ggsave("PBMC_umap_clusters.pdf", p, width = 6, height = 5)

p <- DimPlot(obj, reduction = "umap", group.by = "annotation_broad", pt.size = 0.3)
ggplot2::ggsave("PBMC_umap_annotation_broad.pdf", p, width = 8, height = 8)

p <- DimPlot(obj, reduction = "umap", group.by = "annotation_detailed", pt.size = 0.3)
ggplot2::ggsave("PBMC_umap_annotation_detailed.pdf", p, width = 12, height =12)
 

  
  
  
library(Seurat)
library(ggplot2)
library(scales)  # hue_pal() 用来生成明亮颜色

# 用 annotation_broad 作为身份（可选）
Idents(obj) <- "annotation_broad"

# 明亮配色：用高饱和度的色相轮
lvls <- levels(obj$annotation_broad)              # 所有细胞类型
pal  <- hue_pal(l = 60, c = 100)(length(lvls))    # 亮一些（提高 chroma）

# 画图 + 文字标签
p <- DimPlot(
  obj,
  reduction = "umap",
  group.by  = "annotation_broad",
  cols      = pal,
  label     = TRUE,       # 在簇中心加标签
  repel     = TRUE,       # 标签自动避让
  label.size= 4,
  pt.size   = 1
) + theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 16)
    )

ggsave("PBMC_umap_annotation_broad_labeled.pdf", p, width = 8, height = 8)
  

  
p <- DimPlot(
  obj,
  reduction = "umap",
  group.by  = "annotation_broad",
  cols      = pal,
  label     = F,       # 在簇中心加标签
  repel     = TRUE,       # 标签自动避让
  label.size= 4,
  pt.size   = 1
) + theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 16)
    )

  
  
# 需要 ragg，若未安装：install.packages("ragg")
ggsave(
  "PBMC_umap_annotation_broad_nonlabeled.png",
  p,
  device = ragg::agg_png,   # 抗锯齿更好
  width  = 8, height = 8, units = "in",
  dpi    = 600, bg = "white"
)
  
 



p <- DimPlot(
  obj,
  reduction = "umap",
  group.by  = "annotation_detailed",
  cols      = pal,
  label     = T,       # 在簇中心加标签
  repel     = TRUE,       # 标签自动避让
  label.size= 4,
  pt.size   = 1
) + theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 16)
    )

  
  
# 需要 ragg，若未安装：install.packages("ragg")
ggsave(
  "PBMC_umap_annotation_detailed_nonlabeled.png",
  p,
  device = ragg::agg_png,   # 抗锯齿更好
  width  = 8, height = 8, units = "in",
  dpi    = 600, bg = "white"
)
  
  
  
p <- DimPlot(
  obj,
  reduction = "umap",
  group.by  = "annotation_detailed",
  cols      = pal,
  label     = F,       # 在簇中心加标签
  repel     = TRUE,       # 标签自动避让
  label.size= 4,
  pt.size   = 1
) + theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 16)
    )

  
  
# 需要 ragg，若未安装：install.packages("ragg")
ggsave(
  "PBMC_umap_annotation_detailed_nonlabeled.png",
  p,
  device = ragg::agg_png,   # 抗锯齿更好
  width  = 8, height = 8, units = "in",
  dpi    = 600, bg = "white"
)  
