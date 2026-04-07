#2025-08-06
##------All 8900 PBMC cells
library(Seurat)
df <- readRDS("10X_PBMC.rds")
# Extract the data matrix from the RNA assay
rna_matrix <- GetAssayData(df, assay = "RNA", layer = "counts")
#rna_matrix <- df@assays$RNA@counts
# Save as CSV file
write.csv(as.matrix(rna_matrix), file = "10X_PBMC_RNA_matrix.csv",quote=F, row.names = TRUE)
##------All 3313 monocytes vs 3280 T cells PBMC cells
##------All 3313 monocytes vs 3280 T+cells PBMC cells
##------All 3313 monocytes vs 3280 T+cells PBMC cells
##------All 3313 monocytes vs 3280 T+cells PBMC cells
library(Seurat)
single_cell_8900 <- readRDS("10X_PBMC.rds")
Idents(single_cell_8900) <- single_cell_8900$cell_type
#single_cell$cell_type2[which(single_cell$cell_type2 %in% c("CD14+ Monocytes","FCGR3A+ Monocytes"))] <- "Monocytes"
cell_type2<- as.character(single_cell_8900$cell_type)
cell_type2[which(cell_type2 %in% c("CD14+ Monocytes","FCGR3A+ Monocytes"))] <- "Monocytes"
cell_type2[which(cell_type2 %in% c("CD4+ T cells","CD8+ T cells"))] <- "T_cells"
single_cell_8900$cell_type2 <- cell_type2
Idents(single_cell_8900) <- single_cell_8900$cell_type2 
table(Idents(single_cell_8900))
#--Extract subset:
subset_cells <- subset(single_cell_8900, idents = c("T_cells", "Monocytes"))
table(Idents(subset_cells))
##Re-clustering analysis:
DefaultAssay(subset_cells)<-"RNA"
subset_cells <- NormalizeData(subset_cells, normalization.method = "LogNormalize", scale.factor = 10000)
subset_cells <- FindVariableFeatures(subset_cells, selection.method = "vst", nfeatures = 3000)
subset_cells <- ScaleData(subset_cells, verbose = FALSE)
subset_cells <- RunPCA(subset_cells, npcs = 20, verbose = FALSE)
subset_cells  <- FindNeighbors(subset_cells, dims = 1:10)
subset_cells  <- FindClusters(object = subset_cells, resolution = 0.2, graph.name="wsnn")
subset_cells <- RunUMAP(subset_cells, reduction = "pca", dims = 1:20,reduction.name = 'umap.rna', reduction.key = 'UMAP.RNA_')
pal <- c("#A13B46","#67ADB7")
DefaultAssay(subset_cells)<-"RNA"
pdf("subset_cells_cell_type2_6.5Kcells_v2.pdf")
DimPlot(subset_cells,reduction = "umap.rna",group.by = "cell_type2",label = T,cols=pal, pt.size = 0.1, repel = T)
dev.off()
##--subset
DefaultAssay(subset_cells)<-"RNA"
pdf("subset_cells_cell_type_6.5Kcells_v2.pdf")
DimPlot(subset_cells,reduction = "umap.rna",group.by = "cell_type",label = T,cols=pal, pt.size = 0.1, repel = T)
dev.off()
#saveRDS(subset_cells, file="10X_PBMC_subset_6.5K.rds")
###----Extract matrix for scRBP
subset_cells_matrix <- GetAssayData(subset_cells, assay = "RNA", layer = "counts")
 
# Save as CSV file
write.csv(as.matrix(subset_cells_matrix), file = "10X_PBMC_subset_6.5K.csv", quote=F, row.names = TRUE)
##_-----Save in feather format
library(arrow)
library(tibble)
# subset_cells_matrix: should be a gene x cell matrix
df <- as.data.frame(as.matrix(subset_cells_matrix))
# Write row names into the first column "gene"
df <- rownames_to_column(df, var = "gene")
# Optional: match case with RBP list (if RBPs are uppercase)
# df$gene <- toupper(df$gene)
# Save as feather (note: use write_feather, not write.feather)
write_feather(df, "10X_PBMC_subset_6.5K.feather")
