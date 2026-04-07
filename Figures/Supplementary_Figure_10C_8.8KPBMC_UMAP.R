#2025-08-07


##------All 10393 PBMC cells
library(Seurat)

df <- readRDS("pbmc_merged.rds")

DefaultAssay(df)<-"RNA"
# Extract the data matrix from the RNA assay
rna_matrix <- GetAssayData(df, assay = "RNA", layer = "counts")
#rna_matrix <- df@assays$RNA@counts

# Save as CSV file
write.csv(as.matrix(rna_matrix), file = "10X_pbmc_merged_10Kcells_matrix.csv",quote=F, row.names = TRUE)



#---8.8K----
##------All 3546 monocytes vs 5231 CD8T+cells PBMC cells
##------All 3546 monocytes vs 5231 CD8T+cells PBMC cells
##------All 3546 monocytes vs 5231 CD8T+cells PBMC cells
##------All 3546 monocytes vs 5231 CD8T+cells PBMC cells
library(Seurat)

single_cell_10K <- readRDS("pbmc_merged.rds")


Idents(single_cell_10K) <- single_cell_10K$merged_celltype
 

cell_type2<- as.character(single_cell_10K$merged_celltype)
#cell_type2[which(cell_type2 %in% c("CD14+ Monocytes","FCGR3A+ Monocytes"))] <- "Monocytes"
cell_type2[which(cell_type2 %in% c("CD4+ T cells","CD8+ T cells"))] <- "T_cells"
single_cell_10K$merged_celltype2 <- cell_type2



Idents(single_cell_10K) <- single_cell_10K$merged_celltype2 

table(Idents(single_cell_10K))


###_-----Annotation----new cell type subpopulations
single_cell_10K <- annotate_monocyte_subtypes(
  single_cell_10K, celltype_col="merged_celltype", monocyte_label="Monocytes",
  assay="SCT", delta=0.25, min_expr=0.1, out_col="merged_celltype_refined"
)
table(single_cell_10K$merged_celltype_refined)

cell_names <- single_cell_10K$merged_celltype_refined 
cell_chr   <- as.character(cell_names)
cell_chr[cell_chr == "Monocytes_ambiguous"] <- "CD14+ Monocytes"
single_cell_10K$merged_celltype_refined <- factor(cell_chr)



#--Extract subset: this operation helps me check the clustering of newly assigned cell subpopulations
subset_cells <- subset(single_cell_10K, idents = c("Monocytes"))
table(Idents(subset_cells))


subset_cells <- annotate_monocyte_subtypes(
  subset_cells, celltype_col="merged_celltype", monocyte_label="Monocytes",
  assay="SCT", delta=0.25, min_expr=0.1, out_col="merged_celltype_refined"
)
table(subset_cells$merged_celltype_refined)



#--Extract subset:
subset_cells <- subset(single_cell_10K, idents = c("T_cells", "Monocytes"))
table(Idents(subset_cells))


#1---Create---cell_to_6celltypes_10K.csv file, used for RSS and TRS calculation
dd<-as.data.frame(single_cell_10K$merged_celltype2)
dd$cell <- row.names(dd)
dd2 <- dd[,c(2,1)]
colnames(dd2) <- c("cell","cell_type")
write.csv(dd2, file="cell_to_6celltypes_10K.csv", quote=F, row.names=F)


#2---Create---cell_to_2celltypes_8.8K.csv file, used for RSS and TRS calculation
df<-as.data.frame(subset_cells$merged_celltype2)
df$cell <- row.names(df)
df2 <- dd[,c(2,1)]
colnames(df2) <- c("cell","cell_type")
write.csv(df2, file="cell_to_6celltypes_8.8K.csv", quote=F, row.names=F)




##Re-clustering analysis:
DefaultAssay(subset_cells)<-"RNA"
subset_cells <- NormalizeData(subset_cells, normalization.method = "LogNormalize", scale.factor = 10000)
subset_cells <- FindVariableFeatures(subset_cells, selection.method = "vst", nfeatures = 3000)
subset_cells <- ScaleData(subset_cells, verbose = FALSE)
subset_cells <- RunPCA(subset_cells, npcs = 20, verbose = FALSE)
subset_cells  <- FindNeighbors(subset_cells, dims = 1:20)
subset_cells  <- FindClusters(object = subset_cells, resolution = 0.2, graph.name="wsnn")

subset_cells <- RunUMAP(subset_cells, reduction = "pca", dims = 1:20,reduction.name = 'umap.rna', reduction.key = 'UMAP.RNA_')


pal <- c("#A13B46","#67ADB7")
DefaultAssay(subset_cells)<-"RNA"
pdf("subset_cells_cell_type2_8.8Kcells.pdf")
DimPlot(subset_cells,reduction = "umap.rna",group.by = "merged_celltype2",label = T,cols=pal, pt.size = 0.1, repel = T)
dev.off()

pal <- c("#A13B46","#67ADB7")
DefaultAssay(subset_cells)<-"RNA"
pdf("subset_cells_cell_type2_8.7Kcells_UMAP.pdf")
DimPlot(subset_cells,reduction = "umap",group.by = "merged_celltype2",label = T,cols=pal, pt.size = 0.1, repel = T)
dev.off()

pal <- c("#A13B46","#67ADB7","lightblue","orange")
DefaultAssay(subset_cells)<-"RNA"
pdf("subset_cells_cell_type_8.8Kcells.pdf")
DimPlot(subset_cells,reduction = "umap.rna",group.by = "merged_celltype",label = T,cols=pal, pt.size = 0.1, repel = T)
dev.off()

#saveRDS(subset_cells, file="10X_PBMC_subset_8.8K.rds")


###----Extract matrix for scRBP
subset_cells_matrix <- GetAssayData(subset_cells, assay = "RNA", layer = "counts")
 
# Save as CSV file
write.csv(as.matrix(subset_cells_matrix), file = "10X_PBMC_subset_8.8K.csv", quote=F, row.names = TRUE)


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
write_feather(df, "10X_PBMC_subset_8.8K.feather")
