
##Note---需要对tomtom.tsv 文件进行处理--

grep -v "^#" tomtom.tsv > tomtom_v2.tsv

tomtom_results <- read.csv("tomtom_v2.tsv", header = TRUE, sep = "\t")



## Remove motifs without similar motifs (q-value < 10-5)
##------------------------------------------------------------
import pandas as pd

# Load your data from a CSV file (adjust the file path and format as needed)
file_path = 'tomtom_v2.tsv'  # Replace with the actual file path

# If the file is in CSV format
df = pd.read_csv(file_path,sep='\t')

# If the file is in TSV format (you can comment out the CSV line and use this one if needed)
# df = pd.read_table(file_path)

# Step 1: Remove rows where Query_ID and Target_ID are the same
df_filtered = df[df["Query_ID"] != df["Target_ID"]]

# Step 2: Group by Query_ID and check if any q-value in the group is less than 10^-5
# If none of the q-values for a given Query_ID is less than 10^-5, remove all its rows
grouped = df_filtered.groupby("Query_ID").filter(lambda group: (group["q-value"] < 1e-5).any())

# Display or save the filtered DataFrame
grouped.to_csv('tomtom_v4_remove_Query_ID_10e-3_motifs.csv',index=False)

# To save the result to a CSV file, uncomment the following line:
# grouped.to_csv('filtered_motif_data.csv', index=False)

##------------------------------------------------------------




library(data.table)
library(ggplot2)
library(Seurat)

##-----backup packages----
library(FactoMineR)
library(factoextra)
library(umap)
library(Rtsne)
##------------------------


# (1)---------读取 TOMTOM 的输出

tomtom_results <- read.csv("tomtom_v3.csv", header = TRUE, sep = "\t")


# (2)---------获取所有的 motifs 列表
# 1. Convert the tomtom_results data to a data.table for faster operations
setDT(tomtom_results)
motifs <- unique(c(tomtom_results$Query_ID, tomtom_results$Target_ID))


# (3)---------初始化一个空的相似性矩阵，使用负无穷(-Inf)作为默认值

# 3. Initialize the similarity matrix with default values of -Inf
similarity_matrix <- matrix(-Inf, nrow = length(motifs), ncol = length(motifs), dimnames = list(motifs, motifs))


# 4. Use data.table's efficient processing to update the similarity matrix
# Only considering q-values <= 1 and using -log10(q-value) + 1e-45 as the similarity metric
valid_results <- tomtom_results[q.value <= 1]

# Convert q-values to similarity scores
valid_results[, similarity := -log10(q.value) + 1e-45]

# 5. Fill the similarity matrix efficiently using matrix indexing
# Faster way using matrix row/column indexing for batch updates
idx_query <- match(valid_results$Query_ID, motifs)
idx_target <- match(valid_results$Target_ID, motifs)

# Assign similarities to both directions in the matrix (query <-> target and target <-> query)
similarity_matrix[cbind(idx_query, idx_target)] <- valid_results$similarity
similarity_matrix[cbind(idx_target, idx_query)] <- valid_results$similarity

# Replace remaining -Inf values with the minimum similarity
min_similarity <- min(similarity_matrix[similarity_matrix != -Inf])
similarity_matrix[similarity_matrix == -Inf] <- min_similarity

# Ensure no NA or infinite values
similarity_matrix[is.infinite(similarity_matrix)] <- 0
similarity_matrix[is.na(similarity_matrix)] <- 0


##这一步结束，可以把similarity_matrix数据存起来

saveRDS(similarity_matrix, file="similarity_matrix.rds")
saveRDS(similarity_matrix, file="similarity_matrix_10e-3_motifs.rds")

#--Part II------------------------------------------------------------------
#1) 加载Seurat 库 并创建Seurat 对象
# 将相似性矩阵转换为 Seurat 对象
seurat_obj <- CreateSeuratObject(counts = similarity_matrix)

# 检查 Seurat 对象的基本信息
#print(seurat_obj)

#2）标准化数据 (Scaling)
#Seurat 的标准化通常通过 ScaleData 函数进行。这里我们将对相似性矩阵进行标准化处理。

# 对数据进行标准化处理
#Normalize the data before scaling:
seurat_obj <- NormalizeData(seurat_obj)

#Scaling
seurat_obj <- ScaleData(seurat_obj)

# 查看标准化后的数据
#head(seurat_obj@assays$RNA$scale.data)


#3）PCA降维分析
#Seurat 使用 RunPCA 来进行主成分分析 (PCA)。

# 运行 PCA
seurat_obj <- RunPCA(seurat_obj, features = rownames(similarity_matrix),npcs=100,seed.use=42, verbose=FALSE)

# 查看前几个主成分的载荷
#print(seurat_obj[["pca"]], dims = 1:2, nfeatures = 5)

##查看PCA---a ranking of principle components based on the percentage of variance explained by each one.
pdf("seurat_ElbowPlot.pdf", width = 10, height = 8)
ElbowPlot(seurat_obj)
dev.off()


###remove duplicates
pca_embeddings <- Embeddings(seurat_obj, "pca")
duplicated_cells <- duplicated(pca_embeddings)

#sum(duplicated_cells)
# Remove duplicated cells
seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[!duplicated_cells])


###neighbors and clusters-------------------------------
#10138 motifs： 274 Singletons identified and 151 final clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:100)
seurat_obj <- FindClusters(seurat_obj, resolution = 15)


#4）UMAP 降维分析
# 运行 UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:100)  # 根据需要选择PCA维度
 
#5）TSNE降维分析
# 运行 t-SNE
seurat_obj <- RunTSNE(seurat_obj, dims = 1:100) # 根据需要选择PCA维度


table(seurat_obj$seurat_clusters)

df <- as.data.frame(seurat_obj$seurat_clusters)

write.csv(df,file="motifs_with_clusters_remove_motifs_without_similarity.csv",quote = F)

#saving data

saveRDS(seurat_obj, file="seurat_obj_10138motifs_remove_motifs_without_similarity.rds")

saveRDS(seurat_obj, file="seurat_obj_10138motifs_remove_10e-3_all_unsimilarity_motifs.rds")


#6）保存PCA和UMAP结果
pdf("seurat_PCA_plot.pdf", width = 10, height = 8)
DimPlot(seurat_obj, reduction = "pca", group.by = "orig.ident") +
  ggtitle("PCA of Motif Similarity Matrix")
dev.off()

pdf("seurat_UMAP_plot.pdf", width = 10, height = 8)
DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident") +
  ggtitle("UMAP of Motif Similarity Matrix")
dev.off()

pdf("seurat_UMAP_plot_103clusters_10e-3.pdf", width = 10, height = 8)
DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("UMAP of Motif Similarity Matrix")
dev.off()


pdf("seurat_TSNE_plot_103clusters_10e-3.pdf", width = 10, height = 8)
DimPlot(seurat_obj, reduction = "tsne", group.by = "seurat_clusters") +
  ggtitle("UMAP of Motif Similarity Matrix")
dev.off()


# 可视化 t-SNE 结果
DimPlot(seurat_obj, reduction = "tsne", group.by = "orig.ident") + 
  ggtitle("t-SNE of Motif Similarity Matrix")


#7）提取降维坐标
#将降维后的坐标提取出来供进一步分析或作图

# 提取 PCA 坐标
pca_coords <- Embeddings(seurat_obj, reduction = "pca")

# 提取 UMAP 坐标
umap_coords <- Embeddings(seurat_obj, reduction = "umap")

# 查看前几行的降维坐标
head(pca_coords)
head(umap_coords)







setDT(tomtom_results2)
motifs <- unique(c(tomtom_results2$Query_ID, tomtom_results2$Target_ID))


# (3)---------初始化一个空的相似性矩阵，使用负无穷(-Inf)作为默认值

# 3. Initialize the similarity matrix with default values of -Inf
similarity_matrix <- matrix(-Inf, nrow = length(motifs), ncol = length(motifs), dimnames = list(motifs, motifs))


# 4. Use data.table's efficient processing to update the similarity matrix
# Only considering q-values <= 1 and using -log10(q-value) + 1e-45 as the similarity metric
valid_results <- tomtom_results2[q.value <= 1]

# Convert q-values to similarity scores
valid_results[, similarity := -log10(q.value) + 1e-45]

# 5. Fill the similarity matrix efficiently using matrix indexing
# Faster way using matrix row/column indexing for batch updates
idx_query <- match(valid_results$Query_ID, motifs)
idx_target <- match(valid_results$Target_ID, motifs)

# Assign similarities to both directions in the matrix (query <-> target and target <-> query)
similarity_matrix[cbind(idx_query, idx_target)] <- valid_results$similarity
similarity_matrix[cbind(idx_target, idx_query)] <- valid_results$similarity

# Replace remaining -Inf values with the minimum similarity
min_similarity <- min(similarity_matrix[similarity_matrix != -Inf])
similarity_matrix[similarity_matrix == -Inf] <- min_similarity

# Ensure no NA or infinite values
similarity_matrix[is.infinite(similarity_matrix)] <- 0
similarity_matrix[is.na(similarity_matrix)] <- 0


