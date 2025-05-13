library(data.table)
library(ggplot2)
library(Seurat)

##-----Backup packages----
library(FactoMineR)
library(factoextra)
library(umap)
library(Rtsne)
##------------------------


# (1)---------Read TOMTOM output

# Load the TOMTOM results as a data.table for faster operations
tomtom_results <- read.csv("tomtom_v3.csv", header = TRUE, sep = "\t")


# (2)---------Get the list of all motifs
# Convert the data to a data.table for efficient processing
setDT(tomtom_results)
motifs <- unique(c(tomtom_results$Query_ID, tomtom_results$Target_ID))


# (3)---------Initialize an empty similarity matrix with -Inf as the default value

# Create a similarity matrix with -Inf as the default value
similarity_matrix <- matrix(-Inf, nrow = length(motifs), ncol = length(motifs), dimnames = list(motifs, motifs))


# (4)---------Fill the similarity matrix with valid q-values

# Filter the results to include only q-values less than or equal to 1
valid_results <- tomtom_results[q.value <= 1]

# Convert q-values to similarity scores using -log10(q-value) + 1e-45 as the metric
valid_results[, similarity := -log10(q.value) + 1e-45]

# (5)---------Efficiently populate the similarity matrix

# Use matrix indexing for efficient batch updates
idx_query <- match(valid_results$Query_ID, motifs)
idx_target <- match(valid_results$Target_ID, motifs)

# Assign similarity scores for both directions (query <-> target and target <-> query)
similarity_matrix[cbind(idx_query, idx_target)] <- valid_results$similarity
similarity_matrix[cbind(idx_target, idx_query)] <- valid_results$similarity

# Replace remaining -Inf values with the minimum similarity
min_similarity <- min(similarity_matrix[similarity_matrix != -Inf])
similarity_matrix[similarity_matrix == -Inf] <- min_similarity

# Ensure no NA or infinite values remain
similarity_matrix[is.infinite(similarity_matrix)] <- 0
similarity_matrix[is.na(similarity_matrix)] <- 0


## Save the similarity matrix for later use

saveRDS(similarity_matrix, file="similarity_matrix.rds")
saveRDS(similarity_matrix, file="similarity_matrix_10e-3_motifs.rds")


#--Part II------------------------------------------------------------------
# (6)---------Create Seurat object

# Create a Seurat object from the similarity matrix
seurat_obj <- CreateSeuratObject(counts = similarity_matrix)

# (7)---------Data normalization and scaling

# Normalize the data before scaling
seurat_obj <- NormalizeData(seurat_obj)

# Scale the data
seurat_obj <- ScaleData(seurat_obj)


# (8)---------PCA for dimensionality reduction

# Run PCA with 100 components
seurat_obj <- RunPCA(seurat_obj, features = rownames(similarity_matrix), npcs=100, seed.use=42, verbose=FALSE)

# Save the Elbow plot to visualize the variance explained by each principal component
pdf("seurat_ElbowPlot.pdf", width = 10, height = 8)
ElbowPlot(seurat_obj)
dev.off()


# (9)---------Remove duplicate embeddings

# Extract PCA embeddings
pca_embeddings <- Embeddings(seurat_obj, "pca")

# Identify duplicate cells
duplicated_cells <- duplicated(pca_embeddings)

# Remove duplicated cells
seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[!duplicated_cells])


# (10)---------Find neighbors and clusters

# Find neighbors and clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:100)
seurat_obj <- FindClusters(seurat_obj, resolution = 15)


# (11)---------Run UMAP for dimensionality reduction

# Run UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:100)  # Adjust the number of PCA dimensions as needed


# (12)---------Run t-SNE for dimensionality reduction

# Run t-SNE
seurat_obj <- RunTSNE(seurat_obj, dims = 1:100) # Adjust the number of PCA dimensions as needed


# (13)---------Save clustering results

# Save the cluster information to a CSV file
df <- as.data.frame(seurat_obj$seurat_clusters)
write.csv(df, file="motifs_with_clusters_remove_motifs_without_similarity.csv", quote = FALSE)

# Save the Seurat object for future use
saveRDS(seurat_obj, file="seurat_obj_10138motifs_remove_motifs_without_similarity.rds")
saveRDS(seurat_obj, file="seurat_obj_10138motifs_remove_10e-3_all_unsimilarity_motifs.rds")


# (14)---------Save PCA and UMAP plots

# Save the PCA plot
pdf("seurat_PCA_plot.pdf", width = 10, height = 8)
DimPlot(seurat_obj, reduction = "pca", group.by = "orig.ident") +
  ggtitle("PCA of Motif Similarity Matrix")
dev.off()

# Save the UMAP plot
pdf("seurat_UMAP_plot.pdf", width = 10, height = 8)
DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident") +
  ggtitle("UMAP of Motif Similarity Matrix")
dev.off()

# Save UMAP plot with cluster information
pdf("seurat_UMAP_plot_103clusters_10e-3.pdf", width = 10, height = 8)
DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("UMAP of Motif Similarity Matrix")
dev.off()

# Save t-SNE plot with cluster information
pdf("seurat_TSNE_plot_103clusters_10e-3.pdf", width = 10, height = 8)
DimPlot(seurat_obj, reduction = "tsne", group.by = "seurat_clusters") +
  ggtitle("t-SNE of Motif Similarity Matrix")
dev.off()


# (15)---------Extract PCA and UMAP coordinates

# Extract PCA coordinates
pca_coords <- Embeddings(seurat_obj, reduction = "pca")

# Extract UMAP coordinates
umap_coords <- Embeddings(seurat_obj, reduction = "umap")

# Display the first few rows of PCA and UMAP coordinates
head(pca_coords)
head(umap_coords)
