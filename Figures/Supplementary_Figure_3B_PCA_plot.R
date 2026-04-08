# 1. Load required packages
library(ggplot2)
library(ggforce)  # for geom_mark_ellipse or geom_mark_hull
library(dplyr)
# Assumes pca_df already contains PCA results, including PC1, PC2, and method columns
# 2. Clustering: e.g., k-means (can also use cutree)
set.seed(42)
k <- 4  # Number of clusters, adjustable
km <- kmeans(pca_df[, c("PC1", "PC2")], centers = k)
pca_df$cluster <- as.factor(km$cluster)
pdf("module_construction_PCA_15conditions_with_ellipse.pdf", width = 8, height = 8)
# 3. Plot (with enclosing ellipses)
ggplot(pca_df, aes(x = PC1, y = PC2, label = method)) +
  geom_point(aes(color = cluster), size = 3) +
  geom_text(size = 2, vjust = -1.2, check_overlap = TRUE) +
  ggforce::geom_mark_ellipse(aes(fill = cluster, label = cluster), alpha = 0.15, show.legend = FALSE) +
  theme_minimal() +
  labs(title = "PCA Clustering of Method Similarity (Jaccard)",
       color = "Cluster", fill = "Cluster")
# Close PDF output
dev.off()


##-----Alternative

library(umap)
library(ggplot2)
# Ensure matrix is numeric type
jaccard_df <- as.data.frame(jaccard_matrix)
jaccard_df <- jaccard_df[order(rownames(jaccard_df)), order(colnames(jaccard_df))]
# --- PCA ---
pca <- prcomp(jaccard_df, scale. = TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$method <- rownames(jaccard_df)
pdf("module_construction_PCA_15conditions.pdf", width = 6, height = 6)
ggplot(pca_df, aes(x = PC1, y = PC2, label = method)) +
  geom_point(size = 3) +
  geom_text(size = 2, vjust = -0.8, check_overlap = TRUE) +
  theme_minimal() +
  labs(title = "PCA of Method Similarity (Jaccard)")
# Close PDF output
dev.off()

# --- UMAP ---
set.seed(42)
umap_result <- umap(jaccard_df)
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$method <- rownames(jaccard_df)
pdf("module_construction_UMAP_15conditions.pdf", width = 6, height = 6)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, label = method)) +
  geom_point(size = 3) +
  geom_text(size = 2, vjust = -0.8, check_overlap = TRUE) +
  theme_minimal() +
  labs(title = "UMAP of Method Similarity (Jaccard)")
# Close PDF output
dev.off()
