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
