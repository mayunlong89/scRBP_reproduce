library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(readr)
library(pheatmap)
library(tidyr)
library(tibble)


# 1. Get all file names
files <- list.files(pattern = "*200K_grn_network_.*csv")

# 2. Read all files and build a named list
networks <- lapply(files, function(file) {
  df <- read_csv(file, show_col_types = FALSE)
  df %>%
    group_by(RBP) %>%
    summarise(targets = list(unique(Gene)), .groups = "drop") %>%
    deframe()
})
names(networks) <- files

# 3. Define function to compute Jaccard similarity for each RBP
jaccard <- function(set1, set2) {
  inter <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  if (union == 0) return(NA) else return(inter / union)
}

# 4. Compute average Jaccard for each pair of methods (per RBP)
method_pairs <- combn(files, 2, simplify = FALSE)
jaccard_matrix <- matrix(NA, nrow = length(files), ncol = length(files),
                         dimnames = list(files, files))

for (pair in method_pairs) {
  f1 <- pair[1]; f2 <- pair[2]
  rbps <- intersect(names(networks[[f1]]), names(networks[[f2]]))
  scores <- sapply(rbps, function(rbp) {
    jaccard(networks[[f1]][[rbp]], networks[[f2]][[rbp]])
  }, USE.NAMES = FALSE)
  avg_jaccard <- mean(scores, na.rm = TRUE)
  jaccard_matrix[f1, f2] <- jaccard_matrix[f2, f1] <- avg_jaccard
}

diag(jaccard_matrix) <- 1  # Self-comparison equals 1


# Extract original column names
old_names <- colnames(jaccard_matrix)

# Remove suffix
short_names <- gsub("_200K_grn_network_IM0005_seed42\\.csv", "", old_names)

# Replace row and column names
rownames(jaccard_matrix) <- short_names
colnames(jaccard_matrix) <- short_names



# 5. Visualize heatmap
# Create PDF file
pdf("module_construction_heatmap_15conditions.pdf", width = 10, height = 10)

# Output plot
pheatmap(jaccard_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         main = "Mean RBP-wise Jaccard Similarity Across Methods",
         angle_col = 45)

# Close PDF output
dev.off()




# ========== Additional section ========== Regulon size
# Compute regulon size (i.e., number of targets) for each RBP in each method
regulon_sizes <- lapply(names(networks), function(method) {
  rbp_targets <- networks[[method]]
  data.frame(
    Method = gsub("_200K_grn_network_IM0005_seed42\\.csv", "", method),
    RBP = names(rbp_targets),
    Size = sapply(rbp_targets, length)
  )
}) %>% bind_rows()


# Visualization: boxplot of regulon sizes
pdf("regulon_size_distribution_boxplot_15conditions.pdf", width = 10, height = 8)

ggplot(regulon_sizes, aes(x = reorder(Method, Size, median), y = Size)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue") +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.4) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Regulon Size Distribution Across Methods",
       x = "Method", y = "Number of Target Genes per RBP")

dev.off()


pdf("regulon_size_distribution_boxplot_15conditions_nopoints.pdf", width = 10, height = 8)

ggplot(regulon_sizes, aes(x = reorder(Method, Size, median), y = Size)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Regulon Size Distribution Across Methods",
       x = "Method", y = "Number of Target Genes per RBP")+
	ylim(1,3000)

dev.off()


### Classify and sort by cluster
### Classify and sort by cluster
### Classify and sort by cluster


method_cluster <- tibble::tibble(
  Method = c(
    "top500perRBP", "percent75", "top25perTarget", "top1000perRBP", "percent50", "top50perTarget",  # Cluster 1
    "top300perRBP", "top200perRBP", "percent90",                                                     # Cluster 2
    "top10perTarget", "top5perTarget", "percent95",                                                  # Cluster 3
    "top100perRBP", "top3perTarget",                                                                 # Cluster 4
    "top50perRBP"                                                                                     # Cluster 5
  ),
  Cluster = c(
    rep("Cluster 1", 6),
    rep("Cluster 2", 3),
    rep("Cluster 3", 3),
    rep("Cluster 4", 2),
    "Cluster 5"
  )
)


# Merge cluster information
regulon_sizes_clustered <- regulon_sizes %>%
  left_join(method_cluster, by = "Method")

# Set Method as factor, ordered by cluster and within-cluster method order
# Reverse order so top500perRBP appears at the top
regulon_sizes_clustered$Method <- factor(
  regulon_sizes_clustered$Method,
  levels = rev(method_cluster %>% arrange(Cluster) %>% pull(Method))
)


# Plot
pdf("regulon_size_distribution_boxplot_clustered.pdf", width = 10, height = 8)

ggplot(regulon_sizes_clustered, aes(x = Method, y = Size)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue") +
  coord_flip() +
  theme_classic() +
  labs(title = "Regulon Size Distribution (Clustered Methods)",
       x = "Method (Cluster Ordered)", y = "Number of Target Genes per RBP") +
  ylim(1, 3000)

dev.off()

pdf("regulon_size_distribution_boxplot_clustered_facet.pdf", width = 10, height = 8)

ggplot(regulon_sizes_clustered, aes(x = Method, y = Size)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue") +
  coord_flip() +
  facet_grid(Cluster ~ ., scales = "free_y", space = "free_y") +
  theme_classic () +
  labs(title = "Regulon Size Distribution by Cluster",
       x = "Method", y = "Number of Target Genes per RBP") +
  ylim(1, 3000)

dev.off()
