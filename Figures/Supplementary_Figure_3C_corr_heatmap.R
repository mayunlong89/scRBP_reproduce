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

