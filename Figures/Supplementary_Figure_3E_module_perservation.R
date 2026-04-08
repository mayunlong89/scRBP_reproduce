library(tidyverse)

# Get file names
files_100K <- list.files(pattern = "_100K_grn_network_.*csv$")
files_200K <- list.files(pattern = "_200K_grn_network_.*csv$")
method_names <- gsub("_\\d+K_grn_network_.*", "", files_100K)

# Define Jaccard function
jaccard <- function(x, y) {
  inter <- length(intersect(x, y))
  union <- length(union(x, y))
  if (union == 0) return(NA) else return(inter / union)
}

# Store Jaccard scores for all methods and all RBPs
all_jaccard_scores <- data.frame()

for (method in method_names) {
  file_100K <- files_100K[grepl(method, files_100K)]
  file_200K <- files_200K[grepl(method, files_200K)]
  
  df_100K <- read_csv(file_100K, show_col_types = FALSE)
  df_200K <- read_csv(file_200K, show_col_types = FALSE)
  
  rbp_targets_100K <- df_100K %>% group_by(RBP) %>% summarise(targets = list(Gene), .groups = "drop") %>% deframe()
  rbp_targets_200K <- df_200K %>% group_by(RBP) %>% summarise(targets = list(Gene), .groups = "drop") %>% deframe()
  
  shared_rbps <- intersect(names(rbp_targets_100K), names(rbp_targets_200K))
  
  # Compute Jaccard for each RBP
  rbp_jaccard <- map_dbl(shared_rbps, ~ jaccard(rbp_targets_100K[[.x]], rbp_targets_200K[[.x]]))
  
  # Organize into a long data frame
  temp_df <- data.frame(
    Method = method,
    RBP = shared_rbps,
    Jaccard = rbp_jaccard
  )
  
  all_jaccard_scores <- bind_rows(all_jaccard_scores, temp_df)
}




# Visualization: preservation boxplot
pdf("module_construction_Persevation_15conditions_100K_200K.pdf", width = 8, height = 8)

ggplot(all_jaccard_scores, aes(x = reorder(Method, Jaccard, median), y = Jaccard)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue") +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +
  theme_minimal() +
  coord_flip() +
  labs(title = "RBP Module Preservation Across Methods (100K vs 200K)",
       y = "Jaccard Similarity per RBP", x = "Method")

# Close PDF output
dev.off()




# Step 1: Compute proportion of RBPs with Jaccard > 0.5 for each method
high_quality_ratio <- all_jaccard_scores %>%
  group_by(Method) %>%
  summarise(Proportion_GT_0.5 = mean(Jaccard > 0.5, na.rm = TRUE))

# Step 2: Add text annotation on the boxplot

library(ggplot2)

# To align y-axis coordinates (for annotation on boxplot), compute max value per group for label placement
label_positions <- all_jaccard_scores %>%
  group_by(Method) %>%
  summarise(MaxJ = max(Jaccard, na.rm = TRUE))

# Merge the two data frames
label_data <- left_join(high_quality_ratio, label_positions, by = "Method")

# Plot
ggplot(all_jaccard_scores, aes(x = reorder(Method, Jaccard, median), y = Jaccard)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue") +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.4) +
  geom_text(data = label_data,
            aes(x = Method, y = MaxJ + 0.05, 
                label = paste0(">", "0.5: ", round(Proportion_GT_0.5 * 100), "%")),
            size = 3.2, hjust = 0) +
  theme_minimal() +
  coord_flip() +
  ylim(0, 1.1) +
  labs(title = "RBP Module Preservation (Jaccard > 0.5 Proportion)",
       y = "Jaccard Similarity", x = "Method")




###-----Adjusted regulon size for module preservation

# Regress Jaccard ~ log(ModuleSize) separately for each method
library(dplyr)
library(broom)





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

# Compute regulon size (i.e., number of targets) for each RBP in each method
regulon_sizes <- lapply(names(networks), function(method) {
  rbp_targets <- networks[[method]]
  data.frame(
    Method = gsub("_200K_grn_network_IM0005_seed42\\.csv", "", method),
    RBP = names(rbp_targets),
    Size = sapply(rbp_targets, length)
  )
}) %>% bind_rows()





# 1. Merge the two data frames
merged_df <- all_jaccard_scores %>%
  inner_join(regulon_sizes, by = c("Method", "RBP"))

# 2. Regress Jaccard ~ log1p(Size) within each method group, and compute residuals
merged_df <- merged_df %>%
  group_by(Method) %>%
  mutate(Jaccard_resid = resid(lm(Jaccard ~ log1p(Size)))) %>%
  ungroup()



# Cluster information (same as before)
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
merged_df <- merged_df %>%
  left_join(method_cluster, by = "Method")

# Set plotting order by Cluster (reversed so top500perRBP appears at the top)
merged_df$Method <- factor(
  merged_df$Method,
  levels = rev(method_cluster %>% arrange(Cluster) %>% pull(Method))
)



pdf("module_construction_Persevation_15conditions_100K_200K_residual_jaccard_clustered.pdf", width = 8, height = 8)

ggplot(merged_df, aes(x = Method, y = Jaccard_resid)) +
  geom_boxplot(fill = "lightgreen", outlier.shape = NA) +
  coord_flip() +
  theme_classic() +
  labs(title = "Size-Corrected Module Preservation Across Methods",
       x = "Method (Cluster Ordered)", y = "Residual Jaccard (corrected for regulon size)") +
ylim(0,0.25)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

dev.off()





library(ggplot2)
# Visualization: preservation boxplot
pdf("module_construction_Persevation_15conditions_100K_200K_residual_jaccard.pdf", width = 8, height = 8)

ggplot(merged_df, aes(x = reorder(Method, Jaccard_resid, median), y = Jaccard_resid)) +
  geom_boxplot(fill = "lightgreen", outlier.shape = NA) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Size-Corrected Module Preservation Across Methods",
       x = "Method", y = "Residual Jaccard (corrected for regulon size)")+
 geom_hline(yintercept = 0, linetype = "dashed", color = "red")

# Close PDF output
dev.off()



avg_resid <- merged_df %>%
  group_by(Method) %>%
  summarise(Mean_Residual = mean(Jaccard_resid, na.rm = TRUE))

pdf("module_construction_Persevation_15conditions_100K_200K_residual_jaccard_barplot.pdf", width = 8, height = 8)

ggplot(avg_resid, aes(x = reorder(Method, Mean_Residual), y = Mean_Residual)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_classic() +
  labs(title = "Mean Residual Jaccard per Method",
       x = "Method", y = "Mean Residual (Jaccard ~ log(Size))")


# Close PDF output
dev.off()
