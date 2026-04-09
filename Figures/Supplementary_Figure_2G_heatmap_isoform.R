
##----2025-06-20
#2025-06-06
getwd()
setwd("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/00-motif_target_isoform_ranking_databases")
# Read data
df <- read.csv("./Supplementary_Figure_S2G/spearman_correlation_heatmap.csv")
# Remove the first column (RBP names)
rank_matrix <- df[ , -1]
# Convert all columns to numeric (skipping the RBP name column)
rank_numeric <- as.data.frame(lapply(rank_matrix, function(x) as.numeric(as.character(x))))
# Compute Spearman correlation (NAs handled automatically)
cor_matrix<-cor(rank_numeric, method = "spearman", use = "pairwise.complete.obs")
# Visualize
library(pheatmap)
pheatmap(cor_matrix, main = "Spearman Correlation Between Ranking Methods")
# Generate color scheme
my_colors <- colorRampPalette(c("white","#F6DFC0","#FEA087","#B2182B"))(100)
# Plot
pheatmap(cor_matrix,
         color = my_colors,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = TRUE,        # display values in cells
         number_format = "%.2f",        # keep two decimal places
         main = "Spearman Correlation Between Ranking Methods")
