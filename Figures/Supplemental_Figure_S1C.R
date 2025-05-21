
#2025-05-16

setwd("/Users/mayunlong/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/00-Figures/")

# Load necessary libraries
library(ggplot2)
library(readr)

# Load the data
data_file <- "Sup_Figure_S1C.csv"  # Replace with your actual file path
df <- read.csv(data_file)

# Set the minimum value explicitly
min_info_content <- min(df$Information_content)
max_info_content <- max(df$Information_content)


##---------categorized by 13 resources
# Create the violin + boxplot
ggplot(df, aes(x = Database, y = Information_content, fill = Database)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", alpha = 0.8) +
  labs(title = "Motif Information Content Distribution by Database", 
       x = "Database", y = "Information Content") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  scale_fill_brewer(palette = "Set3") +
  ylim(min_info_content, max_info_content)  # Set y-axis limits
