#2025-05-16

setwd("/Users/mayunlong/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/00-Figures/")


# Load necessary libraries
library(ggplot2)
library(reshape2)
library(ggrepel)


data <- read.csv("Sup_Figure_S1A.csv")

# convert data format
data_melted <- melt(data, id.vars = c("Years", "Database"), 
                    measure.vars = c("Clustered_motifs", "Unclustered_motifs"), 
                    variable.name = "Motif_Type", 
                    value.name = "Count")


data_melted <- melt(data, id.vars = c("Years", "Database"), 
                    measure.vars = c("Clustered_motifs", "Unclustered_motifs"), 
                    variable.name = "Motif_Type", 
                    value.name = "Count")

# set color
colors <- c("Clustered_motifs" = "#66c2a5", "Unclustered_motifs" = "#fc8d62")

# plot
ggplot(data_melted, aes(x = reorder(Database, Count), y = Count, fill = Motif_Type)) + 
  geom_bar(stat = "identity") + 
  #geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3, color = "black") +
  labs(title = "Motif Collection by Database", 
       x = "Motif Collection", 
       y = "Number of Motifs") + 
  theme_minimal(base_size = 14) + 
  theme(axis.text.y = element_text(size = 12, hjust = 1)) + 
  scale_fill_manual(values = colors, name = "Exclusiveness") + 
  coord_flip()


