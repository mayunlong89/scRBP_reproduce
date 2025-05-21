
#2025-05-13

setwd("/Users/mayunlong/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/00-Figures/")


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(ggridges)

##----All motifs-----(1)
# Load the data
data_file <- "Sup_Figure_S1G.csv"   
df <- read.csv(data_file)

# Count the number of motifs per RBP
motif_counts2 <- table(df$RBPs)

# Convert table into data.frame
motif_df <- data.frame(RBP = names(motif_counts2), Motif_Count = as.numeric(motif_counts2))

# Prepare the data for fitting
count_df <- as.data.frame(table(motif_df$Motif_Count))
colnames(count_df) <- c("Motif_Count", "Frequency")
count_df$Motif_Count <- as.numeric(as.character(count_df$Motif_Count))

# Create the histogram with fitted line
library(ggplot2)
ggplot(count_df, aes(x = Motif_Count, y = Frequency)) +
  geom_bar(stat = "identity", fill = "#56B4E9", color = "black", alpha = 0.8, width = 0.9) +
  geom_smooth(method = "loess", se = FALSE, color = "darkred", size = 1.2, span = 0.3) +
  labs(title = "Motif Count Distribution per RBP",
       x = "Motif Count", 
       y = "Number of RBPs") +
  theme_minimal(base_size = 14)


###----clustered motifs-----(2)

df_clusteredmotifs <- df[which(df$Clustered_motifs=="Yes"),]
length(unique(df_clusteredmotifs$MotifIDs))
length(unique(df_clusteredmotifs$MotifIDs))/length(unique(df_clusteredmotifs$RBPs))



# Count the number of motifs per RBP
motif_counts <- table(df_clusteredmotifs$RBPs)
motif_df <- data.frame(RBP = names(motif_counts), Motif_Count = as.numeric(motif_counts))
# Prepare the data for fitting
count_df <- as.data.frame(table(motif_df$Motif_Count))
colnames(count_df) <- c("Motif_Count", "Frequency")
count_df$Motif_Count <- as.numeric(as.character(count_df$Motif_Count))

# Create the histogram with fitted line
ggplot(count_df, aes(x = Motif_Count, y = Frequency)) +
  geom_bar(stat = "identity", fill = "#56B4E9", color = "black", alpha = 0.8, width = 0.9) +
  geom_smooth(method = "loess", se = FALSE, color = "darkred", size = 1.2, span = 0.3) +
  labs(title = "Motif Count Distribution per RBP",
       x = "Motif Count", 
       y = "Number of RBPs") +
  theme_minimal(base_size = 14)




###----unclustered motifs-----(3)

df_unclusteredmotifs <- df[which(df$Clustered_motifs=="No"),]
length(unique(df_unclusteredmotifs$MotifIDs))
length(unique(df_unclusteredmotifs$MotifIDs))/length(unique(df_unclusteredmotifs$RBPs))



# Count the number of motifs per RBP
motif_counts <- table(df_unclusteredmotifs$RBPs)
motif_df <- data.frame(RBP = names(motif_counts), Motif_Count = as.numeric(motif_counts))
# Prepare the data for fitting
count_df <- as.data.frame(table(motif_df$Motif_Count))
colnames(count_df) <- c("Motif_Count", "Frequency")
count_df$Motif_Count <- as.numeric(as.character(count_df$Motif_Count))

# Create the histogram with fitted line
ggplot(count_df, aes(x = Motif_Count, y = Frequency)) +
  geom_bar(stat = "identity", fill = "#56B4E9", color = "black", alpha = 0.8, width = 0.9) +
  geom_smooth(method = "loess", se = FALSE, color = "darkred", size = 1.2, span = 0.3) +
  labs(title = "Motif Count Distribution per RBP",
       x = "Motif Count", 
       y = "Number of RBPs") +
  theme_minimal(base_size = 14)





