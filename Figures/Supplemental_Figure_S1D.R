

#2025-05-17


setwd("/Users/mayunlong/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/00-Figures/")

# Load necessary libraries
library(ggplot2)
library(readr)

# Load the data
data_file <- "Sup_Figure_S1D.csv"  # Replace with your actual file path
df <- read.csv(data_file)

head(df)

# Summary of Shannon entropy
summary(df$Shannon_entropy)
# Determine optimal threshold based on quantile
quantile(df$Shannon_entropy, probs = c(0.2, 0.5, 0.8))

# Plot the density of Shannon entropy
ggplot(df, aes(x=Shannon_entropy)) +
  geom_density(fill="lightblue", alpha=0.7) +
  geom_vline(xintercept=10, linetype="dashed", color="red", size=0.8) +
  #geom_vline(xintercept=8, linetype="dashed", color="blue", size=0.8) +
  labs(title="Shannon Entropy Distribution",
       x="Shannon Entropy (bits)",
       y="Density") +
  theme_minimal()


#df_ATtRACT <- df[which(df$Database=="ATtRACT"),]
df_ATtRACT <- df[which(df$Old_motifIDs!="new_ID"),]

# Plot the density of Shannon entropy
ggplot(df_ATtRACT, aes(x=Shannon_entropy)) +
  geom_density(fill="lightblue", alpha=0.7) +
  geom_vline(xintercept=10, linetype="dashed", color="red", size=0.8) +
  #geom_vline(xintercept=8, linetype="dashed", color="blue", size=0.8) +
  labs(title="Shannon Entropy Distribution",
       x="Shannon Entropy (bits)",
       y="Density") +
  theme_minimal()

length(df_ATtRACT$MotifIDs[which(df_ATtRACT$Shannon_entropy >10)])/length(df_ATtRACT$MotifIDs) 





#df_ATtRACT <- df[which(df$Database=="ATtRACT"),]
df_new <- df[which(df$Old_motifIDs=="new_ID" & df$Clustered_motifs == "No"),]

# Plot the density of Shannon entropy
ggplot(df_new, aes(x=Shannon_entropy)) +
  geom_density(fill="lightblue", alpha=0.7) +
  geom_vline(xintercept=10, linetype="dashed", color="red", size=0.8) +
  #geom_vline(xintercept=8, linetype="dashed", color="blue", size=0.8) +
  labs(title="Shannon Entropy Distribution",
       x="Shannon Entropy (bits)",
       y="Density") +
  theme_minimal()

length(df_new$MotifIDs[which(df_new$Shannon_entropy >10)])/length(df_new$MotifIDs) 




