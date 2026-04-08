
#2025-05-16

setwd("/Users/mayunlong/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/00-Figures/")

# Load necessary libraries
library(ggplot2)
 

# Load the data
data_file <- "Sup_Figure_S1C.csv"  # Replace with your actual file path
df <- read.csv(data_file)

 
# Set color
df$Color_Group <- ifelse(df$Information_content < 5 | df$Log_Odds_Score < 5, "gray", "#4DBBD5FF")
 
ggplot(df, aes(x = Information_content, y = Log_Odds_Score, color = Color_Group)) + 
  geom_point(alpha = 0.7, size = 2) + 
  geom_vline(xintercept = 5, linetype = "dashed", color = "blue", size = 0.8) + 
  geom_hline(yintercept = 5, linetype = "dashed", color = "blue", size = 0.8) + 
  geom_smooth(data = subset(df, Information_content >= 5 & Log_Odds_Score >= 5), 
              method = "lm", color = "darkred", linetype = "dashed", se = FALSE, size = 0.8) + 
  labs(title = "Information Content vs. Log_Odds_Score (Unclustered Motifs)", 
       x = "Information Content (bits)", 
       y = "Log_Odds_Score (bits)") + 
  scale_color_identity() +
  theme_classic()

# Calculate correlation and significance
cor_test <- cor.test(df$Information_content, df$Log_Odds_Score, method="pearson")
# Print correlation results
print(cor_test)


