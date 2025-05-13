
# Load necessary libraries
library(ggplot2)

# Create the data frame
data <- data.frame(
  Years = c(2006, 2011, 2013, 2013, 2014, 2015, 2016, 2017, 2020, 2022, 2023, 2023, 2023),
  Database = c("AEDB", "RBPDB", "SpliceAid-F", "CISBP-RNA", "RBPmap", "DoRiNA3", "ATtRACT", 
               "MotifMap-RNA", "oRNAmEnt", "POSTAR3", "ENCORI", "RBP-Tar", "GEO-2024"),
  Number = c(108, 91, 18, 153, 223, 638, 1195, 244, 340, 4374, 9937, 3497, 4744)
)

# Generate the bar plot
ggplot(data, aes(x = reorder(Database, -Number), y = Number, fill = as.factor(Years))) +
  geom_bar(stat = "identity") +
  labs(title = "Number of RBPs in Various Databases", x = "Database", y = "Number") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete(name = "Year") +
  coord_flip()
