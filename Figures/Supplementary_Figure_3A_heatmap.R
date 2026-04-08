#2025-11-11
library(pheatmap)

# Read matrix (same as before)
epr <- read.csv("01_GRN_assessment_EPR_summary.csv", header = TRUE, sep = ",", check.names = FALSE)
rownames(epr) <- epr[[1]]
epr <- as.matrix(epr[ , -1, drop=FALSE])

# Gradient (values > 1 become brighter)
pal_hi <- colorRampPalette(c("#8C2981", "#DE4968", "#FE9F6D", "#FDE725"))

# Segmented breaks: one segment [min, 1] + N segments [1, max]
n_hi   <- 100
maxlim <- max(1.4, max(epr, na.rm = TRUE))         # Adjust upper limit as needed
minlim <- min(epr, na.rm = TRUE)
breaks <- c(minlim, seq(1, maxlim, length.out = n_hi + 1))  # Total n_hi+1 segments (including one [min,1] segment)
cols   <- c("black", pal_hi(n_hi))                          # 1 black segment + n_hi gradient segments

pheatmap(epr,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = cols, breaks = breaks,
         display_numbers = TRUE, number_format = "%.2f",
         fontsize_number = 9, border_color = NA,
         main = "EPR (higher = better)"
)
