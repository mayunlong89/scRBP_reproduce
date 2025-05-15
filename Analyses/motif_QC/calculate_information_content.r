#2024-10-29 
#Calculate the information content of each motif within 10138 refined motifs


# Load the universalmotif library
library(universalmotif)

# Path to your MEME format file
meme_file_path <- "01_merged_3_10138motifs_all_databases_refined.meme"

# Read in the MEME file as a list of motifs
motifs <- read_meme(meme_file_path)

# Function to calculate information content for a motif
calculate_ic <- function(ppm) {
  ic <- 0
  for (i in 1:ncol(ppm)) {
    col_ic <- 2 + sum(ifelse(ppm[, i] > 0, ppm[, i] * log2(ppm[, i]), 0))
    ic <- ic + col_ic
  }
  return(ic)
}

# Initialize an empty data frame to store results
results <- data.frame(Motif_ID = character(), Information_Content = numeric(), stringsAsFactors = FALSE)

# Loop through each motif, calculate IC, and store results
for (motif in motifs) {
  motif_id <- motif@name  # Extract the motif ID from the motif object
  ppm <- convert_type(motif, "PPM")@motif  # Convert to PPM (Position Probability Matrix)
  ic <- calculate_ic(ppm)  # Calculate information content manually
  results <- rbind(results, data.frame(Motif_ID = motif_id, Information_Content = ic))
}

# Save results to a file
write.csv(results, file = "motif_information_content.csv", row.names = FALSE,quote=F)

# Print message when done
print("Information content calculation complete. Results saved to motif_information_content.csv.")



