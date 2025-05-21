# Load the universalmotif library
library(universalmotif)

# Path to your MEME format file
meme_file_path <- "04_merged_all_25044motifs_all_databases.meme"

# Read in the MEME file as a list of motifs
motifs <- read_meme(meme_file_path)

# Function to calculate Shannon entropy for a motif
calculate_entropy <- function(ppm) {
  # Calculate Shannon entropy for each position
  entropy <- 0
  for (i in 1:ncol(ppm)) {
    col_entropy <- -sum(ifelse(ppm[, i] > 0, ppm[, i] * log2(ppm[, i]), 0))
    entropy <- entropy + col_entropy
  }
  return(entropy)
}

# Initialize an empty data frame to store results
results <- data.frame(Motif_ID = character(), Shannon_Entropy = numeric(), stringsAsFactors = FALSE)

# Loop through each motif, calculate Shannon entropy, and store results
for (motif in motifs) {
  motif_id <- motif@name  # Extract the motif ID from the motif object
  ppm <- convert_type(motif, "PPM")@motif  # Convert to PPM (Position Probability Matrix)
  entropy <- calculate_entropy(ppm)  # Calculate Shannon entropy
  results <- rbind(results, data.frame(Motif_ID = motif_id, Shannon_Entropy = entropy))
}

# Save results to a file
write.csv(results, file = "04_merged_all_25044motifs_all_databases_shannon_entropy.csv", row.names = FALSE, quote=FALSE)

# Print message when done
print("Shannon entropy calculation complete. Results saved to motif_shannon_entropy.csv.")
