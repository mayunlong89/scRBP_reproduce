# Load the universalmotif library
library(universalmotif)

# Function to calculate information content for a motif
calculate_ic <- function(ppm) {
  ic <- 0
  for (i in 1:ncol(ppm)) {
    col_ic <- 2 + sum(ifelse(ppm[, i] > 0, ppm[, i] * log2(ppm[, i]), 0))
    ic <- ic + col_ic
  }
  return(ic)
}

# Check command line arguments for input and output file paths
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Please provide exactly two arguments: input MEME file path and output CSV file path.")
}

input_file_path <- args[1]  # Path to the input MEME format file
output_file_path <- args[2]  # Path to the output CSV file

# Read in the MEME file as a list of motifs
motifs <- read_meme(input_file_path)

# Initialize an empty data frame to store results
results <- data.frame(Motif_ID = character(), Information_Content = numeric(), stringsAsFactors = FALSE)

# Loop through each motif, calculate IC, and store results
for (motif in motifs) {
  motif_id <- motif@name  # Extract the motif ID from the motif object
  ppm <- convert_type(motif, "PPM")@motif  # Convert to PPM (Position Probability Matrix)
  ic <- calculate_ic(ppm)  # Calculate information content manually
  results <- rbind(results, data.frame(Motif_ID = motif_id, Information_Content = ic))
}

# Save results to the specified output file
write.csv(results, file = output_file_path, row.names = FALSE, quote = FALSE)

# Print message when done
cat("Information content calculation complete. Results saved to", output_file_path, "\n")

