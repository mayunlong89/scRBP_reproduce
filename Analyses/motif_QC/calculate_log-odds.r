
# Load the universalmotif library
library(universalmotif)

# Path to your MEME format file
meme_file_path <- "04_merged_all_25044motifs_all_databases.meme"

# Read in the MEME file as a list of motifs
motifs <- read_meme(meme_file_path)

# Set background probability (uniform distribution)
background_prob <- 0.25

# Function to calculate log-odds scores for a motif
calculate_log_odds <- function(ppm) {
  # Calculate log-odds matrix (PWM)
  pwm_matrix <- log2(ppm / background_prob)
  
  # Calculate total log-odds score as the maximum possible score for the motif
  log_odds_score <- sum(apply(pwm_matrix, 2, max))
  
  return(log_odds_score)
}

# Initialize an empty data frame to store results
results <- data.frame(Motif_ID = character(), Log_Odds_Score = numeric(), stringsAsFactors = FALSE)

# Loop through each motif, calculate log-odds score, and store results
for (motif in motifs) {
  motif_id <- motif@name  # Extract the motif ID from the motif object
  ppm <- convert_type(motif, "PPM")@motif  # Convert to PPM (Position Probability Matrix)
  log_odds <- calculate_log_odds(ppm)  # Calculate log-odds score
  results <- rbind(results, data.frame(Motif_ID = motif_id, Log_Odds_Score = log_odds))
}

# Save results to a file
write.csv(results, file = "04_merged_all_25044motifs_all_databases_log_odds_scores.csv", row.names = FALSE, quote=FALSE)

# Print message when done
print("Log-odds score calculation complete. Results saved to motif_log_odds_scores.csv.")


