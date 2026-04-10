# Inspect regulon sizes from a GMT file

# Select one GMT file to analyze (uncomment the desired region)
# gmt_file <- "trait_15Kcells_CDS_REAL_PLUS_NULLS.entrez.gmt"
# gmt_file <- "trait_15Kcells_Introns_REAL_PLUS_NULLS.entrez.gmt"
# gmt_file <- "trait_15Kcells_5UTR_REAL_PLUS_NULLS.entrez.gmt"
gmt_file <- "trait_15Kcells_3UTR_REAL_PLUS_NULLS.entrez.gmt"

# Read all lines from the GMT file
lines <- readLines(gmt_file)

# Split each line by tab
splitted <- strsplit(lines, "\t", fixed = TRUE)

# Extract regulon names (column 1)
regulon_id <- vapply(splitted, `[`, character(1), 1L)

# Count genes per regulon (total columns minus 2 for ID and description)
n_genes <- vapply(splitted, function(x) max(0L, length(x) - 2L), integer(1))

# Combine into a summary data frame
regulon_sizes <- data.frame(
  regulon = regulon_id,
  n_genes = n_genes,
  stringsAsFactors = FALSE
)

# Check overall distribution
head(regulon_sizes)
summary(regulon_sizes$n_genes)
quantile(regulon_sizes$n_genes, probs = c(0, .25, .5, .75, .9, .95, .99, 1))
