
#@ Step 4------trimming PWMs

##----------------------------------------------------------------------------------------------
#location:

universal_loc=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/01_snakemake_pipeline/


Rscript $universal_loc/universalmotif.r streme.txt trimmed_meme_file.meme




##########---------------------detailed codes-----------------------------------------------------------

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("universalmotif")

#load package
library(universalmotif)

#read STREME generated PWM;

streme_pwm <- read_meme("streme.txt")

#removing any peripheral positions with IC (<= 0.1) and rounding down base probabilities <= 0.025;

# define trim and rounding function 
trim_and_round <- function(motif, min_ic = 0.1, min_prob = 0.025) {
  # trim PWM motif，remove IC <= min_ic 
  trimmed_motif <- trim_motifs(motif, min.ic = min_ic)
  
  # rounding down，set the value of base probalities <= min_prob to 0
  pwm <- trimmed_motif@motif
  pwm[pwm <= min_prob] <- 0
  trimmed_motif@motif <- pwm
  
  return(trimmed_motif)
}

# trimmed motifs and rounding down base probabilities 
trimmed_motifs <- lapply(streme_pwm, trim_and_round, min_ic = 0.1, min_prob = 0.025)

#save trimmed PWM files

write_meme(trimmed_motifs, "trimmed_meme_file.meme")
