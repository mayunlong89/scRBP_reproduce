#!/bin/bash
#SBATCH --job-name=cbust_3utr      # Set the job name
#SBATCH --output=cbust_3utr.log    # Set the log file for standard output
#SBATCH --error=cbust_3utr.log     # Set the log file for error messages
#SBATCH --mem=50G                  # Set the maximum memory allocation
#SBATCH --time=30-00:00:00         # Set the maximum runtime (30 days)
#SBATCH --cpus-per-task=10         # Set the number of CPUs per task

# Load the conda environment
source ~/.bashrc   # This step is required if .bashrc contains conda initialization
conda activate pyscenic

# Set file paths
fasta_ref=/mnt/isilon/gandal_lab/mayl/reference/transcript_regions
motifs=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/ATtRACT/ATtRACT/01_cluster_buster

# Run Cluster-Buster
cbust -g 35 -l -f 5 $motifs/combined_cbust_motifs.txt $fasta_ref/3UTR_transcript_regions.fa > output_3utr.bed
