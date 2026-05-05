#!/bin/bash

# Check the input parameter
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

# Input and output directory
input_dir="$1"
output_dir="$2"

# Tool path
snakemake_loc="/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/01_snakemake_pipeline/"

# Ensure the output directory exists
mkdir -p "$output_dir"

# Loop through all the .ihbcp files in the input directory
for ihbcp_file in "$input_dir"/*.ihbcp; do

    # Extract filename without extension
    base_name=$(basename "$ihbcp_file" .ihbcp)

    # Define the name of the output '.meme' file
    meme_file="${output_dir}/${base_name}.meme"

    # Define log file path

    log_file="${output_dir}/${base_name}.log"



    # Run the script
    python "${snakemake_loc}/bamm2pwm.py" "$ihbcp_file" "$meme_file" > "$log_file" 2>&1

done

