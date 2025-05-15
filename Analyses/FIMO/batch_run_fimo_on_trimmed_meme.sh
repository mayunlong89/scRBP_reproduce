#!/bin/bash

# Check if the required number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_dir> <output_dir_base> <ref_dir> <bgfile_dir>"
    exit 1
fi

# Get the input directory, output directory base path, reference sequence path, and background file path from command line arguments
input_dir="$1"
output_dir_base="$2"
ref_dir="$3"
bgfile_dir="$4"

# Set the reference sequence file and background file paths
ref_fasta="$ref_dir/gene_sequence_v44_GCRh38.fasta"
bg_file="$bgfile_dir/background.meme"

# Loop through each trimmed MEME file in the input directory
for meme_file in "$input_dir"/*_trimmed.meme; do
    # Extract the filename without path and extension
    filename=$(basename -- "$meme_file")
    name="${filename%.*}"

    # Set the output directory path
    output_dir="$output_dir_base/fimo_output_${name}"

    # Use srun to submit the FIMO task to the cluster node and run it in the background
    srun --mem=6G --time=5-99:00:00 --pty fimo --verbosity 2 --max-stored-scores 1000000 --thresh 0.001 --bgfile "$bg_file" --oc "$output_dir" "$meme_file" "$ref_fasta" &

done

# Wait for all background tasks to complete
wait

echo "All FIMO scan tasks for trimmed MEME files have been submitted to cluster nodes and completed."

#------------------------------------------------------------------------------------------------------

# & symbol: Allows each srun command to run in the background, enabling simultaneous task submission.
# wait command: The wait command ensures that the script waits for all background tasks to complete before proceeding. This ensures the script only exits after all FIMO tasks have finished.
