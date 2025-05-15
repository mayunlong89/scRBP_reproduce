import os
import subprocess

# Set directories
ppm_dir = "./cluster_pwm_for_cbuster"  # Directory containing 1,575 PPM files
output_dir = "./5UTR"                  # Directory for output files
fasta_file = "/mnt/isilon/gandal_lab/mayl/reference/transcript_regions/5UTR_transcript_regions.fa"  # Input FASTA file
cbust_command = "cbust"                # Path to the Cluster-Buster command

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Iterate over all PPM files
for ppm_file in os.listdir(ppm_dir):
    if ppm_file.endswith(".ppm"):  # Process only .ppm files
        cluster_name = os.path.splitext(ppm_file)[0]  # Remove the file extension to get the cluster name
        ppm_path = os.path.join(ppm_dir, ppm_file)
        output_path = os.path.join(output_dir, f"{cluster_name}_output.txt")
        
        # Construct the Cluster-Buster command
        command = f"{cbust_command} -g 20 -l -f 5 -c 3 {ppm_path} {fasta_file} > {output_path}"
        
        # Execute the command
        try:
            subprocess.run(command, shell=True, check=True)
            print(f"Processed: {ppm_file} -> {output_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error processing {ppm_file}: {e}")
