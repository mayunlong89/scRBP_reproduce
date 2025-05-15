import os
import subprocess

def batch_process_motifs(input_dir, pwm_file, output_dir, script_path):
    """
    Batch process motif ID files by running extract_subset_motifs_v2.py for each motif file.
    
    Args:
        input_dir (str): Directory containing the motif ID files (e.g., Cluster29_motifID.txt).
        pwm_file (str): Path to the PPM file containing all motifs.
        output_dir (str): Directory where the subset motifs files will be saved.
        script_path (str): Path to the Python script extract_subset_motifs_v2.py.
    """
    # Create the output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Iterate over each motif ID file in the input directory
    for motif_file in os.listdir(input_dir):
        if motif_file.endswith(".txt"):  # Process only .txt files
            cluster_name = os.path.splitext(motif_file)[0]  # Extract the cluster name (e.g., Cluster29)
            input_motif_file = os.path.join(input_dir, motif_file)
            output_file = os.path.join(output_dir, f"{cluster_name}_motifs.ppm")
            
            # Construct the command
            command = [
                "python", script_path,
                "-p", pwm_file,
                "-m", input_motif_file,
                "-o", output_file
            ]
            
            # Run the command
            print(f"Processing {motif_file} -> {output_file}")
            subprocess.run(command, check=True)
    
    print("Batch processing completed.")


if __name__ == "__main__":
    # Specify paths and directories
    input_dir = "clustered_motifs"  # Directory containing motif ID files
    pwm_file = "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/08_HOMER/matrix-clustering_stand-alone/01_RBP_616RBPs_25044motifs/01_8758motifs_clustered/04_merged_all_25044motifs_all_databases.pwm"  # Path to the PPM file
    output_dir = "subset_motifs"  # Directory for the output subset motif files
    script_path = "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/08_HOMER/matrix-clustering_stand-alone/01_RBP_616RBPs_25044motifs/01_8758motifs_clustered/extract_subset_motifs_v2.py"  # Path to the Python script for extracting motifs
    
    # Run batch processing
    batch_process_motifs(input_dir, pwm_file, output_dir, script_path)
