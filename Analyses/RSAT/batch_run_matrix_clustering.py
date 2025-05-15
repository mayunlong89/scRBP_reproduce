import os
import subprocess

def batch_run_matrix_clustering(input_dir, output_dir, rsat_script_path):
    """
    Batch run RSAT matrix-clustering.R for each input file.
    
    Args:
        input_dir (str): Directory containing input RSAT TF motif files.
        output_dir (str): Directory to save refined cluster results.
        rsat_script_path (str): Path to RSAT matrix-clustering.R script.
    """
    # Ensure the input directory exists
    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    
    # Ensure the RSAT script path exists
    if not os.path.exists(rsat_script_path):
        raise FileNotFoundError(f"RSAT script does not exist: {rsat_script_path}")

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    for tf_file in os.listdir(input_dir):
        if tf_file.endswith(".txt"):  # Process only .tf files
            input_file = os.path.join(input_dir, tf_file)  # Full path to the .tf file
            cluster_name = os.path.splitext(tf_file)[0]  # Extract cluster name without extension
            cluster_output_dir = os.path.join(output_dir, cluster_name)  # Create output directory for the cluster
            os.makedirs(cluster_output_dir, exist_ok=True)
            
            # Construct the command
            command = [
                "Rscript", rsat_script_path,
                "-i", input_file,
                "-o", cluster_output_dir,
                "-w", "1"
            ]
            
            try:
                # Run the command
                print(f"Running matrix-clustering for {tf_file}")
                subprocess.run(command, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error occurred while processing {tf_file}: {e}")
                continue  # Continue processing the next file
    
    print("Batch matrix-clustering completed.")

if __name__ == "__main__":
    # Specify paths and directories
    input_dir = "cluster_txt_files"  # Directory containing input RSAT TF motif files
    output_dir = "refined_clusters"  # Directory for saving refined clusters
    rsat_script_path = "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/08_HOMER/matrix-clustering_stand-alone/matrix-clustering.R"
    
    # Run batch processing
    batch_run_matrix_clustering(input_dir, output_dir, rsat_script_path)


