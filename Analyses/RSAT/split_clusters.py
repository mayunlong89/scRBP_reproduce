import os
import pandas as pd

def split_clusters_to_files(input_file, output_dir):
    """
    Split a file containing motif IDs and clusters into separate files by cluster.
    
    Args:
        input_file (str): Path to the input file (tab-delimited with MotifID and Cluster columns).
        output_dir (str): Directory where the output files will be saved.
    """
    # Create output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Load the data
    data = pd.read_csv(input_file, sep="\t")
    
    # Check if required columns exist
    if "MotifID" not in data.columns or "ClusterID" not in data.columns:
        raise ValueError("Input file must contain 'MotifID' and 'ClusterID' columns.")
    
    # Group by ClusterID and save each cluster's MotifID to a separate file
    for cluster_id, group in data.groupby("ClusterID"):
        output_file = os.path.join(output_dir, f"{cluster_id}_motifID.txt")
        group["MotifID"].to_csv(output_file, index=False, header=False)
        print(f"Cluster '{cluster_id}' motifs saved to: {output_file}")


if __name__ == "__main__":
    # Specify input file and output directory
    input_file = "motifs_clusters_10-3_resolution15.txt"  # Replace with your input file
    output_dir = "clustered_motifs"  # Replace with your desired output directory

    # Run the function
    split_clusters_to_files(input_file, output_dir)
