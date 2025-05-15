import os

def generate_cluster_txt_files(input_dir, output_dir):
    """
    Generate .txt files for each _rsat.tf file with specific content format.
    
    Args:
        input_dir (str): Directory containing the input _rsat.tf files.
        output_dir (str): Directory where the output .txt files will be saved.
    """
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Iterate over all _rsat.tf files in the input directory
    for rsat_file in os.listdir(input_dir):
        if rsat_file.endswith("_rsat.tf"):  # Only process _rsat.tf files
            cluster_name = os.path.splitext(rsat_file)[0]  # Get file name without extension
            input_path = os.path.join(input_dir, rsat_file)
            output_file = os.path.join(output_dir, f"{cluster_name}.txt")
            
            # Prepare content for the .txt file
            content = f"{input_path}\t{cluster_name}\ttf\n"
            
            # Write content to the new .txt file
            with open(output_file, "w") as f:
                f.write(content)
            
            print(f"Generated: {output_file}")

    print("All .txt files have been generated.")


if __name__ == "__main__":
    # Specify input and output directories
    input_dir = "standard_rsat_motifs"  # Directory containing the _rsat.tf files
    output_dir = "cluster_txt_files"  # Directory to save the .txt files
    
    # Run the script
    generate_cluster_txt_files(input_dir, output_dir)
