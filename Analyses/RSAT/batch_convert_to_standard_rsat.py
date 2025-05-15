import os
import subprocess

def batch_convert_to_standard_rsat(input_dir, output_dir, transfac_script_path):
    """
    Batch convert TF format motif files to standard RSAT format using transfac_to_standard_rsat.py.
    
    Args:
        input_dir (str): Directory containing input TF motif files.
        output_dir (str): Directory where the converted RSAT files will be saved.
        transfac_script_path (str): Path to the transfac_to_standard_rsat.py script.
    """
    # Create the output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Iterate over each TF file in the input directory
    for tf_file in os.listdir(input_dir):
        if tf_file.endswith(".tf"):  # Process only .tf files
            input_file = os.path.join(input_dir, tf_file)
            output_file = os.path.join(output_dir, tf_file.replace(".tf", "_rsat.tf"))
            
            # Construct the command
            command = [
                "python", transfac_script_path,
                "-i", input_file,
                "-o", output_file
            ]
            
            # Run the command
            print(f"Converting {tf_file} to {output_file}")
            subprocess.run(command, check=True)
    
    print("Batch conversion to standard RSAT format completed.")


if __name__ == "__main__":
    # Specify paths and directories
    input_dir = "tf_motifs"  # Directory containing input TF motif files
    output_dir = "standard_rsat_motifs"  # Directory for the converted RSAT files
    transfac_script_path = "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/08_HOMER/matrix-clustering_stand-alone/transfac_to_standard_rsat.py"  # Path to the transfac_to_standard_rsat.py script
    
    # Run batch conversion
    batch_convert_to_standard_rsat(input_dir, output_dir, transfac_script_path)
