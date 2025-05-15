import os
import subprocess

def batch_convert_meme_to_tf(input_dir, output_dir, rsat_script_path):
    """
    Batch convert MEME format motif files to TF format using RSAT convert-matrix.R.
    
    Args:
        input_dir (str): Directory containing input MEME motif files.
        output_dir (str): Directory where the converted TF files will be saved.
        rsat_script_path (str): Path to the RSAT convert-matrix.R script.
    """
    # Create the output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Iterate over each MEME file in the input directory
    for meme_file in os.listdir(input_dir):
        if meme_file.endswith(".meme"):  # Process only .meme files
            input_file = os.path.join(input_dir, meme_file)
            output_file = os.path.join(output_dir, meme_file.replace(".meme", ".tf"))
            
            # Construct the command
            command = [
                "Rscript", rsat_script_path,
                "-i", input_file,
                "--from", "meme",
                "--to", "tf",
                "--output_file", output_file
            ]
            
            # Run the command
            print(f"Converting {meme_file} to {output_file}")
            subprocess.run(command, check=True)
    
    print("Batch conversion to TF format completed.")


if __name__ == "__main__":
    # Specify paths and directories
    input_dir = "meme_motifs"  # Directory containing input MEME motif files
    output_dir = "tf_motifs"  # Directory for the converted TF files
    rsat_script_path = "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/08_HOMER/matrix-clustering_stand-alone/convert-matrix.R"  # Path to RSAT convert-matrix.R script
    
    # Run batch conversion
    batch_convert_meme_to_tf(input_dir, output_dir, rsat_script_path)
