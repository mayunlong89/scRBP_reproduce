import os
import subprocess

def batch_convert_to_meme(input_dir, output_dir, chen2meme_path="chen2meme"):
    """
    Batch convert PPM format motif files to MEME format using chen2meme.
    
    Args:
        input_dir (str): Directory containing input PPM motif files.
        output_dir (str): Directory where the converted MEME files will be saved.
        chen2meme_path (str): Path to the chen2meme executable. Default is "chen2meme".
    """
    # Create the output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Iterate over each PPM file in the input directory
    for ppm_file in os.listdir(input_dir):
        if ppm_file.endswith(".ppm"):  # Process only .ppm files
            input_file = os.path.join(input_dir, ppm_file)
            output_file = os.path.join(output_dir, ppm_file.replace(".ppm", ".meme"))
            
            # Construct the command
            command = f"{chen2meme_path} {input_file} > {output_file}"
            
            # Run the command
            print(f"Converting {ppm_file} to {output_file}")
            subprocess.run(command, shell=True, check=True)
    
    print("Batch conversion to MEME format completed.")


if __name__ == "__main__":
    # Specify paths and directories
    input_dir = "subset_motifs"  # Directory containing input PPM motif files
    output_dir = "meme_motifs"  # Directory for the converted MEME files
    chen2meme_path = "chen2meme"  # Path to the chen2meme executable
    
    # Run batch conversion
    batch_convert_to_meme(input_dir, output_dir, chen2meme_path)
