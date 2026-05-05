import os
import subprocess
import argparse

# Function to perform HOMER motif analysis on all BED files in a given directory
def run_homer_on_bed_files(bed_files_directory):
    # Define the genome version and parameters for HOMER
    genome_version = "hg38"
    size_option = "given"
    lengths = "6,8,10"

    # Check if the directory exists
    if not os.path.exists(bed_files_directory):
        print(f"Directory {bed_files_directory} does not exist.")
        return

    # Loop over each file in the directory
    for filename in os.listdir(bed_files_directory):
        if filename.endswith("_filtered.txt"):  # Filter to process only your bed files
            rbp_name = filename.split("_filtered.txt")[0]
            bed_file_path = os.path.join(bed_files_directory, filename)
            output_directory = os.path.join(bed_files_directory, f"output_directory_{rbp_name}")

            # Create the command to run HOMER for the current file
            homer_command = [
                "findMotifsGenome.pl", 
                bed_file_path, 
                genome_version, 
                output_directory, 
                "-size", size_option, 
                "-len", lengths
            ]

            # Run the HOMER command
            print(f"Running HOMER for {rbp_name}...")
            subprocess.run(homer_command)

    print("Batch processing completed.")

# Main function to accept directory path as an argument
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run HOMER motif analysis on all filtered BED files.")
    parser.add_argument("directory", help="The directory path containing the BED files.")
    
    args = parser.parse_args()
    
    # Call the function with the provided directory argument
    run_homer_on_bed_files(args.directory)

