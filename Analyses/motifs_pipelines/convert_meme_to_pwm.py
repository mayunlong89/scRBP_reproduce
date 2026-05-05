import sys

def convert_meme_to_ppm(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        in_motif = False
        for line in infile:
            # Identify the motif line
            if line.startswith("MOTIF"):
                parts = line.split()
                motif_id = parts[1]  # Extract the simplified motif ID
                outfile.write(f">{motif_id}\n")
                in_motif = True
            elif in_motif and line.strip():
                # Skip the `letter-probability matrix` line
                if not line.startswith("letter-probability matrix"):
                    outfile.write(line)
            elif not line.strip():  # Stop processing if a blank line is found
                in_motif = False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    convert_meme_to_ppm(input_file, output_file)
    print(f"Converted MEME file '{input_file}' to PPM file '{output_file}' without 'letter-probability matrix' lines.")

