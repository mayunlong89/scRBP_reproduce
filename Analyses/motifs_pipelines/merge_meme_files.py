import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Merge multiple MEME files into one.')
    parser.add_argument('output_file', help='The output MEME file.')
    parser.add_argument('input_files', nargs='+', help='Input MEME files to be merged.')
    return parser.parse_args()

def merge_meme_files(output_file, input_files):
    with open(output_file, 'w') as outfile:
        outfile.write("MEME version 4\n\n")  # Write the MEME header
        alphabet_written = False
        background_written = False

        for file in input_files:
            with open(file, 'r') as infile:
                in_motif_section = False
                for line in infile:
                    if line.startswith('ALPHABET='):
                        if not alphabet_written:
                            outfile.write(line)
                            alphabet_written = True
                    elif line.startswith('Background letter frequencies'):
                        if not background_written:
                            outfile.write(line)
                            background_written = True
                    elif line.startswith('MOTIF'):
                        in_motif_section = True
                        outfile.write(line)
                    elif in_motif_section:
                        outfile.write(line)
                    elif line.strip() == '':
                        continue

if __name__ == '__main__':
    args = parse_arguments()
    merge_meme_files(args.output_file, args.input_files)
