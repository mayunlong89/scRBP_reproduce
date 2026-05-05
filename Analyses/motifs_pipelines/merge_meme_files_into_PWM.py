#!/usr/bin/env python

import os
import sys
import argparse
import glob

def main():
    description = 'Merge multiple MEME files into one PWM file with each motif title prefixed by ">"'
    parser = argparse.ArgumentParser(description=description)
    
    parser.add_argument('idir', help='input directory containing MEME files')
    parser.add_argument('output_file', help='specify the output file path')
    
    args = parser.parse_args(sys.argv[1:])
    
    idir = args.idir
    output_file = args.output_file
    
    merge_meme_files(idir, output_file)


def merge_meme_files(idir, output_file):
    meme_files = glob.glob(os.path.join(idir, "*.meme"))
    nmotif = 0
    
    with open(output_file, "w") as fw:
        for meme_file in meme_files:
            with open(meme_file) as fh:
                motif_section = False
                for line in fh:
                    line = line.strip()
                    if line.startswith("MOTIF"):
                        nmotif += 1
                        motif_section = True
                        print(f">{line.split()[1]}", file=fw)
                    elif line.startswith("letter-probability matrix"):
                        motif_section = True
                        continue  # Skip this line
                    elif motif_section:
                        if line == "":
                            motif_section = False
                            print(file=fw)
                        else:
                            print(line, file=fw)

    print(f"Combined motifs from {len(meme_files)} MEME files into {output_file}")


if __name__ == "__main__":
    main()
