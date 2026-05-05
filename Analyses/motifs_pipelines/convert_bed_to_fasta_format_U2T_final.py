#!/usr/bin/env python

import sys

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Initialize a dictionary to keep track of the name counts
    name_count = {}

    with open(output_file, 'w') as fw:
        with open(input_file, 'r') as fr:
            for line in fr:
                parts = line.strip().split()
                if len(parts) < 2:
                    continue
                
                name = parts[0]
                sequence = parts[1].replace('U', 'T')

                # Update name count
                if name in name_count:
                    name_count[name] += 1
                else:
                    name_count[name] = 1
                
                new_name = f"{name}.{name_count[name]}"

                # Write to output file in FASTA format
                fw.write(f">{new_name}\n")
                fw.write(f"{sequence}\n")
    
    print(f"Conversion to FASTA format completed. Output file: {output_file}")

if __name__ == "__main__":
    main()
