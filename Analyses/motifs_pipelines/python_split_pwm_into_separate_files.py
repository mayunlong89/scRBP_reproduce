import os
import argparse

def split_pwm_file(input_file, output_dir, num_files):
    with open(input_file, 'r') as f:
        content = f.read().split('>')[1:]  # Split content by '>' and remove the first empty element
        content = ['>' + motif for motif in content]  # Add '>' back to each motif

    num_motifs = len(content)
    motifs_per_file = num_motifs // num_files
    extra_motifs = num_motifs % num_files

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    start = 0
    for i in range(num_files):
        end = start + motifs_per_file + (i < extra_motifs)
        with open(os.path.join(output_dir, f'motifs_part_{i+1}.pwm'), 'w') as f:
            f.write(''.join(content[start:end]))
        start = end

def create_parser():
    parser = argparse.ArgumentParser(description='Split a PWM file into multiple smaller files.')
    parser.add_argument('input_file', type=str, help='Path to the input PWM file.')
    parser.add_argument('output_dir', type=str, help='Directory to save the output PWM files.')
    parser.add_argument('num_files', type=int, help='Number of output PWM files to create.')
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()

    split_pwm_file(args.input_file, args.output_dir, args.num_files)

if __name__ == '__main__':
    main()

