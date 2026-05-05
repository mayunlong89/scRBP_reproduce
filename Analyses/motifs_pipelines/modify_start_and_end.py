
#2024-07-25
# 要将 BED 文件中的 start 和 end 位置保持一致，并且确保 start 总是在 end 之前，可以编写一个简单的 Python 脚本。
# 这个脚本会读取 BED 文件，检查每一行中的 start 和 end，如果发现 start 大于 end，则交换它们。

import argparse

def create_parser():
    parser = argparse.ArgumentParser(description="Normalize start and end positions in a BED file.")
    parser.add_argument('input_bed', help="Input BED file")
    parser.add_argument('output_bed', help="Output BED file")
    return parser

def normalize_bed(input_bed, output_bed):
    with open(input_bed, 'r') as infile, open(output_bed, 'w') as outfile:
        for line in infile:
            if line.startswith('track') or line.startswith('browser'):
                # Skip header lines
                outfile.write(line)
                continue

            fields = line.strip().split('\t')
            if len(fields) < 6:
                # Skip lines that don't have enough columns (at least 6 columns needed for strand info)
                continue

            chrom, start, end, name, score, strand = fields[:6]
            start, end = int(start), int(end)

            if start > end:
                start, end = end, start

            # Write normalized line to the output file
            outfile.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    normalize_bed(args.input_bed, args.output_bed)

