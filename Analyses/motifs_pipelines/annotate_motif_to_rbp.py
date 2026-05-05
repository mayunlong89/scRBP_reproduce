import pandas as pd
import argparse

def annotate_motif_to_rbp(motif_id_list_file, rbp_motif_mapping_file, output_file):
    # Load RBP-motif mapping file
    rbp_motif_mapping = pd.read_csv(rbp_motif_mapping_file, sep='\t')

    # Load motif ID list file
    motif_id_list = pd.read_csv(motif_id_list_file, header=None, names=['MotifID'])

    # Merge the motif ID list with the RBP-motif mapping file to handle many-to-many relationships
    annotated_rbp = pd.merge(motif_id_list, rbp_motif_mapping, on='MotifID', how='left')

    # Group by MotifID and aggregate corresponding RBPs into a list
    annotated_rbp_grouped = annotated_rbp.groupby('MotifID').agg({'RBP': lambda x: ','.join(set(x.dropna()))}).reset_index()

    # Output the result to a file
    annotated_rbp_grouped.to_csv(output_file, sep='\t', index=False)

    print(f"Annotation results have been saved to '{output_file}'")

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Annotate motif IDs to corresponding RBPs.")
    parser.add_argument("motif_id_list_file", type=str, help="Path to the motif ID list file")
    parser.add_argument("rbp_motif_mapping_file", type=str, help="Path to the RBP-motif mapping file")
    parser.add_argument("output_file", type=str, help="Path to the output file")

    args = parser.parse_args()

    # Call the function to perform the annotation
    annotate_motif_to_rbp(args.motif_id_list_file, args.rbp_motif_mapping_file, args.output_file)

