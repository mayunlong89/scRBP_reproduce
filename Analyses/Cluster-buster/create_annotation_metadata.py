import os

# Set the base directory containing all the .txt files
base_dir = "./"
output_file = "create_annotation_file.tsv"

# Initialize variables
sequence_id = 1  # Start sequence ID counter
annotations = []  # Store annotation data

# Iterate over all .txt files in the directory
for file_name in sorted(os.listdir(base_dir)):
    if file_name.endswith(".txt"):
        # Get the cluster ID by removing ".txt" from the file name
        cluster_id = file_name.replace(".txt", "")

        # Read motifs from the file
        file_path = os.path.join(base_dir, file_name)
        with open(file_path, "r") as f:
            motifs = [line.strip() for line in f if line.strip()]

        # Create annotations for each motif
        for motif in motifs:
            annotations.append([sequence_id, motif, cluster_id])
            sequence_id += 1

# Write annotations to the output file
with open(output_file, "w") as f:
    f.write("Sequence_ID\tMotif_ID\tCluster_ID\n")  # Header
    for annotation in annotations:
        f.write("\t".join(map(str, annotation)) + "\n")

print(f"Annotation file created: {output_file}")
