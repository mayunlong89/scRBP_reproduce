
`2025-05-15`

## RSAT matrix clustering for subclustering these motifs

### 1.1 activiating conda environment

```bash

conda activate rsat_matrix_clustering

```

### 1.2 install pacakges

```bash

#version 1: Installing rsat-core

conda install -c bioconda rsat-core


#version 2: matrix_clustering stand-alone (current study based on this version)
#https://github.com/jaimicore/matrix-clustering_stand-alone

Rscript matrix-clustering.R                           \
  -i data/OCT4_datasets/OCT4_motif_table.txt          \
  -o results/OCT4_motifs_clusters/OCT4_motif_analysis \
  -w 1


```

### 1.3 format conversion
```bash

#1.3.1 get .tf format for rsat matrix clustering analyses

convert-matrix -from meme -to tf input.meme > output_rsat.tf #used for rsat matrix clustering

#1.3.2 transfac_rsat convert to meme

rsat=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/08_HOMER/matrix-clustering_stand-alone

Rscript $rsat/convert-matrix.R                    \
  -i Root_motifs.tf   \
  --from tf --to meme                   \
  --output_file Root_motifs.meme



```

### 1.4 Main running
```bash
#Step 1 extract Motifs ID in each main cluster
#'motifs_clusters_10-3_resolution15.txt'--motif-cluster annotation
cat motifs_clusters_10-3_resolution15.txt | awk '$3=="Cluster29"' | cut -f 1 > cluster29_motifID.txt

# Step 2 use the 'extract_subset_motifs_v2.py' code to subset the .pwm file from '04_merged_all_25044motifs_all_databases.pwm'

python extract_subset_motifs_v2.py -p 04_merged_all_25044motifs_all_databases.pwm \
                         -m Cluster29_motifID.txt \
                         -o subset_cluster29_motifs.ppm

# Step 3 use 'chen2meme' transform .ppm to .meme file

chen2meme subset_cluster29_motifs.ppm > subset_cluster29_motifs.meme

# Step 4 use 'convert-matrix.R' to convert meme format into transfac(tf) format
rsat=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/08_HOMER/matrix-clustering_stand-alone

Rscript $rsat/convert-matrix.R                    \
  -i subset_cluster29_motifs.meme    \
  --from meme --to tf                     \
  --output_file subset_cluster29_motifs.tf


# Step 5 use 'transfac_to_standard_rsat.py' to convert transfac format to standard rsat format
# This step is very important, and the script is in-house generated.
rsat=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/08_HOMER/matrix-clustering_stand-alone

python $rsat/transfac_to_standard_rsat.py -i subset_cluster29_motifs.tf -o subset_cluster29_motifs_rsat.tf


#put 'subset_cluster29_motifs_rsat.tf' into the 'subset_cluster29_motifs.txt' file
#Format: test/subset_cluster29_motifs_rsat.tf MEME tf


# Step 6 use 'matrix-clustering.R' to cluster analysis the motifs
rsat=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/08_HOMER/matrix-clustering_stand-alone
Rscript $rsat/matrix-clustering.R -i ./test/subset_cluster29_motifs.txt -o ./test/subset_cluster29 -w 1


rsat=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/08_HOMER/matrix-clustering_stand-alone
Rscript $rsat/matrix-clustering.R -i ./test/combined_clusters.txt -o ./test/combined_clusters -w 1


rsat=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/08_HOMER/matrix-clustering_stand-alone
Rscript $rsat/matrix-clustering.R -i ./all_motifs/allmotifs.txt -o ./all_motifs/combined_clusters -w 1


##-----final version results: 'refined_clusters'



# Step 7 use '--radial_tree' parameter for generating plots

# 「RBP_motifs_annotation.txt」 ##This file is required as an RBP-motifs classification annotation file.

##--It will be used in downstream analyses, including phylogenetic tree construction for the 25,044 motifs.
##Additionally, it will support RBP functional annotation for future studies.

Rscript matrix-clustering.R                           \
  -i data/JASPAR_2022/Jaspar_nematodes_motifs_tab.txt \
  -o results/JASPAR_nematodes_radial/JASPAR_nematodes \
  -a data/JASPAR_2022/JASPAR_nematodes_metadata.txt   \
  --radial_tree TRUE                                  \
  --title JASPAR_CORE_nematodes                       \
  -w 8
  
rsat=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/08_HOMER/matrix-clustering_stand-alone
Rscript $rsat/matrix-clustering.R -i ./test/subset_cluster29_motifs.txt -o ./test/subset_cluster29 --radial_tree TRUE  -a RBP_motifs_annotation.txt   -w 1

```

### 1.5 RSAT for batch run 164 initial clusters based on top 100 PCs from motif-by-motif similarity 

```bash

#Working location

/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/08_HOMER/matrix-clustering_stand-alone/01_RBP_616RBPs_25044motifs/01_8758motifs_clustered/clusters


#Step 1
python split_clusters.py

#Step 2
python batch_process_motifs.py

#Step 3
python batch_convert_to_meme.py

#Step 4
python batch_convert_meme_to_tf.py

#Step 5
python batch_convert_to_standard_rsat.py

#Step 6
python generate_cluster_txt_files.py

#Step 7
python batch_run_matrix_clustering.py


```




### 1.6 RSAT result folder structure

```bash
###---RAST matrix-clustering

results
├── *_motifs
│   ├── individual_motifs_with_gaps
│   │   └── Two files per motif (direct and reverse orientation) in transfac format, these motifs are already aligned (may contain gaps).
│   │
│   ├── motifs_sep_by_cluster   (each folder contains the motifs belonging to a cluster)
│   │   └── Cluster_01
│   │   └── Cluster_02
│   │   └── ...
│   │   └── Cluster_N
│   │
│   └── root_motifs
│       └── Root_motifs.tf  (also referred as archetype motifs)
│
├── *_plots
│   ├── Clusters_vs_reference_contingency_table.pdf
│   └── Heatmap_clusters.pdf
│   
├── *_tables
│   ├── alignment_table.tab
│   ├── clusters.tab
│   ├── distance_table.tab
│   ├── pairwise_motif_comparison.tab
│   └── summary_table.tab
│
└── *_trees
    ├── tree.json
    ├── tree.newick
    └── tree.RData

```













































