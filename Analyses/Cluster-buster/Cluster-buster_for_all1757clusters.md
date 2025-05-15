
### 1.1 Note
```R
##------------------------------------------------------------------------------------
# 2024-12-20
##------------------------------------------------------------------------------------
# For 8578 clustered motifs to obtain 164 initial clusters and 1575 refined sub-clusters, we leveraged Cluster-Buster to run clustered motifs three strategies:
# 1) Use all motifs in each subculster from 1575 subclusters to run Cluster-buster method;
# 2) Use archetype motif in each subcluster to run Cluster-buster method;
# 3) Use all 8578 motifs to individually run Cluster-buster method;

# For other singletons, we also individually run Cluster-buster method, which is same with the strategy 3.

# Meanwhile, for benchmark, we also use all 25044 motifs to individually run Fimo and HOMER2. 
##------------------------------------------------------------------------------------
```

### 1.2 running for 5UTR
```bash
##------------------------------------------------------------------------------------
# 2024-12-11
##------------------------------------------------------------------------------------

##`PCA-Seurat-RSAT Workflow for Motif Clustering`
#1) Clustering 25,044 Motifs
Using PCA-Seurat-RSAT, we obtained 1,575 refined clusters from 25,044 motifs.

These clusters were further analyzed using Cluster-Buster to identify regulatory modules.

#2) Extracting Motif IDs for Each Cluster
Motif IDs for each cluster have been extracted using custom scripts.

The resulting files are saved as Cluster9_cluster_01.txt.

#3) Working Directory
Location:
/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/08_HOMER/matrix-clustering_stand-alone/01_RBP_616RBPs_25044motifs/01_8758motifs_clustered/clusters/refined_clusters/

#4) RSAT Matrix Clustering
Using RSAT matrix clustering, we refined the 164 main clusters into 1,575 refined clusters.

# Each main cluster generates a key result file named Cluster99_motifID_motifs_rsat_motifs.


mkdir -p zzz_all_motifs
cd zzz_all_motifs

# remove '_motifID_motifs_rsat_motifs'
# convert Cluster9_motifID_motifs_rsat_motifs >> Cluster9
python convert_dir_names.py 

# convert Cluster9 >> Cluster9_cluster_10
python convert_subdir_names.py

# output: 'merge_all_motifs_cluster'
python merge_all_clusters_into_oneDir.py

# process the cluster-motif .txt files
# folder: Cluster9_cluster_43 >> Cluster9_cluster_43.txt
cd merge_all_motifs_cluster
bash process_motifs_ID.sh 

# move cluster-motifs into cluster_ID folder
cp *.txt /mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/03_cluster_buster/cluster_ID/


#calculate the count of motifs in each cluster
cd /mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/03_cluster_buster/cluster_ID/

python process_count_motifs.py


#--use 'create_annotation_metadata.py' to generate a annotation file
#Output: create_annotation_file.tsv
#Columns: Sequence_ID, Motif_ID, Cluster_ID
 
 python create_annotation_metadata.py

#2), clustered motif .ppm file extracted from all 25044 motifs.
#using 'extract_pwm_from_allMotifs.py' to extract motif.ppm

python extract_pwm_from_allMotifs.py -i 04_merged_all_25044motifs_all_databases.pwm -l ./cluster_ID/Cluster9_cluster_01.txt -o ./cluster_pwm/Cluster9_cluster_01.ppm


#3), using `convert_pwm_to_cluster.py` transform .ppm file scale to the format of Cluster-buster

python convert_pwm_to_cluster.py


#4) ------------Cluster-buster methods
#-g 35 default
#After considering the RBP motifs are smaller than TFs' motifs, we here used -g 20 for identifying clusters of motifs

cbust  -g 20 -l -f 5 -c 3 ./cluster_pwm_for_cbuster/Cluster9_cluster_01.ppm /mnt/isilon/gandal_lab/mayl/reference/transcript_regions/5UTR_transcript_regions.fa > ./Cluster9_cluster_01_output.txt

#batch--running
# see codes in file 'batch_clusterbuster_5UTR.py'
python batch_clusterbuster_5UTR.py

```

### 1.3 Clustered motifs---> cluster annotation metadata.py

```bash
#2024-12-22

#Working location:
/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/03_cluster_buster/cluster_ID/

#output: create_annotation_file.tsv

python create_annotation_metadata.py

```







