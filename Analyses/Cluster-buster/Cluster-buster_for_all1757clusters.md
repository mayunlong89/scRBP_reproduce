
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

#Each main cluster generates a key result file named Cluster99_motifID_motifs_rsat_motifs.
```
