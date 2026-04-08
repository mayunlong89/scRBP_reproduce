### For CDS (Gene and Isoform), 5'UTR (Gene and Isoform), 3'UTR (Gene and Isoform), Intron (Gene), we used this version.

```bash
python convert_fimo_to_feather.py --input_file fimo_best_site_narrowPeak_targets_rankedByPvalue_ranked.tsv --input_dir /path/to/input/dir --output_dir /path/to/output/dir
```

For more details on `convert_fimo_to_feather.py`, please refer to see [here](https://github.com/mayunlong89/scRBP_reproduce/blob/main/Analyses/convert_fimo_to_feather/convert_fimo_to_feather.py)


### Example

```bash
#location:-----Introns--
#location:-----Introns--

#----	Remove version numbers from ENST transcript IDs.
final_intron=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings

awk 'BEGIN{OFS="\t"} {sub(/\.[0-9]+$/, "", $6); print}' $final_intron/merged_all_20746motifs_for_final_motif_isoforms_rankings_Introns_ranked.bed > $final_intron/merged_all_20746motifs_for_final_motif_isoforms_rankings_Introns_ranked_noVersion.bed


#------convert .bed into .feather
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/CDS_cluster_buster

final_intron=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings
bed_intron=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings

mkdir -p $final_intron/Cluster_Buster_matrix_Introns_isoform_final_20476_v2


srun --mem=850G --cpus-per-task=20 --time=30-00:00:00 \
  --export=ALL,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1,NUMEXPR_MAX_THREADS=1 \
  --pty python $TOOL/convert_fimo_to_feather.py \
  --input_file $bed_intron/merged_all_20746motifs_for_final_motif_isoforms_rankings_Introns_ranked_noVersion.bed \
  --input_dir $bed_intron \
  --output_dir $final_intron/Cluster_Buster_matrix_Introns_isoform_final_20476_v2

```

### For Intron (isoform), we use an enhanced version due to the file is so big.
```bash

TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/CDS_cluster_buster

final_intron=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings
bed_intron=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/10_final_version_motif_rankings

mkdir -p $final_intron/Cluster_Buster_matrix_Introns_isoform_final_20476_v2
srun --mem=850G --cpus-per-task=20 --time=30-00:00:00 \
  --export=ALL,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1,NUMEXPR_MAX_THREADS=1 \
  --pty python $TOOL/convert_fimo_to_feather_enhanced.py \
    --input_file merged_all_20746motifs_for_final_motif_isoforms_rankings_Introns_ranked_noVersion.bed \
    --input_dir $bed_intron \
    --output_dir $final_intron/Cluster_Buster_matrix_Introns_isoform_final_20746_v2

````
For more details on `convert_fimo_to_feather_enhanced`, please refer to see [here](https://github.com/mayunlong89/scRBP_reproduce/blob/main/Analyses/convert_fimo_to_feather/convert_fimo_to_feather_enhanced.py)




