#Location: 
matrixID=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/CISBP_RNA/Human_2024_07_05/TomTom_M354_0.6

#compared with streme-identified motifs
#location:

location=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/MotifMap-rna/archived/test_for_M354_05_Motif_fimo_pipeline/M354_0.6_motifs_dir

tomtom -oc tomtom_output_dir $location/trimmed_meme_file.meme   $matrixID/M354_0.6.meme
