##2024-10-31
#-----------------------------------
#@ Step 1, due to  the number of Motifs is so large in ATtRACT, we split the human_pwm_chen2meme_trimmed.meme file into 100 separate files

universal_loc=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/01_snakemake_pipeline/

#Split 10 files
python $universal_loc/python_split_pwm_into_separate_files.py ./01_merged_3_10138motifs_motifs_name_replaced.ppm ./pwm_splites_10 10

#Split 100 files
python $universal_loc/python_split_pwm_into_separate_files.py ./01_merged_3_10138motifs_motifs_name_replaced.ppm ./pwm_splites 100


#@ Step 2， homer2 find 
###------分成了100个文件， 下面是单独分析每个文件----------------------------------------------
fasta_ref=/mnt/isilon/gandal_lab/mayl/reference/transcript_regions

homer2 find -i $fasta_ref/5UTR_transcript_regions.fa -o homer2_5UTR_regions_targets_1.txt -m ./motifs_part_1.pwm -mscore -p 1 -strand both 
