##2024-10-31
#-----------------------------------
#@ Step 1, due to  the number of Motifs is so large in ATtRACT, we split the human_pwm_chen2meme_trimmed.meme file into 100 separate files

universal_loc=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/01_snakemake_pipeline/

#Split 10 files
python $universal_loc/python_split_pwm_into_separate_files.py ./01_merged_3_10138motifs_motifs_name_replaced.ppm ./pwm_splites_10 10

#Split 100 files
python $universal_loc/python_split_pwm_into_separate_files.py ./01_merged_3_10138motifs_motifs_name_replaced.ppm ./pwm_splites 100


#@ Step 2， homer2 find 
###------splites into 100 files----------------------------------------------
fasta_ref=/mnt/isilon/gandal_lab/mayl/reference/transcript_regions

homer2 find -i $fasta_ref/5UTR_transcript_regions.fa -o homer2_5UTR_regions_targets_1.txt -m ./motifs_part_1.pwm -mscore -p 1 -strand both 


###-----5UTR----------------------------------------------------------------------
srun --mem=20G --time=10-00:00:00 --pty homer2 find -i /mnt/isilon/gandal_lab/mayl/reference/transcript_regions/5UTR_transcript_regions.fa -o homer2_5UTR_regions_targets_motifs_part_1.txt -m ../motifs_part_1.pwm -mscore -p 1 -strand both 


srun --mem=20G --time=10-00:00:00 --pty homer2 find -i /mnt/isilon/gandal_lab/mayl/reference/transcript_regions/5UTR_transcript_regions.fa -o homer2_5UTR_regions_targets_motifs_part_2.txt -m ../motifs_part_2.pwm -mscore -p 1 -strand both 


srun --mem=20G --time=10-00:00:00 --pty homer2 find -i /mnt/isilon/gandal_lab/mayl/reference/transcript_regions/5UTR_transcript_regions.fa -o homer2_5UTR_regions_targets_motifs_part_3.txt -m ../motifs_part_3.pwm -mscore -p 1 -strand both 


srun --mem=100G --time=10-00:00:00 --pty homer2 find -i /mnt/isilon/gandal_lab/mayl/reference/transcript_regions/5UTR_transcript_regions.fa -o homer2_5UTR_regions_targets_motifs_part_10.txt -m ../motifs_part_10.pwm -mscore -p 1 -strand both 
