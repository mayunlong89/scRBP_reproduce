#2024-10-31

#working location:

/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs/002_fimo_reuslts

#========================================
#@ Step 1, due to  the number of Motifs is so large in ATtRACT, we split the human_pwm_chen2meme_trimmed.meme file into 100 separate files
universal_loc=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/01_snakemake_pipeline/


python $universal_loc/python_split_pwm_into_separate_files.py ../01_merged_3_10138motifs_all_databases.pwm ./pwm_splites 100

#========================================
#@ Step 2 using batch_convert bash code to use 'chen2meme' to convert .pwm file into .meme format
universal_loc=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/01_snakemake_pipeline/

bash $universal_loc/convert_pwm_to_meme.sh ./pwm_splites ./pwm_convert_chen2meme_splites


#========================================
#@ Step 3 using universalmotif to refine the meme files

universal_loc=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/01_snakemake_pipeline/




bash $universal_loc/trim_meme_files.sh ./pwm_convert_chen2meme_splites ./pwm_convert_chen2meme_trimmed_splites

#universal processed files
Rscript $universal_loc/universalmotif.r human_pwm_chen2meme.meme human_pwm_chen2meme_trimmed.meme


#========================================
#@ Step 4 using Fimo to search motifs' targets across the genome.
#Using FIMO to scan input FASTA sequence with background

#ref=/mnt/isilon/gandal_lab/mayl/reference/
#fimo --verbosity 2 --max-stored-scores 1000000 --thresh 0.001 --oc fimo_output_dir final_trimmed_pwm.txt $ref/gene_sequence_v44_GCRh38.fasta

#2. based on background model to run FIMO for scanning input FASTA sequence
ref=/mnt/isilon/gandal_lab/mayl/reference/
bgfile=/mnt/isilon/gandal_lab/mayl/reference/

fimo --verbosity 2 --max-stored-scores 1000000 --thresh 0.001 --bgfile $bgfile/background.meme --oc fimo_output_dir_bg trimmed_refined_bamm2meme_file.meme $ref/gene_sequence_v44_GCRh38.fasta



###5UTR--motifs_part_1 to motifs_part_10
ref=/mnt/isilon/gandal_lab/mayl/reference/
bgfile=/mnt/isilon/gandal_lab/mayl/reference/

srun --mem=120G --pty --time=10-00:00:00 fimo --verbosity 2 --max-stored-scores 1000000 --thresh 0.001 --bgfile $bgfile/background.meme --oc fimo_output_dir_bg_part_10 motifs_part_10_trimmed.meme $ref/transcript_regions/5UTR_transcript_regions.fa


##3UTR motifs_part_1 to motifs_part_10
ref=/mnt/isilon/gandal_lab/mayl/reference/
bgfile=/mnt/isilon/gandal_lab/mayl/reference/

srun --mem=120G --pty --time=10-00:00:00 fimo --verbosity 2 --max-stored-scores 1000000 --thresh 0.001 --bgfile $bgfile/background.meme --oc fimo_output_dir_bg_part_1_3UTR motifs_part_1_trimmed.meme $ref/transcript_regions/3UTR_transcript_regions.fa



###CDS motifs_part 1 to motif part 10 

ref=/mnt/isilon/gandal_lab/mayl/reference/
bgfile=/mnt/isilon/gandal_lab/mayl/reference/

srun --mem=350G --pty --time=10-00:00:00 fimo --verbosity 2 --max-stored-scores 1000000 --thresh 0.001 --bgfile $bgfile/background.meme --oc fimo_output_dir_bg_part_1_CDS ../motifs_part_1_trimmed.meme $ref/transcript_regions/CDS_transcript_regions.fa




###Intron motifs_part 1 to motif part 10  

srun --mem=350G --pty --time=15-00:00:00 fimo --verbosity 2 --max-stored-scores 1000000 --thresh 0.001 --bgfile $bgfile/background.meme --oc fimo_output_dir_bg_part_1_introns ../motifs_part_1_trimmed.meme $ref/transcript_regions/introns_transcript_regions.fa


