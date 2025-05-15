#2. based on background model to run FIMO
ref=/mnt/isilon/gandal_lab/mayl/reference/
bgfile=/mnt/isilon/gandal_lab/mayl/reference/

srun --mem=80G --time=4-99:00:00 --pty fimo --verbosity 2 --max-stored-scores 1000000 --thresh 0.001 --bgfile $bgfile/background.meme --oc fimo_output_dir_bg trimmed_refined_bamm2meme_file.meme $ref/gene_sequence_v44_GCRh38.fasta

