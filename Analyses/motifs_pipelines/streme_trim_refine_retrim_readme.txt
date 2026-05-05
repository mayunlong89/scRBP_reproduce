###--------------------------------
output_prefix/
  ├── motifs/
  │   ├── file1_motifs_dir/
  │   └── file2_motifs_dir/
  ├── trimmed_motifs/
  │   ├── file1_trimmed_meme_file.meme
  │   └── file2_trimmed_meme_file.meme
  ├── refined_motifs/
  │   ├── file1_refine/
  │   └── file2_refine/
  ├── merged_motifs/
  │   ├── file1_merged_meme2pwm.pwm
  │   └── file2_merged_meme2pwm.pwm
  ├── file1_refine_bamm2meme.meme
  └── file1_trimmed_refined_bamm2meme_file.meme
##-----------------------------------------------------------------------------------------------------

#test example
bash streme_trim_refine_retrim.sh 05_test_combinded_fasta 05_test_motifs

