
#2025-10-13

See the code `plot_trs_null.py` in [here](https://github.com/mayunlong89/scRBP_reproduce/blob/main/Analyses/plot_trs_null.py)

```bash
#---6.6K cells  


#DDX21
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_6.5Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_6.5K_100runs_TRS_ct_monocyte_count_3UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp DDX21 \
  --cell-type Monocytes \
  --out DDX21_Monocytes.TRS_null_6.5K_monocyte_count_3UTR.pdf
  
  
 
#DDX21
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_6.5Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_6.5K_100runs_TRS_ct_monocyte_count_5UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp DDX21 \
  --cell-type Monocytes \
  --out DDX21_Monocytes.TRS_null_6.5K_monocyte_count_5UTR.pdf
  
 
 
 
#DDX3X
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_6.5Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_6.5K_100runs_TRS_ct_monocyte_count_3UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp DDX3X \
  --cell-type Monocytes \
  --out DDX3X_Monocytes.TRS_null_6.5K_monocyte_count_3UTR.pdf
  
  
 
#DDX3X
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_6.5Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_6.5K_100runs_TRS_ct_monocyte_count_5UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp DDX3X \
  --cell-type Monocytes \
  --out DDX3X_Monocytes.TRS_null_6.5K_monocyte_count_5UTR.pdf
  



#TRA2A
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_6.5Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_6.5K_100runs_TRS_ct_monocyte_count_5UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp TRA2A \
  --cell-type Monocytes \
  --out TRA2A_Monocytes.TRS_null_6.5K_monocyte_count_5UTR.pdf
  
 
#FAM120A
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_6.5Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_6.5K_100runs_TRS_ct_monocyte_count_5UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp FAM120A \
  --cell-type Monocytes \
  --out FAM120A_Monocytes.TRS_null_6.5K_monocyte_count_5UTR.pdf
  
 
 
 #QKI
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_6.5Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_6.5K_100runs_TRS_ct_monocyte_count_5UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp QKI \
  --cell-type Monocytes \
  --out QKI_Monocytes.TRS_null_6.5K_monocyte_count_5UTR.pdf
  
 

 #YTHDF3
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_6.5Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_6.5K_100runs_TRS_ct_monocyte_count_5UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp YTHDF3 \
  --cell-type Monocytes \
  --out YTHDF3_Monocytes.TRS_null_6.5K_monocyte_count_5UTR.pdf
  
 


#CELF1
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_6.5Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_6.5K_100runs_TRS_ct_monocyte_count_3UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp CELF1 \
  --cell-type Monocytes \
  --out CELF1_Monocytes.TRS_null_6.5K_monocyte_count_3UTR.pdf
  
  
#YBX3
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_6.5Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_6.5K_100runs_TRS_ct_monocyte_count_Introns_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp YBX3 \
  --cell-type Monocytes \
  --out YBX3_Monocytes.TRS_null_6.5K_monocyte_count_Introns.pdf
    
  
#RBM47
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_6.5Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_6.5K_100runs_TRS_ct_monocyte_count_Introns_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp RBM47 \
  --cell-type Monocytes \
  --out RBM47_Monocytes.TRS_null_6.5K_monocyte_count_Introns.pdf
    
    
 
 #AKAP1
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_6.5Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_6.5K_100runs_TRS_ct_monocyte_count_5UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp AKAP1 \
  --cell-type T_cells \
  --out AKAP1_Tcells.TRS_null_6.5K_monocyte_count_5UTR.pdf
    
    
 
 
 
 
#---------8.8K cells

#2025-10-13

  

#DDX21
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_8.7Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_8.7K_100runs_TRS_ct_monocyte_count_3UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp DDX21 \
  --cell-type Monocytes \
  --out DDX21_Monocytes.TRS_null_8.7K_monocyte_count_3UTR.pdf
  
  
 
#DDX21
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_8.7Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_8.7K_100runs_TRS_ct_monocyte_count_5UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp DDX21 \
  --cell-type Monocytes \
  --out DDX21_Monocytes.TRS_null_8.7K_monocyte_count_5UTR.pdf
  
 
 
 
#DDX3X
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_8.7Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_8.7K_100runs_TRS_ct_monocyte_count_3UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp DDX3X \
  --cell-type Monocytes \
  --out DDX3X_Monocytes.TRS_null_8.7K_monocyte_count_3UTR.pdf
  
  
 
#DDX3X
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_8.7Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_8.7K_100runs_TRS_ct_monocyte_count_5UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp DDX3X \
  --cell-type Monocytes \
  --out DDX3X_Monocytes.TRS_null_8.7K_monocyte_count_5UTR.pdf
  



#TRA2A
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_8.7Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_8.7K_100runs_TRS_ct_monocyte_count_5UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp TRA2A \
  --cell-type Monocytes \
  --out TRA2A_Monocytes.TRS_null_8.7K_monocyte_count_5UTR.pdf
  
 
#FAM120A
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_8.7Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_8.7K_100runs_TRS_ct_monocyte_count_5UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp FAM120A \
  --cell-type Monocytes \
  --out FAM120A_Monocytes.TRS_null_8.7K_monocyte_count_5UTR.pdf
  
 
 
 #QKI
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_8.7Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_8.7K_100runs_TRS_ct_monocyte_count_5UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp QKI \
  --cell-type Monocytes \
  --out QKI_Monocytes.TRS_null_8.7K_monocyte_count_5UTR.pdf
  
 

 #YTHDF3
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_8.7Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_8.7K_100runs_TRS_ct_monocyte_count_5UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp YTHDF3 \
  --cell-type Monocytes \
  --out YTHDF3_Monocytes.TRS_null_8.7K_monocyte_count_5UTR.pdf
  
 


#CELF1
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_8.7Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_8.7K_100runs_TRS_ct_monocyte_count_3UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp CELF1 \
  --cell-type Monocytes \
  --out CELF1_Monocytes.TRS_null_8.7K_monocyte_count_3UTR.pdf
  
  
#YBX3
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_8.7Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_8.7K_100runs_TRS_ct_monocyte_count_Introns_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp YBX3 \
  --cell-type Monocytes \
  --out YBX3_Monocytes.TRS_null_8.7K_monocyte_count_Introns.pdf
    
  
#RBM47
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_8.7Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_8.7K_100runs_TRS_ct_monocyte_count_Introns_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp RBM47 \
  --cell-type Monocytes \
  --out RBM47_Monocytes.TRS_null_8.7K_monocyte_count_Introns.pdf
    
    
 
 #AKAP1
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
OUT=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/01_Long_read_single_cell_data/10X_PBMC_discovery/03_8.7Kcell_100runs_results/01_monocyte_count_trs

python $Tool/plot_trs_null.py \
  --npz $OUT/10X_PBMC_subset_8.7K_100runs_TRS_ct_monocyte_count_5UTR_NULLS_distribution_v15.ct.nulls_perRBP.npz \
  --rbp AKAP1 \
  --cell-type T_cells \
  --out AKAP1_Tcells.TRS_null_8.7K_monocyte_count_5UTR.pdf
    
```
 
 
 
 
 
