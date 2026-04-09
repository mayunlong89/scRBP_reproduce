#!/bin/bash

######---------------------------------------------------------------------
######---------------------------------------------------------------------
#---AD
######---------------------------------------------------------------------
######---------------------------------------------------------------------

#paths
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
ras_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/01_AD
rgs_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/01_AD
trs_output=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/01_AD/02_assess_lambda_factor
cell2celltype=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD

#lambdas and regions
lambdas=(1 0.05 0.1 0.3 0.5 0.8 2 3 5 10)
regions=(3UTR 5UTR CDS Introns)

echo "Running scRBP_trs_v14.py" >&2

# Loop over lambdas and regions
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "starting ${region}; lambda=${lambda}" >&2
    
    srun --pty --time=10-00:00:00 --cpus-per-task=16 --mem=120G \
      --export=ALL,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1,NUMEXPR_MAX_THREADS=1 \
    python $Tool/scRBP_trs_v14.py \
      --mode ct \
      --ras $ras_data/10X_PBMC_subset_15K_ras_sc_${region}_REAL_PLUS_NULLS.symbol.loom \
      --rgs-csv $rgs_data/AD_processed_magma_results_15Kcells_${region}_ct_entrez_fromMatrix.gsa_RGS.csv \
      --celltypes-csv $cell2celltype/PD_anno_RNAonly_subset_morethan500cell_15Kcells_cell_to_celltype.csv \
      --out-prefix  $trs_output/AD_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda} \
      --lambda-penalty ${lambda} \
      --rgs-score mlog10p \
      --min_cells_pert_ct 25 \
      --q-cap-ras 0.99 \
      --q-cap-rgs 0.99 \
      --do-fdr 1

  done
done

cell_types=("ODC" "N")

################## getEnergyScore.R
echo "Running getEnergyScore.R" >&2

#loop over lambdas, regions, and cell types
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    for cell_type in "${cell_types[@]}"; do
      
      #replace spaces with underscores for naming convention
      cell_type_clean=$(echo "$cell_type" | tr ' ' '_')
      
      echo "Region: ${region}, Lambda: ${lambda}, Cell type: ${cell_type}" >&2
      
      Rscript $Tool/getEnergyScore.R \
        --input "$trs_output/AD_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv" \
        --score_col z_TRS \
        --celltype_col cell_type \
        --regulon_col RBP \
        --target "${cell_type}" \
        -o "$trs_output/${cell_type_clean}_vs_others_summary_${region}_${lambda}"
            
    done
  done
done

################## getSigRegulonCount.R
echo "Running getSigRegulonCount.R" >&2

for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "Region: ${region}, Lambda: ${lambda}" >&2
    
    Rscript $Tool/getSigRegulonCount.R \
      -i $trs_output/AD_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv \
      -o "$trs_output/q01_${region}_lambda_${lambda}" \
      -a 0.05
        
  done
done

######---------------------------------------------------------------------
######---------------------------------------------------------------------
#---PD
######---------------------------------------------------------------------
######---------------------------------------------------------------------
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
ras_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/02_PD
rgs_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/02_PD
trs_output=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/02_PD/02_assess_lambda_factor
cell2celltype=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD

#lambdas and regions
lambdas=(1 0.05 0.1 0.3 0.5 0.8 2 3 5 10)
regions=(3UTR 5UTR CDS Introns)

echo "Running scRBP_trs_v14.py" >&2

# Loop over lambdas and regions
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "starting ${region}; lambda=${lambda}" >&2
    
    srun --pty --time=10-00:00:00 --cpus-per-task=16 --mem=120G \
      --export=ALL,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1,NUMEXPR_MAX_THREADS=1 \
    python $Tool/scRBP_trs_v14.py \
      --mode ct \
      --ras $ras_data/10X_PBMC_subset_15K_ras_sc_${region}_REAL_PLUS_NULLS.symbol.loom \
      --rgs-csv $rgs_data/PD_processed_magma_results_15Kcells_${region}_ct_entrez_fromMatrix.gsa_RGS.csv \
      --celltypes-csv $cell2celltype/PD_anno_RNAonly_subset_morethan500cell_15Kcells_cell_to_celltype.csv \
      --out-prefix  $trs_output/PD_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda} \
      --lambda-penalty ${lambda} \
      --rgs-score mlog10p \
      --min_cells_pert_ct 25 \
      --q-cap-ras 0.99 \
      --q-cap-rgs 0.99 \
      --do-fdr 1

  done
done

cell_types=("ODC" "N")

################## getEnergyScore.R
echo "Running getEnergyScore.R" >&2

#loop over lambdas, regions, and cell types
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    for cell_type in "${cell_types[@]}"; do
      
      #replace spaces with underscores for naming convention
      cell_type_clean=$(echo "$cell_type" | tr ' ' '_')
      
      echo "Region: ${region}, Lambda: ${lambda}, Cell type: ${cell_type}" >&2
      
      Rscript $Tool/getEnergyScore.R \
        --input "$trs_output/PD_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv" \
        --score_col z_TRS \
        --celltype_col cell_type \
        --regulon_col RBP \
        --target "${cell_type}" \
        -o "$trs_output/${cell_type_clean}_vs_others_summary_${region}_${lambda}"
            
    done
  done
done

################## getSigRegulonCount.R
echo "Running getSigRegulonCount.R" >&2

for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "Region: ${region}, Lambda: ${lambda}" >&2
    
    Rscript $Tool/getSigRegulonCount.R \
      -i $trs_output/PD_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv \
      -o "$trs_output/q01_${region}_lambda_${lambda}" \
      -a 0.05
        
  done
done

######---------------------------------------------------------------------
######---------------------------------------------------------------------
#---EAA
######---------------------------------------------------------------------
######---------------------------------------------------------------------
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
ras_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/03_epigenetic_age_acceleration_EAA
rgs_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/03_epigenetic_age_acceleration_EAA
trs_output=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/03_epigenetic_age_acceleration_EAA/02_assess_lambda_factor
cell2celltype=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD

#lambdas and regions
lambdas=(1 0.05 0.1 0.3 0.5 0.8 2 3 5 10)
regions=(3UTR 5UTR CDS Introns)

echo "Running scRBP_trs_v14.py" >&2

# Loop over lambdas and regions
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "starting ${region}; lambda=${lambda}" >&2
    
    srun --pty --time=10-00:00:00 --cpus-per-task=16 --mem=120G \
      --export=ALL,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1,NUMEXPR_MAX_THREADS=1 \
    python $Tool/scRBP_trs_v14.py \
      --mode ct \
      --ras $ras_data/10X_PBMC_subset_15K_ras_sc_${region}_REAL_PLUS_NULLS.symbol.loom \
      --rgs-csv $rgs_data/EAA_processed_magma_results_15Kcells_${region}_ct_entrez_fromMatrix.gsa_RGS.csv \
      --celltypes-csv $cell2celltype/PD_anno_RNAonly_subset_morethan500cell_15Kcells_cell_to_celltype.csv \
      --out-prefix  $trs_output/EAA_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda} \
      --lambda-penalty ${lambda} \
      --rgs-score mlog10p \
      --min_cells_pert_ct 25 \
      --q-cap-ras 0.99 \
      --q-cap-rgs 0.99 \
      --do-fdr 1

  done
done

cell_types=("ODC" "N")

################## getEnergyScore.R
echo "Running getEnergyScore.R" >&2

#loop over lambdas, regions, and cell types
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    for cell_type in "${cell_types[@]}"; do
      
      #replace spaces with underscores for naming convention
      cell_type_clean=$(echo "$cell_type" | tr ' ' '_')
      
      echo "Region: ${region}, Lambda: ${lambda}, Cell type: ${cell_type}" >&2
      
      Rscript $Tool/getEnergyScore.R \
        --input "$trs_output/EAA_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv" \
        --score_col z_TRS \
        --celltype_col cell_type \
        --regulon_col RBP \
        --target "${cell_type}" \
        -o "$trs_output/${cell_type_clean}_vs_others_summary_${region}_${lambda}"
            
    done
  done
done

################## getSigRegulonCount.R
echo "Running getSigRegulonCount.R" >&2

for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "Region: ${region}, Lambda: ${lambda}" >&2
    
    Rscript $Tool/getSigRegulonCount.R \
      -i $trs_output/EAA_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv \
      -o "$trs_output/q01_${region}_lambda_${lambda}" \
      -a 0.05
        
  done
done

######---------------------------------------------------------------------
######---------------------------------------------------------------------
#---90thEL
######---------------------------------------------------------------------
######---------------------------------------------------------------------
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
ras_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/04_extreme_longevity_90thEL
rgs_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/04_extreme_longevity_90thEL
trs_output=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/04_extreme_longevity_90thEL/02_assess_lambda_factor
cell2celltype=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD

#lambdas and regions
lambdas=(1 0.05 0.1 0.3 0.5 0.8 2 3 5 10)
regions=(3UTR 5UTR CDS Introns)

echo "Running scRBP_trs_v14.py" >&2

# Loop over lambdas and regions
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "starting ${region}; lambda=${lambda}" >&2
    
    srun --pty --time=10-00:00:00 --cpus-per-task=16 --mem=120G \
      --export=ALL,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1,NUMEXPR_MAX_THREADS=1 \
    python $Tool/scRBP_trs_v14.py \
      --mode ct \
      --ras $ras_data/10X_PBMC_subset_15K_ras_sc_${region}_REAL_PLUS_NULLS.symbol.loom \
      --rgs-csv $rgs_data/90thEL_processed_magma_results_15Kcells_${region}_ct_entrez_fromMatrix.gsa_RGS.csv \
      --celltypes-csv $cell2celltype/PD_anno_RNAonly_subset_morethan500cell_15Kcells_cell_to_celltype.csv \
      --out-prefix  $trs_output/90thEL_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda} \
      --lambda-penalty ${lambda} \
      --rgs-score mlog10p \
      --min_cells_pert_ct 25 \
      --q-cap-ras 0.99 \
      --q-cap-rgs 0.99 \
      --do-fdr 1

  done
done

cell_types=("ODC" "N")

################## getEnergyScore.R
echo "Running getEnergyScore.R" >&2

#loop over lambdas, regions, and cell types
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    for cell_type in "${cell_types[@]}"; do
      
      #replace spaces with underscores for naming convention
      cell_type_clean=$(echo "$cell_type" | tr ' ' '_')
      
      echo "Region: ${region}, Lambda: ${lambda}, Cell type: ${cell_type}" >&2
      
      Rscript $Tool/getEnergyScore.R \
        --input "$trs_output/90thEL_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv" \
        --score_col z_TRS \
        --celltype_col cell_type \
        --regulon_col RBP \
        --target "${cell_type}" \
        -o "$trs_output/${cell_type_clean}_vs_others_summary_${region}_${lambda}"
            
    done
  done
done

################## getSigRegulonCount.R
echo "Running getSigRegulonCount.R" >&2

for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "Region: ${region}, Lambda: ${lambda}" >&2
    
    Rscript $Tool/getSigRegulonCount.R \
      -i $trs_output/90thEL_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv \
      -o "$trs_output/q01_${region}_lambda_${lambda}" \
      -a 0.05
        
  done
done


######---------------------------------------------------------------------
######---------------------------------------------------------------------
#---99thEL
######---------------------------------------------------------------------
######---------------------------------------------------------------------
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
ras_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/05_extreme_longevity_99thEL
rgs_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/05_extreme_longevity_99thEL
trs_output=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/05_extreme_longevity_99thEL/02_assess_lambda_factor
cell2celltype=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD

#lambdas and regions
lambdas=(1 0.05 0.1 0.3 0.5 0.8 2 3 5 10)
regions=(3UTR 5UTR CDS Introns)

echo "Running scRBP_trs_v14.py" >&2

# Loop over lambdas and regions
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "starting ${region}; lambda=${lambda}" >&2
    
    srun --pty --time=10-00:00:00 --cpus-per-task=16 --mem=120G \
      --export=ALL,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1,NUMEXPR_MAX_THREADS=1 \
    python $Tool/scRBP_trs_v14.py \
      --mode ct \
      --ras $ras_data/10X_PBMC_subset_15K_ras_sc_${region}_REAL_PLUS_NULLS.symbol.loom \
      --rgs-csv $rgs_data/99thEL_processed_magma_results_15Kcells_${region}_ct_entrez_fromMatrix.gsa_RGS.csv \
      --celltypes-csv $cell2celltype/PD_anno_RNAonly_subset_morethan500cell_15Kcells_cell_to_celltype.csv \
      --out-prefix  $trs_output/99thEL_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda} \
      --lambda-penalty ${lambda} \
      --rgs-score mlog10p \
      --min_cells_pert_ct 25 \
      --q-cap-ras 0.99 \
      --q-cap-rgs 0.99 \
      --do-fdr 1

  done
done

cell_types=("ODC" "N")

################## getEnergyScore.R
echo "Running getEnergyScore.R" >&2

#loop over lambdas, regions, and cell types
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    for cell_type in "${cell_types[@]}"; do
      
      #replace spaces with underscores for naming convention
      cell_type_clean=$(echo "$cell_type" | tr ' ' '_')
      
      echo "Region: ${region}, Lambda: ${lambda}, Cell type: ${cell_type}" >&2
      
      Rscript $Tool/getEnergyScore.R \
        --input "$trs_output/99thEL_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv" \
        --score_col z_TRS \
        --celltype_col cell_type \
        --regulon_col RBP \
        --target "${cell_type}" \
        -o "$trs_output/${cell_type_clean}_vs_others_summary_${region}_${lambda}"
            
    done
  done
done

################## getSigRegulonCount.R
echo "Running getSigRegulonCount.R" >&2

for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "Region: ${region}, Lambda: ${lambda}" >&2
    
    Rscript $Tool/getSigRegulonCount.R \
      -i $trs_output/99thEL_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv \
      -o "$trs_output/q01_${region}_lambda_${lambda}" \
      -a 0.05
        
  done
done

######---------------------------------------------------------------------
######---------------------------------------------------------------------
#---HS
######---------------------------------------------------------------------
######---------------------------------------------------------------------
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
ras_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/06_healthspan_HS
rgs_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/06_healthspan_HS
trs_output=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/06_healthspan_HS/02_assess_lambda_factor
cell2celltype=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD

#lambdas and regions
lambdas=(1 0.05 0.1 0.3 0.5 0.8 2 3 5 10)
regions=(3UTR 5UTR CDS Introns)

echo "Running scRBP_trs_v14.py" >&2

# Loop over lambdas and regions
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "starting ${region}; lambda=${lambda}" >&2
    
    srun --pty --time=10-00:00:00 --cpus-per-task=16 --mem=120G \
      --export=ALL,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1,NUMEXPR_MAX_THREADS=1 \
    python $Tool/scRBP_trs_v14.py \
      --mode ct \
      --ras $ras_data/10X_PBMC_subset_15K_ras_sc_${region}_REAL_PLUS_NULLS.symbol.loom \
      --rgs-csv $rgs_data/HS_processed_magma_results_15Kcells_${region}_ct_entrez_fromMatrix.gsa_RGS.csv \
      --celltypes-csv $cell2celltype/PD_anno_RNAonly_subset_morethan500cell_15Kcells_cell_to_celltype.csv \
      --out-prefix  $trs_output/HS_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda} \
      --lambda-penalty ${lambda} \
      --rgs-score mlog10p \
      --min_cells_pert_ct 25 \
      --q-cap-ras 0.99 \
      --q-cap-rgs 0.99 \
      --do-fdr 1

  done
done

cell_types=("ODC" "N")

################## getEnergyScore.R
echo "Running getEnergyScore.R" >&2

#loop over lambdas, regions, and cell types
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    for cell_type in "${cell_types[@]}"; do
      
      #replace spaces with underscores for naming convention
      cell_type_clean=$(echo "$cell_type" | tr ' ' '_')
      
      echo "Region: ${region}, Lambda: ${lambda}, Cell type: ${cell_type}" >&2
      
      Rscript $Tool/getEnergyScore.R \
        --input "$trs_output/HS_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv" \
        --score_col z_TRS \
        --celltype_col cell_type \
        --regulon_col RBP \
        --target "${cell_type}" \
        -o "$trs_output/${cell_type_clean}_vs_others_summary_${region}_${lambda}"
            
    done
  done
done

################## getSigRegulonCount.R
echo "Running getSigRegulonCount.R" >&2

for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "Region: ${region}, Lambda: ${lambda}" >&2
    
    Rscript $Tool/getSigRegulonCount.R \
      -i $trs_output/HS_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv \
      -o "$trs_output/q01_${region}_lambda_${lambda}" \
      -a 0.05
        
  done
done

######---------------------------------------------------------------------
######---------------------------------------------------------------------
#---FI
######---------------------------------------------------------------------
######---------------------------------------------------------------------
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
ras_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/07_frailty_index_FI
rgs_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/07_frailty_index_FI
trs_output=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/07_frailty_index_FI/02_assess_lambda_factor
cell2celltype=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD

#lambdas and regions
lambdas=(1 0.05 0.1 0.3 0.5 0.8 2 3 5 10)
regions=(3UTR 5UTR CDS Introns)

echo "Running scRBP_trs_v14.py" >&2

# Loop over lambdas and regions
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "starting ${region}; lambda=${lambda}" >&2
    
    srun --pty --time=10-00:00:00 --cpus-per-task=16 --mem=120G \
      --export=ALL,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1,NUMEXPR_MAX_THREADS=1 \
    python $Tool/scRBP_trs_v14.py \
      --mode ct \
      --ras $ras_data/10X_PBMC_subset_15K_ras_sc_${region}_REAL_PLUS_NULLS.symbol.loom \
      --rgs-csv $rgs_data/FI_processed_magma_results_15Kcells_${region}_ct_entrez_fromMatrix.gsa_RGS.csv \
      --celltypes-csv $cell2celltype/PD_anno_RNAonly_subset_morethan500cell_15Kcells_cell_to_celltype.csv \
      --out-prefix  $trs_output/FI_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda} \
      --lambda-penalty ${lambda} \
      --rgs-score mlog10p \
      --min_cells_pert_ct 25 \
      --q-cap-ras 0.99 \
      --q-cap-rgs 0.99 \
      --do-fdr 1

  done
done

cell_types=("ODC" "N")

################## getEnergyScore.R
echo "Running getEnergyScore.R" >&2

#loop over lambdas, regions, and cell types
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    for cell_type in "${cell_types[@]}"; do
      
      #replace spaces with underscores for naming convention
      cell_type_clean=$(echo "$cell_type" | tr ' ' '_')
      
      echo "Region: ${region}, Lambda: ${lambda}, Cell type: ${cell_type}" >&2
      
      Rscript $Tool/getEnergyScore.R \
        --input "$trs_output/FI_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv" \
        --score_col z_TRS \
        --celltype_col cell_type \
        --regulon_col RBP \
        --target "${cell_type}" \
        -o "$trs_output/${cell_type_clean}_vs_others_summary_${region}_${lambda}"
            
    done
  done
done

################## getSigRegulonCount.R
echo "Running getSigRegulonCount.R" >&2

for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "Region: ${region}, Lambda: ${lambda}" >&2
    
    Rscript $Tool/getSigRegulonCount.R \
      -i $trs_output/FI_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv \
      -o "$trs_output/q01_${region}_lambda_${lambda}" \
      -a 0.05
        
  done
done

######---------------------------------------------------------------------
######---------------------------------------------------------------------
#---PL
######---------------------------------------------------------------------
######---------------------------------------------------------------------
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
ras_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/08_parent_lifespans_PL
rgs_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/08_parent_lifespans_PL
trs_output=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/08_parent_lifespans_PL/02_assess_lambda_factor
cell2celltype=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD

#lambdas and regions
lambdas=(1 0.05 0.1 0.3 0.5 0.8 2 3 5 10)
regions=(3UTR 5UTR CDS Introns)

echo "Running scRBP_trs_v14.py" >&2

# Loop over lambdas and regions
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "starting ${region}; lambda=${lambda}" >&2
    
    srun --pty --time=10-00:00:00 --cpus-per-task=16 --mem=120G \
      --export=ALL,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1,NUMEXPR_MAX_THREADS=1 \
    python $Tool/scRBP_trs_v14.py \
      --mode ct \
      --ras $ras_data/10X_PBMC_subset_15K_ras_sc_${region}_REAL_PLUS_NULLS.symbol.loom \
      --rgs-csv $rgs_data/PL_processed_magma_results_15Kcells_${region}_ct_entrez_fromMatrix.gsa_RGS.csv \
      --celltypes-csv $cell2celltype/PD_anno_RNAonly_subset_morethan500cell_15Kcells_cell_to_celltype.csv \
      --out-prefix  $trs_output/PL_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda} \
      --lambda-penalty ${lambda} \
      --rgs-score mlog10p \
      --min_cells_pert_ct 25 \
      --q-cap-ras 0.99 \
      --q-cap-rgs 0.99 \
      --do-fdr 1

  done
done

cell_types=("ODC" "N")

################## getEnergyScore.R
echo "Running getEnergyScore.R" >&2

#loop over lambdas, regions, and cell types
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    for cell_type in "${cell_types[@]}"; do
      
      #replace spaces with underscores for naming convention
      cell_type_clean=$(echo "$cell_type" | tr ' ' '_')
      
      echo "Region: ${region}, Lambda: ${lambda}, Cell type: ${cell_type}" >&2
      
      Rscript $Tool/getEnergyScore.R \
        --input "$trs_output/PL_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv" \
        --score_col z_TRS \
        --celltype_col cell_type \
        --regulon_col RBP \
        --target "${cell_type}" \
        -o "$trs_output/${cell_type_clean}_vs_others_summary_${region}_${lambda}"
            
    done
  done
done

################## getSigRegulonCount.R
echo "Running getSigRegulonCount.R" >&2

for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "Region: ${region}, Lambda: ${lambda}" >&2
    
    Rscript $Tool/getSigRegulonCount.R \
      -i $trs_output/PL_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv \
      -o "$trs_output/q01_${region}_lambda_${lambda}" \
      -a 0.05
        
  done
done

######---------------------------------------------------------------------
######---------------------------------------------------------------------
#---mvAge
######---------------------------------------------------------------------
######---------------------------------------------------------------------
Tool=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/10blood_cell_traits_GWAS/10_lymp_count
ras_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/09_mvAge
rgs_data=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/09_mvAge
trs_output=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD/z_GRNBoost2_PD_scRBP_50times/z_scRBP_results_for_9agingtraits_based_on_15Kcells/09_mvAge/02_assess_lambda_factor
cell2celltype=/mnt/isilon/gandal_lab/mayl/02_single-cell_data_ctDRTF/03_brain_PD_AD

#lambdas and regions
lambdas=(1 0.05 0.1 0.3 0.5 0.8 2 3 5 10)
regions=(3UTR 5UTR CDS Introns)

echo "Running scRBP_trs_v14.py" >&2

# Loop over lambdas and regions
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "starting ${region}; lambda=${lambda}" >&2
    
    srun --pty --time=10-00:00:00 --cpus-per-task=16 --mem=120G \
      --export=ALL,OMP_NUM_THREADS=1,MKL_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1,NUMEXPR_MAX_THREADS=1 \
    python $Tool/scRBP_trs_v14.py \
      --mode ct \
      --ras $ras_data/10X_PBMC_subset_15K_ras_sc_${region}_REAL_PLUS_NULLS.symbol.loom \
      --rgs-csv $rgs_data/mvAge_processed_magma_results_15Kcells_${region}_ct_entrez_fromMatrix.gsa_RGS.csv \
      --celltypes-csv $cell2celltype/PD_anno_RNAonly_subset_morethan500cell_15Kcells_cell_to_celltype.csv \
      --out-prefix  $trs_output/mvAge_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda} \
      --lambda-penalty ${lambda} \
      --rgs-score mlog10p \
      --min_cells_pert_ct 25 \
      --q-cap-ras 0.99 \
      --q-cap-rgs 0.99 \
      --do-fdr 1

  done
done

cell_types=("ODC" "N")

################## getEnergyScore.R
echo "Running getEnergyScore.R" >&2

#loop over lambdas, regions, and cell types
for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    for cell_type in "${cell_types[@]}"; do
      
      #replace spaces with underscores for naming convention
      cell_type_clean=$(echo "$cell_type" | tr ' ' '_')
      
      echo "Region: ${region}, Lambda: ${lambda}, Cell type: ${cell_type}" >&2
      
      Rscript $Tool/getEnergyScore.R \
        --input "$trs_output/mvAge_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv" \
        --score_col z_TRS \
        --celltype_col cell_type \
        --regulon_col RBP \
        --target "${cell_type}" \
        -o "$trs_output/${cell_type_clean}_vs_others_summary_${region}_${lambda}"
            
    done
  done
done

################## getSigRegulonCount.R
echo "Running getSigRegulonCount.R" >&2

for lambda in "${lambdas[@]}"; do
  for region in "${regions[@]}"; do
    
    echo "Region: ${region}, Lambda: ${lambda}" >&2
    
    Rscript $Tool/getSigRegulonCount.R \
      -i $trs_output/mvAge_processed_magma_results_15Kcells_${region}_REAL_PLUS_NULLS_lambda_${lambda}.ct.TRS_long.csv \
      -o "$trs_output/q01_${region}_lambda_${lambda}" \
      -a 0.05
        
  done
done
