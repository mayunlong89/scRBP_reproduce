
# This is script for generating regulon-wise association scores from GWAS using MAGMA.
#### 2025-08-11



## Generate MAGMA-based gene_set results
```bash
# 1) MAGMA codes for generating disease-relevant association scores for TFs and target genes

#DIRECTORY
export MAGMA_DIR=/share/pub/mayl/MAGMA
export DATA=/share/pub/mayl/MAGMA_test
export OUTPUT=/share/pub/mayl/MAGMA_test

#MAGMA annotation:
# By default, a 10 kb window centered on the TSS of a gene is used.

$MAGMA_DIR/magma \
    --snp-loc  $DATA/GWAS_UKBiobank_summary_final.hg19.location  \
    --annotate window=10,10 --gene-loc $MAGMA_DIR/NCBI37.3.gene.loc \
    --out $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_annotation  

# gene-based association analysis:
$MAGMA_DIR/magma \
    --bfile $MAGMA_DIR/1000G_data/g1000_eur \
    --pval $DATA/GWAS_UKBiobank_summary_final.results_Pval \
    N=13239 \
    --gene-annot   $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_annotation.genes.annot  \
    --out $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_Analysis_P

# gene-set analysis
geneset=/home/he/tools/magma/MSigDB/msigdb_v2022.1.Hs_files_to_download_locally/msigdb_v2022.1.Hs_GMTs/msigdb.v2022.1.Hs.entrez.gmt
$MAGMA_DIR/magma \
    --gene-results GWAS_UKBiobank_summary_final.hg19_SNP_Gene_Analysis_P.raw \
    --set-annot ${geneset} \
    --out GWAS_UKBiobank_summary_final.hg19_geneset


```
