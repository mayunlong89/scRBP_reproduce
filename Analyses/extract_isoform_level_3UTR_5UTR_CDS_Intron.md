
`2025-05-15`
`2024-07-21`

## 1. extract isoform_level 3UTR, 5UTR, CDS, intronic regions

### 1.1 Gene_coordinate_genomicsRanges.R 

```R
# Set working directory
setwd("/Users/mayunlong/Desktop/UPENN/00000------------UPenn-Projects/02-Aimed projects/05-devBrain-RBP-methods/RBP_data/")

install.packages("RMariaDB")

library(GenomicFeatures)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")

# Extract 3'UTR regions ---------------------------------------------
# Build TxDb object from UCSC
txdb <- makeTxDbFromUCSC(genome="hg38", tablename='knownGene')

# Extract 3'UTR information
threeUTRs <- threeUTRsByTranscript(txdb, use.names=TRUE)

# Convert 3'UTR to GRanges object
threeUTRs_gr <- unlist(threeUTRs, use.names=TRUE)

# Create BED file
# Add row names as a column
threeUTRs_gr$name <- names(threeUTRs_gr)

# Assign a unique ID to each row
row_ids <- paste0("ID_", seq_along(threeUTRs_gr))

# Replace row names with unique IDs
names(threeUTRs_gr) <- row_ids

# Convert GRanges object to data frame
threeUTRs_df <- as.data.frame(threeUTRs_gr)

# Format as BED file
bed_data <- data.frame(
  chrom = threeUTRs_df$seqnames,
  start = threeUTRs_df$start - 1,  # 0-based start position for BED files
  end = threeUTRs_df$end,
  name = threeUTRs_df$name,
  score = ".",
  strand = threeUTRs_df$strand
)

# Remove duplicate rows, keep the first occurrence
bed_data_unique <- bed_data[!duplicated(bed_data$name),]

# Save BED file
write.table(bed_data_unique, file="3UTRs_hg38.bed", quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

# Extract 5'UTR regions ---------------------------------------------
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(GenomicRanges)

# Build TxDb object
txdb <- makeTxDbFromUCSC(genome="hg38", tablename="knownGene")

# Extract 5'UTR information
fiveUTRs <- fiveUTRsByTranscript(txdb, use.names=TRUE)

# Convert 5'UTR to GRanges object
fiveUTRs_gr <- unlist(fiveUTRs, use.names=TRUE)

# Add row names as a column
fiveUTRs_gr$name <- names(fiveUTRs_gr)

# Assign a unique ID to each row
row_ids <- paste0("ID_", seq_along(fiveUTRs_gr))

# Replace row names with unique IDs
names(fiveUTRs_gr) <- row_ids

# Convert GRanges object to data frame
fiveUTRs_df <- as.data.frame(fiveUTRs_gr)

# Format as BED file
bed_data <- data.frame(
  chrom = fiveUTRs_df$seqnames,
  start = fiveUTRs_df$start - 1,  # 0-based start position for BED files
  end = fiveUTRs_df$end,
  name = fiveUTRs_df$name,
  score = ".",
  strand = fiveUTRs_df$strand
)

# Remove duplicate rows, keep the first occurrence
bed_data_unique <- bed_data[!duplicated(bed_data$name),]

# Save BED file
write.table(bed_data, file="5UTRs_hg38.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# Extract CDS genomic regions -----------------------------------------
txdb <- makeTxDbFromUCSC(genome="hg38", tablename="knownGene")

# Extract CDS regions
cds <- cdsBy(txdb, by="tx", use.names=TRUE)

# Convert CDS to GRanges object
cds_gr <- unlist(cds, use.names=TRUE)

# Add row names as a column
cds_gr$name <- names(cds_gr)

# Assign a unique ID to each row
row_ids <- paste0("ID_", seq_along(cds_gr))

# Replace row names with unique IDs
names(cds_gr) <- row_ids

# Convert GRanges object to data frame
cds_df <- as.data.frame(cds_gr)

# Format as BED file
bed_data <- data.frame(
  chrom = cds_df$seqnames,
  start = cds_df$start - 1,  # 0-based start position for BED files
  end = cds_df$end,
  name = cds_df$name,
  score = ".",
  strand = cds_df$strand
)

# Remove duplicate rows, keep the first occurrence
bed_data_unique <- bed_data[!duplicated(bed_data$name),]

# Save BED file
write.table(bed_data, file="CDS_coordinate_hg38.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# Extract intron regions ---------------------------------------------
txdb <- makeTxDbFromUCSC(genome="hg38", tablename="knownGene")

# Extract intron regions
introns <- intronsByTranscript(txdb, use.names=TRUE)

# Convert to GRanges format
introns_gr <- unlist(introns, use.names=TRUE)

# Add row names as a column
introns_gr$name <- names(introns_gr)

# Assign a unique ID to each row
row_ids <- paste0("ID_", seq_along(introns_gr))

# Replace row names with unique IDs
names(introns_gr) <- row_ids

# Convert GRanges object to data frame
introns_df <- as.data.frame(introns_gr)

# Format as BED file
bed_data <- data.frame(
  chrom = introns_df$seqnames,
  start = introns_df$start - 1,  # 0-based start position for BED files
  end = introns_df$end,
  name = introns_df$name,
  score = ".",
  strand = introns_df$strand
)

# Remove duplicate rows, keep the first occurrence
bed_data_unique <- bed_data[!duplicated(bed_data$name),]

# Save BED file
write.table(bed_data_unique, file="introns_hg38.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


```

### 1.2 bedtools --getfasta extract sequences of 3UTR, 5UTR, CDS, Intronic regions convert into .fasta files
```bash

ref=/mnt/isilon/gandal_lab/mayl/reference/

# Extract 3'UTR sequences using BED file
bedtools getfasta \
-fi $ref/gene_sequence_v44_GCRh38.fasta \
-bed 3UTR.best_hit.hg38.bed \
-fo 3UTR_motif_binding_sites.fasta \
-name


```



