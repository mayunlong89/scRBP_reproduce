#!/usr/bin/env Rscript

## ============================================================
## Gene_coordinate_longread_isoform_scRBP.R
##
## Purpose:
##   Extract transcript-domain genomic regions from long-read
##   isoform annotations for scRBP motif scanning.
##
## Input:
##   Long-read isoform GTF file
##
## Main outputs per domain:
##   1. *.scan_regions.bed
##      Coordinate-deduplicated BED file for motif scanning.
##
##   2. *.region_to_transcript.tsv
##      Mapping table:
##      region_id -> transcript_id -> gene_id/gene_name/domain.
##
##   3. *.transcript_segments.bed
##      Transcript-segment-level BED preserving isoform-specific fragments.
##
##   4. *.transcript_segments.metadata.tsv
##      Metadata for transcript-specific fragments.
##
## Notes:
##   - BED start is 0-based.
##   - GTF/GenomicRanges start/end are 1-based.
##   - Do NOT deduplicate by transcript_id only.
##   - Use bedtools getfasta -s downstream for strand-aware sequence extraction.
## ============================================================


## -----------------------------
## 1. Load libraries
## -----------------------------

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
})


## -----------------------------
## 2. User parameters
## -----------------------------

## Long-read isoform annotation GTF
GTF_FILE <- "/path/to/your_long_read_isoforms.gtf"

## Output directory
OUT_DIR <- "/path/to/output/longread_scRBP_regions"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Genome label used only for output file names
GENOME_LABEL <- "hg38"

## Keep only standard chromosomes?
## For hg38, this keeps chr1-22/chrX/chrY and also supports 1-22/X/Y style.
KEEP_STANDARD_CHROMS <- TRUE

standard_chroms <- c(
  paste0("chr", 1:22), "chrX", "chrY",
  as.character(1:22), "X", "Y"
)

## Drop regions without clear strand?
## For RBP motif scanning, strand-aware extraction requires + or -.
DROP_UNSTRANDED <- TRUE

## Extract exon regions as an extra domain.
## This is useful because many long-read GTF files only contain exon/transcript
## features and do not contain CDS/UTR annotations.
EXTRACT_EXONS <- TRUE


## -----------------------------
## 3. Helper functions
## -----------------------------

first_non_missing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "."]
  if (length(x) == 0) {
    return(NA_character_)
  } else {
    return(x[1])
  }
}


get_chrom_order <- function(chrom) {
  x <- gsub("^chr", "", chrom)
  ref <- c(as.character(1:22), "X", "Y", "M", "MT")
  idx <- match(x, ref)
  idx[is.na(idx)] <- 1000L
  return(idx)
}


get_feature_count <- function(gtf_dt) {
  if (!"type" %in% names(gtf_dt)) {
    return(data.table(feature = NA_character_, n = nrow(gtf_dt)))
  }
  
  gtf_dt[
    ,
    .(n = .N),
    by = .(feature = type)
  ][order(-n)]
}


make_tx2gene_from_gtf <- function(gtf_dt) {
  
  if (!"transcript_id" %in% names(gtf_dt)) {
    stop("The GTF file does not contain transcript_id. Please check the annotation format.")
  }
  
  if (!"gene_id" %in% names(gtf_dt)) {
    gtf_dt[, gene_id := NA_character_]
  }
  
  if (!"gene_name" %in% names(gtf_dt)) {
    gtf_dt[, gene_name := NA_character_]
  }
  
  if (!"type" %in% names(gtf_dt)) {
    gtf_dt[, type := NA_character_]
  }
  
  tx_raw <- gtf_dt[
    !is.na(transcript_id) & transcript_id != "",
    .(
      transcript_id = as.character(transcript_id),
      gene_id = as.character(gene_id),
      gene_name = as.character(gene_name),
      feature_type = as.character(type)
    )
  ]
  
  ## Prefer transcript records, then exon, then other features.
  tx_raw[
    ,
    feature_priority := fifelse(
      feature_type == "transcript", 1L,
      fifelse(feature_type == "exon", 2L, 3L)
    )
  ]
  
  setorder(tx_raw, transcript_id, feature_priority)
  
  tx2gene <- tx_raw[
    ,
    .(
      gene_id = first_non_missing(gene_id),
      gene_name = first_non_missing(gene_name)
    ),
    by = transcript_id
  ]
  
  tx2gene <- unique(tx2gene, by = "transcript_id")
  
  return(tx2gene)
}


extract_one_domain <- function(grl,
                               label,
                               tx2gene,
                               out_dir,
                               genome_label = "hg38",
                               keep_standard_chroms = TRUE,
                               standard_chroms = NULL,
                               drop_unstranded = TRUE) {
  
  message("\n============================================================")
  message("Processing domain: ", label)
  message("============================================================")
  
  ## Remove empty transcripts
  grl <- grl[elementNROWS(grl) > 0]
  
  if (length(grl) == 0) {
    message("[", label, "] empty. Write empty QC and skip.")
    
    qc <- data.table(
      domain = label,
      n_transcript_fragments = 0,
      n_unique_scan_regions = 0,
      n_unique_transcripts = 0,
      n_unique_genes = 0,
      n_missing_gene_id_fragments = 0
    )
    
    return(qc)
  }
  
  n_per_tx <- elementNROWS(grl)
  gr <- unlist(grl, use.names = FALSE)
  transcript_ids <- rep(names(grl), n_per_tx)
  part_index <- sequence(n_per_tx)
  
  if (any(is.na(transcript_ids)) || any(transcript_ids == "")) {
    warning("[", label, "] Some transcript IDs are missing from GRangesList names.")
  }
  
  fragment_dt <- data.table(
    domain = label,
    chrom = as.character(seqnames(gr)),
    start = start(gr) - 1L,       # BED 0-based start
    end = end(gr),                # BED half-open end
    transcript_id = as.character(transcript_ids),
    part_index = part_index,
    strand = as.character(strand(gr)),
    width = width(gr)
  )
  
  ## Keep standard chromosomes if requested
  if (keep_standard_chroms) {
    n_before <- nrow(fragment_dt)
    fragment_dt <- fragment_dt[chrom %in% standard_chroms]
    n_after <- nrow(fragment_dt)
    message("[", label, "] kept standard chromosomes: ", n_after, " / ", n_before)
  }
  
  ## Drop unstranded regions if requested
  if (drop_unstranded) {
    n_before <- nrow(fragment_dt)
    fragment_dt <- fragment_dt[strand %in% c("+", "-")]
    n_after <- nrow(fragment_dt)
    message("[", label, "] kept stranded regions: ", n_after, " / ", n_before)
  }
  
  if (nrow(fragment_dt) == 0) {
    message("[", label, "] no regions left after chromosome/strand filtering.")
    
    qc <- data.table(
      domain = label,
      n_transcript_fragments = 0,
      n_unique_scan_regions = 0,
      n_unique_transcripts = 0,
      n_unique_genes = 0,
      n_missing_gene_id_fragments = 0
    )
    
    return(qc)
  }
  
  ## Remove exact duplicate fragments within the same transcript.
  ## This is safe. It does NOT remove different CDS/intron fragments from the same transcript.
  fragment_dt <- unique(
    fragment_dt,
    by = c("domain", "chrom", "start", "end", "transcript_id", "strand")
  )
  
  ## Add gene_id and gene_name
  fragment_dt <- merge(
    fragment_dt,
    tx2gene,
    by = "transcript_id",
    all.x = TRUE
  )
  
  ## Coordinate-deduplicated scan regions.
  ## Same genomic sequence is scanned once, even if shared by many isoforms.
  scan_dt <- unique(
    fragment_dt[, .(domain, chrom, start, end, strand)],
    by = c("domain", "chrom", "start", "end", "strand")
  )
  
  scan_dt[, chrom_order := get_chrom_order(chrom)]
  setorder(scan_dt, chrom_order, start, end, strand)
  scan_dt[, chrom_order := NULL]
  
  scan_dt[
    ,
    region_id := sprintf(
      "%s_region_%09d",
      label,
      seq_len(.N)
    )
  ]
  
  ## Add region_id back to transcript fragments
  map_dt <- merge(
    fragment_dt,
    scan_dt,
    by = c("domain", "chrom", "start", "end", "strand"),
    all.x = TRUE
  )
  
  setcolorder(
    map_dt,
    c(
      "region_id",
      "domain",
      "transcript_id",
      "gene_id",
      "gene_name",
      "part_index",
      "chrom",
      "start",
      "end",
      "strand",
      "width"
    )
  )
  
  map_dt[, chrom_order := get_chrom_order(chrom)]
  setorder(map_dt, chrom_order, start, end, transcript_id, strand)
  map_dt[, chrom_order := NULL]
  
  ## Transcript-segment-level ID
  map_dt[
    ,
    transcript_segment_id := sprintf(
      "%s_txseg_%09d",
      label,
      seq_len(.N)
    )
  ]
  
  ## BED for motif scanning: one row per unique genomic coordinate
  scan_bed <- scan_dt[
    ,
    .(
      chrom,
      start,
      end,
      name = region_id,
      score = ".",
      strand
    )
  ]
  
  ## BED preserving transcript-specific segments
  txseg_bed <- map_dt[
    ,
    .(
      chrom,
      start,
      end,
      name = transcript_segment_id,
      score = ".",
      strand
    )
  ]
  
  ## Output file names
  prefix <- file.path(
    out_dir,
    paste0(label, "_", genome_label)
  )
  
  scan_bed_file <- paste0(prefix, ".scan_regions.bed")
  map_file <- paste0(prefix, ".region_to_transcript.tsv")
  txseg_bed_file <- paste0(prefix, ".transcript_segments.bed")
  txseg_meta_file <- paste0(prefix, ".transcript_segments.metadata.tsv")
  
  ## Write files
  fwrite(
    scan_bed,
    scan_bed_file,
    sep = "\t",
    col.names = FALSE,
    quote = FALSE
  )
  
  fwrite(
    map_dt,
    map_file,
    sep = "\t",
    col.names = TRUE,
    quote = FALSE
  )
  
  fwrite(
    txseg_bed,
    txseg_bed_file,
    sep = "\t",
    col.names = FALSE,
    quote = FALSE
  )
  
  fwrite(
    map_dt,
    txseg_meta_file,
    sep = "\t",
    col.names = TRUE,
    quote = FALSE
  )
  
  qc <- data.table(
    domain = label,
    n_transcript_fragments = nrow(map_dt),
    n_unique_scan_regions = nrow(scan_bed),
    n_unique_transcripts = uniqueN(map_dt$transcript_id),
    n_unique_genes = uniqueN(na.omit(map_dt$gene_id)),
    n_missing_gene_id_fragments = sum(is.na(map_dt$gene_id) | map_dt$gene_id == ""),
    scan_bed_file = scan_bed_file,
    region_to_transcript_file = map_file,
    transcript_segments_bed_file = txseg_bed_file,
    transcript_segments_metadata_file = txseg_meta_file
  )
  
  print(qc)
  
  return(qc)
}


## -----------------------------
## 4. Import GTF and build tx2gene map
## -----------------------------

message("Input GTF: ", GTF_FILE)
message("Output directory: ", OUT_DIR)

if (!file.exists(GTF_FILE)) {
  stop("GTF file does not exist: ", GTF_FILE)
}

message("\nImporting GTF...")
gtf_gr <- rtracklayer::import(GTF_FILE)
gtf_dt <- as.data.table(as.data.frame(gtf_gr))

message("\nGTF feature counts:")
feature_count <- get_feature_count(gtf_dt)
print(feature_count)

fwrite(
  feature_count,
  file.path(OUT_DIR, "GTF_feature_counts.tsv"),
  sep = "\t"
)

message("\nBuilding transcript-to-gene mapping...")
tx2gene <- make_tx2gene_from_gtf(gtf_dt)

fwrite(
  tx2gene,
  file.path(OUT_DIR, "isoform_gene_map.tsv"),
  sep = "\t",
  col.names = TRUE,
  quote = FALSE
)

message("Number of transcripts in tx2gene map: ", nrow(tx2gene))
message("Number of genes in tx2gene map: ", uniqueN(na.omit(tx2gene$gene_id)))


## -----------------------------
## 5. Build TxDb from long-read GTF
## -----------------------------

message("\nBuilding TxDb from GTF...")
txdb <- makeTxDbFromGFF(
  GTF_FILE,
  format = "gtf"
)


## -----------------------------
## 6. Extract domains
## -----------------------------

qc_list <- list()

## 5'UTR
message("\nExtracting 5'UTR regions...")
fiveUTR_grl <- fiveUTRsByTranscript(
  txdb,
  use.names = TRUE
)

qc_list[["fiveUTR"]] <- extract_one_domain(
  grl = fiveUTR_grl,
  label = "fiveUTR",
  tx2gene = tx2gene,
  out_dir = OUT_DIR,
  genome_label = GENOME_LABEL,
  keep_standard_chroms = KEEP_STANDARD_CHROMS,
  standard_chroms = standard_chroms,
  drop_unstranded = DROP_UNSTRANDED
)

## 3'UTR
message("\nExtracting 3'UTR regions...")
threeUTR_grl <- threeUTRsByTranscript(
  txdb,
  use.names = TRUE
)

qc_list[["threeUTR"]] <- extract_one_domain(
  grl = threeUTR_grl,
  label = "threeUTR",
  tx2gene = tx2gene,
  out_dir = OUT_DIR,
  genome_label = GENOME_LABEL,
  keep_standard_chroms = KEEP_STANDARD_CHROMS,
  standard_chroms = standard_chroms,
  drop_unstranded = DROP_UNSTRANDED
)

## CDS
message("\nExtracting CDS regions...")
CDS_grl <- cdsBy(
  txdb,
  by = "tx",
  use.names = TRUE
)

qc_list[["CDS"]] <- extract_one_domain(
  grl = CDS_grl,
  label = "CDS",
  tx2gene = tx2gene,
  out_dir = OUT_DIR,
  genome_label = GENOME_LABEL,
  keep_standard_chroms = KEEP_STANDARD_CHROMS,
  standard_chroms = standard_chroms,
  drop_unstranded = DROP_UNSTRANDED
)

## Intron
message("\nExtracting intron regions...")
intron_grl <- intronsByTranscript(
  txdb,
  use.names = TRUE
)

qc_list[["intron"]] <- extract_one_domain(
  grl = intron_grl,
  label = "intron",
  tx2gene = tx2gene,
  out_dir = OUT_DIR,
  genome_label = GENOME_LABEL,
  keep_standard_chroms = KEEP_STANDARD_CHROMS,
  standard_chroms = standard_chroms,
  drop_unstranded = DROP_UNSTRANDED
)

## Exon
if (EXTRACT_EXONS) {
  message("\nExtracting exon regions...")
  exon_grl <- exonsBy(
    txdb,
    by = "tx",
    use.names = TRUE
  )
  
  qc_list[["exon"]] <- extract_one_domain(
    grl = exon_grl,
    label = "exon",
    tx2gene = tx2gene,
    out_dir = OUT_DIR,
    genome_label = GENOME_LABEL,
    keep_standard_chroms = KEEP_STANDARD_CHROMS,
    standard_chroms = standard_chroms,
    drop_unstranded = DROP_UNSTRANDED
  )
}


## -----------------------------
## 7. Save QC summary
## -----------------------------

qc_summary <- rbindlist(qc_list, fill = TRUE)

qc_file <- file.path(
  OUT_DIR,
  paste0("longread_scRBP_region_QC_summary_", GENOME_LABEL, ".tsv")
)

fwrite(
  qc_summary,
  qc_file,
  sep = "\t",
  col.names = TRUE,
  quote = FALSE
)

message("\nDone.")
message("QC summary written to: ", qc_file)
message("Output directory: ", OUT_DIR)
