#2025-12-30 (Only SNP rsIDs can be annotated; indels cannot)

#A. First parse your SNP list into chr / pos / ref / alt
library(data.table)
library(stringr)
library(GenomicRanges)

# Read in a single-column file
x <- fread("snps.txt", header = FALSE)$V1
x <- str_trim(x)

# Allow spaces or underscores in the middle:
# replace all non-alphanumeric characters with underscores, then split
x2 <- gsub("[^A-Za-z0-9]+", "_", x)
parts <- tstrsplit(x2, "_", fixed = TRUE)

df <- data.frame(
  chr = parts[[1]],
  pos = as.integer(parts[[2]]),
  ref = parts[[3]],
  alt = parts[[4]],
  raw = x,
  stringsAsFactors = FALSE
)

# Standardize chr prefix
# (the style will be adjusted again later based on the reference database)
df$chr <- ifelse(grepl("^chr", df$chr, ignore.case = TRUE), df$chr, paste0("chr", df$chr))

# Build GRanges:
# use width = ref length for both SNPs and indels (more robust)
gr_all <- GRanges(
  seqnames = df$chr,
  ranges   = IRanges(start = df$pos, width = nchar(df$ref)),
  strand   = "*"
)
mcols(gr_all)$ref <- df$ref
mcols(gr_all)$alt <- df$alt


#B. First annotate “SNPs (single-base variants)” with rsIDs using SNPlocs
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(GenomeInfoDb)

snp_db <- SNPlocs.Hsapiens.dbSNP155.GRCh38

# Key step: unify seqlevels style (chr20 vs 20)
# This line usually solves the hits = 0 issue directly
seqlevelsStyle(gr_all) <- seqlevelsStyle(snp_db)

# Keep only single-base SNPs
is_snp <- (nchar(df$ref) == 1 & nchar(df$alt) == 1)
gr_snp <- gr_all[is_snp]

# Use the Bioconductor-recommended interface: snpsByOverlaps
library(GenomicFeatures)
snps_hit <- snpsByOverlaps(snp_db, gr_snp)

snps_hit
# The returned object contains RefSNP_id
# (the numeric part after "rs")


# Merge rsIDs back into the original df
# Match by position; for stricter matching, ref/alt can also be compared
# snps_hit positions are single-base, so start = pos
hit_df <- data.frame(
  chr = as.character(seqnames(snps_hit)),
  pos = start(snps_hit),
  rs  = paste0("", snps_hit$RefSNP_id),
  stringsAsFactors = FALSE
)

# Return to the original df
# Note that df chr may be in chr20 format,
# but we already aligned the style with the database once above
# A safer approach: create another chr column in the same style for merging
df2 <- df
df2$chr2 <- df2$chr
# Convert df2$chr2 to the same style as snp_db
tmp_gr <- GRanges(df2$chr2, IRanges(df2$pos, width = 1))
seqlevelsStyle(tmp_gr) <- seqlevelsStyle(snp_db)
df2$chr2 <- as.character(seqnames(tmp_gr))

df2$rsID <- NA_character_
idx <- match(paste(df2$chr2, df2$pos), paste(hit_df$chr, hit_df$pos))
df2$rsID[is_snp] <- hit_df$rs[idx[is_snp]]
