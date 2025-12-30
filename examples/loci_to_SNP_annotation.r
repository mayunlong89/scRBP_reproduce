#2025-12-30 (只能注释SNP rs, 不能注释indel)

#A. 先把你的 SNP list 解析成 chr / pos / ref / alt
library(data.table)
library(stringr)
library(GenomicRanges)

# 读入一列
x <- fread("snps.txt", header = FALSE)$V1
x <- str_trim(x)

# 允许中间有空格或下划线：把所有非字母数字替换成下划线再 split
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

# 统一 chr 前缀（后面会再根据参考数据库自动调 style）
df$chr <- ifelse(grepl("^chr", df$chr, ignore.case = TRUE), df$chr, paste0("chr", df$chr))

# 建 GRanges：SNP/indel 都用宽度=ref长度（更稳）
gr_all <- GRanges(
  seqnames = df$chr,
  ranges   = IRanges(start = df$pos, width = nchar(df$ref)),
  strand   = "*"
)
mcols(gr_all)$ref <- df$ref
mcols(gr_all)$alt <- df$alt


#B. 先把 “SNP(单碱基)” 用 SNPlocs 注释 rsID
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(GenomeInfoDb)

snp_db <- SNPlocs.Hsapiens.dbSNP155.GRCh38

# 关键：统一 seqlevels 风格（chr20 vs 20）
seqlevelsStyle(gr_all) <- seqlevelsStyle(snp_db)   # 这句通常能直接解决 hits=0

# 只取单碱基 SNP
is_snp <- (nchar(df$ref) == 1 & nchar(df$alt) == 1)
gr_snp <- gr_all[is_snp]

# 用 Bioconductor 推荐接口：snpsByOverlaps
library(GenomicFeatures)
snps_hit <- snpsByOverlaps(snp_db, gr_snp)

snps_hit
# 返回对象里有 RefSNP_id（rs号后面的数字）


#把 rsID 合回去（按位置匹配；如要更严谨可再比对 ref/alt）：
# snps_hit 的位点是单碱基，start=pos
hit_df <- data.frame(
  chr = as.character(seqnames(snps_hit)),
  pos = start(snps_hit),
  rs  = paste0("", snps_hit$RefSNP_id),
  stringsAsFactors = FALSE
)

# 回到原 df（注意 df 的 chr 可能是 chr20；但现在我们已经跟 db 风格统一过一次了）
# 保险做法：再建一个同风格的 chr 用于 merge
df2 <- df
df2$chr2 <- df2$chr
# 让 df2 的 chr2 也变成和 snp_db 一样的风格
tmp_gr <- GRanges(df2$chr2, IRanges(df2$pos, width = 1))
seqlevelsStyle(tmp_gr) <- seqlevelsStyle(snp_db)
df2$chr2 <- as.character(seqnames(tmp_gr))

df2$rsID <- NA_character_
idx <- match(paste(df2$chr2, df2$pos), paste(hit_df$chr, hit_df$pos))
df2$rsID[is_snp] <- hit_df$rs[idx[is_snp]]






