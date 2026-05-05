#2024-07-21




########_-------------------------------

setwd("/Users/mayunlong/Desktop/UPENN/00000------------UPenn-Projects/02-Aimed projects/05-devBrain-RBP-methods/RBP_data/")

install.packages("RMariaDB")


library(GenomicFeatures)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")

#提取3-UTR信息：---------------------------------------------
#构建 TxDb 对象
#从 UCSC 提取基因注释信息并构建一个 TxDb 对象：

txdb <- makeTxDbFromUCSC(genome="hg38",tablename = 'knownGene')

#提取3'UTR 信息

threeUTRs <- threeUTRsByTranscript(txdb,use.names=TRUE)

#将3UTR信息转换成GRanges对象

threeUTRs_gr <- unlist(threeUTRs,use.names = TRUE)

#创建BED文件
#将GRanges对象转换成BED文件并保存：
#将行名变成一列
threeUTRs_gr$name <- names(threeUTRs_gr)

#为每一行分配一个唯一的ID
row_ids <- paste0("ID_", seq_along(threeUTRs_gr))

#将行名替换成唯一的ID
names(threeUTRs_gr) <- row_ids

#将GRanges对象转换成data.frame
threeUTRs_df <- as.data.frame(threeUTRs_gr)


#提取需要的列，并格式化BED文件格式：

bed_data <- data.frame(chrom=threeUTRs_df$seqnames,
                       start=threeUTRs_df$start-1, #BED文件的起始位置为：0-based
                       end=threeUTRs_df$end,
                       name=threeUTRs_df$name,
                       score=".",
                       strand = threeUTRs_df$strand
)

#去除重复行，保留第一个出现的记录：
bed_data_unique <- bed_data[!duplicated(bed_data$name),]


#保存BED文件：
write.table(bed_data_unique,file="3UTRs_hg38.bed",quote=FALSE,sep = '\t',
            row.names = FALSE, col.names = FALSE)





#提取5-UTR信息：---------------------------------------------
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(GenomicRanges)


txdb <- makeTxDbFromUCSC(genome="hg38", tablename="knownGene")

fiveUTRs <- fiveUTRsByTranscript(txdb, use.names=TRUE)

fiveUTRs_gr <- unlist(fiveUTRs, use.names=TRUE)

#创建BED文件
#将GRanges对象转换成BED文件并保存：
#将行名变成一列
fiveUTRs_gr$name <- names(fiveUTRs_gr)

#为每一行分配一个唯一的ID
row_ids <- paste0("ID_", seq_along(fiveUTRs_gr))

#将行名替换成唯一的ID
names(fiveUTRs_gr) <- row_ids

#将GRanges对象转换成data.frame
fiveUTRs_df <- as.data.frame(fiveUTRs_gr)


#bed format
bed_data <- data.frame(
  chrom = fiveUTRs_df$seqnames,
  start = fiveUTRs_df$start - 1,  # BED 文件中的起始位置为 0-based
  end = fiveUTRs_df$end,
  name = fiveUTRs_df$name,
  score = ".",
  strand = fiveUTRs_df$strand
)

#去除重复行，保留第一个出现的记录：
bed_data_unique <- bed_data[!duplicated(bed_data$name),]


write.table(bed_data, file="5UTRs_hg38.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)



##CDS genomic regions bed format

txdb <- makeTxDbFromUCSC(genome="hg38", tablename="knownGene")

#提取CDS区域：
cds <- cdsBy(txdb, by="tx", use.names=TRUE)

#将CDS信息转换为GRangs对象
cds_gr <- unlist(cds, use.names=TRUE)



#创建BED文件
#将GRanges对象转换成BED文件并保存：
#将行名变成一列
cds_gr$name <- names(cds_gr)

#为每一行分配一个唯一的ID
row_ids <- paste0("ID_", seq_along(cds_gr))

#将行名替换成唯一的ID
names(cds_gr) <- row_ids

#将GRanges对象转换成data.frame
cds_df <- as.data.frame(cds_gr)


bed_data <- data.frame(
  chrom = cds_df$seqnames,
  start = cds_df$start - 1,  # BED 文件中的起始位置为 0-based
  end = cds_df$end,
  name = cds_df$name,
  score = ".",
  strand = cds_df$strand
)

#去除重复行，保留第一个出现的记录：
bed_data_unique <- bed_data[!duplicated(bed_data$name),]

write.table(bed_data, file="CDS_coordinate_hg38.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)





##-----提取内含子区域-intron

txdb <- makeTxDbFromUCSC(genome="hg38", tablename="knownGene")


# 提取内含子区域
introns <- intronsByTranscript(txdb,use.names=T)

# 转换为 GRanges 格式
introns_gr <- unlist(introns,use.names=T)


# 创建 BED 文件
# 将行名变成一列
introns_gr$name <- names(introns_gr)

# 为每一行分配一个唯一的ID
row_ids <- paste0("ID_", seq_along(introns_gr))

# 将行名替换成唯一的ID
names(introns_gr) <- row_ids

# 将GRanges对象转换成data.frame
introns_df <- as.data.frame(introns_gr)

# 创建 BED 文件格式的数据框
bed_data <- data.frame(
  chrom = introns_df$seqnames,
  start = introns_df$start - 1,  # BED 文件中的起始位置为 0-based
  end = introns_df$end,
  name = introns_df$name,
  score = ".",
  strand = introns_df$strand
)

# 去除重复行，保留第一个出现的记录
bed_data_unique <- bed_data[!duplicated(bed_data$name),]

# 保存为 BED 文件
write.table(bed_data_unique, file="introns_hg38.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)



