
#conda evironment under: pyscenic

# srun --mem=80G --pty Rscript ./CDS_extract_bed_script.r



##code for extract
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(GenomicRanges)



txdb <- makeTxDbFromGFF("/mnt/isilon/gandal_lab/mayl/reference/gencode.v44.primary_assembly.annotation.gtf",format="gtf")


#提取CDS区域：
cds <- cdsBy(txdb, by="tx", use.names=TRUE)

#将CDS信息转换为GRangs对象
cds_gr <- unlist(cds, use.names=TRUE)



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

