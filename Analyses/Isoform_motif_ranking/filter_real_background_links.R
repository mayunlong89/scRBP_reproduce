# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript script.R <enrich_file> <link_file> <real_output> <background_output>")
}

enrich_file <- args[1]
link_file <- args[2]
real_output <- args[3]
background_output <- args[4]

# 读取文件
enrich <- read.csv(enrich_file)
link <- read.csv(link_file)

# 创建唯一键: Motif + RBP
enrich$key <- paste(enrich$Motif, enrich$RBP, sep = "::")
link$key <- paste(link$Motif, link$RBP, sep = "::")

# 筛选真实 link 和背景 link
real_links <- enrich[enrich$key %in% link$key, ]
background_links <- enrich[!(enrich$key %in% link$key), ]

# 保存
write.csv(real_links, real_output, row.names = FALSE)
write.csv(background_links, background_output, row.names = FALSE)
