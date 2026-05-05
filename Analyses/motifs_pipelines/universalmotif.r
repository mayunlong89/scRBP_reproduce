#!/usr/bin/env Rscript

# 加载需要的库
library(universalmotif)

# 读取命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript script.R <input_file> <output_file>")
}

input_file <- args[1]
output_file <- args[2]

# 读取STREME生成的PWM文件
streme_pwm <- read_meme(input_file)

# 定义修剪和舍入函数
trim_and_round <- function(motif, min_ic = 0.1, min_prob = 0.025) {
  # 修剪PWM motif，去除信息含量（IC）<= min_ic的部分
  trimmed_motif <- trim_motifs(motif, min.ic = min_ic)
  
  # 舍入，将基因概率值 <= min_prob的设为0
  pwm <- trimmed_motif@motif
  pwm[pwm <= min_prob] <- 0
  trimmed_motif@motif <- pwm
  
  return(trimmed_motif)
}

# 修剪和舍入概率值
trimmed_motifs <- lapply(streme_pwm, trim_and_round, min_ic = 0.1, min_prob = 0.025)

# 保存修剪后的PWM文件
write_meme(trimmed_motifs, output_file)
