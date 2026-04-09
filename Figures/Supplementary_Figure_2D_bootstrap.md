
### 1. filter_real_background_links.R

```R
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
```

### 2.calc_scaled_jaccard.py
```python
import pandas as pd
import numpy as np
import argparse

# 1. 获取命令行参数
parser = argparse.ArgumentParser(description='Calculate scRBP Jaccard similarity and z-score normalization.')
parser.add_argument('--k562_real', required=True, help='Filtered real links for K562')
parser.add_argument('--hepg2_real', required=True, help='Filtered real links for HepG2')
parser.add_argument('--k562_bg', required=True, help='Background eCLIP targets for K562')
parser.add_argument('--hepg2_bg', required=True, help='Background eCLIP targets for HepG2')
parser.add_argument('--output', required=True, help='Output file for scaled Jaccard similarity')
args = parser.parse_args()

# 2. 读取 scRBP prune 结果
k562_df = pd.read_csv(args.k562_real)
hepg2_df = pd.read_csv(args.hepg2_real)

k562_df['motif_rbp'] = k562_df['Motif'] + '::' + k562_df['RBP']
hepg2_df['motif_rbp'] = hepg2_df['Motif'] + '::' + hepg2_df['RBP']

common_keys = set(k562_df['motif_rbp']).intersection(set(hepg2_df['motif_rbp']))
k562_common = k562_df[k562_df['motif_rbp'].isin(common_keys)].set_index('motif_rbp')
hepg2_common = hepg2_df[hepg2_df['motif_rbp'].isin(common_keys)].set_index('motif_rbp')

def jaccard(set1, set2):
    if not set1 and not set2:
        return 1.0
    return len(set1 & set2) / len(set1 | set2)

sc_jaccard_scores = []
for key in common_keys:
    genes_k562 = set(str(k562_common.loc[key, 'Leading_Edge_Genes']).split(','))
    genes_hepg2 = set(str(hepg2_common.loc[key, 'Leading_Edge_Genes']).split(','))
    jaccard_sim = jaccard(genes_k562, genes_hepg2)
    sc_jaccard_scores.append((key, jaccard_sim))

sc_jaccard_df = pd.DataFrame(sc_jaccard_scores, columns=['motif_rbp', 'scRBP_Jaccard'])

# 3. 读取背景 eCLIP 结果
background_k562 = pd.read_csv(args.k562_bg)
background_hepg2 = pd.read_csv(args.hepg2_bg)

rbps = set(background_k562['RBP']).intersection(set(background_hepg2['RBP']))

background_jaccards = []
for rbp in rbps:
    genes_k562 = set(background_k562[background_k562['RBP'] == rbp]['Gene'])
    genes_hepg2 = set(background_hepg2[background_hepg2['RBP'] == rbp]['Gene'])
    if genes_k562 or genes_hepg2:
        jac = jaccard(genes_k562, genes_hepg2)
        background_jaccards.append(jac)

bg_mean = np.mean(background_jaccards)
bg_std = np.std(background_jaccards)

# 4. Z-score 标准化
sc_jaccard_df['Z_Scaled_Jaccard'] = (sc_jaccard_df['scRBP_Jaccard'] - bg_mean) / bg_std

# 5. 保存
sc_jaccard_df.to_csv(args.output, index=False)
```
### 3. bootstrap_cal.R
```R
# 1. 解析命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript script.R <real_file> <fake_file> <metric> <output_plot> <output_pval>")
}

real_file <- args[1]
fake_file <- args[2]
metric <- args[3]  # "scRBP_Jaccard" or "Z_Scaled_Jaccard"
output_plot <- args[4]
output_pval <- args[5]

# 2. 读取数据
real_df <- read.csv(real_file)
fake_df <- read.csv(fake_file)

# 检查 metric 合法性
if (!(metric %in% c("scRBP_Jaccard", "Z_Scaled_Jaccard"))) {
  stop("Metric must be 'scRBP_Jaccard' or 'Z_Scaled_Jaccard'")
}

# 3. 计算真实均值
real_mean <- mean(real_df[[metric]])

# 4. Bootstrap
set.seed(123)
bootstrap_means <- replicate(1000, {
  fake_sample <- fake_df[sample(nrow(fake_df), nrow(real_df)), ]
  mean(fake_sample[[metric]])
})

# 5. 绘图
library(ggplot2)

boot_df <- data.frame(Mean = bootstrap_means)

p <- ggplot(boot_df, aes(x = Mean)) +
  geom_histogram(bins = 30, fill = "grey", color = "black") +
  geom_vline(xintercept = real_mean, color = "red", linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = paste("Bootstrap Null Distribution of Fake Means (", metric, ")", sep = ""),
       x = paste("Mean", metric), y = "Frequency",
       subtitle = paste("Red dashed line = Real Mean =", round(real_mean, 3)))

ggsave(output_plot, plot = p, width = 6, height = 4)

# 6. 计算 p 值 (双尾)
p_value <- mean(abs(bootstrap_means) >= abs(real_mean))
write(paste("P-value:", p_value), file = output_pval)



```


### 4. run bootstrap_cal.R
```bash
#---2025-04-28---
#bootstrap analysis for all conditions
#we can use original Jaccard similarity scores or Z-scaled Jaccard similarity scores (adjusted Jaccard Simialrity (AJS))
####----------------------------------------------------------------------
####----------------------------------------------------------------------
####--------------Example-------------------------------------------
####----------------------------------------------------------------------
## Original JS:  scRBP_Jaccard ---不常用
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/3UTR
Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/FIMO_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/FIMO_K562_HepG2_scaled_jaccard_similarity_background.csv \
  scRBP_Jaccard \
  $output_DIR/bootstrap_null_distribution_scRBP_Jaccard.pdf \
  $output_DIR/pval_scRBP_Jaccard.txt


## adjusted Jaccard Simialrity (AJS):  Z_Scaled_Jaccard
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/FIMO_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/FIMO_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/bootstrap_null_distribution_scaled.pdf \
  $output_DIR/pval_scaled.txt
####----------------------------------------------------------------------
####----------------------------------------------------------------------
####----------------------------------------------------------------------
####----------------------------------------------------------------------



###----Real condition


####-----3UTR------------
####-----3UTR------------
####-----3UTR------------
####-----3UTR------------
###----FIMO--
# Calculate-Bootstrap-----FIMO---3UTR---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/3UTR

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/FIMO_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/FIMO_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/FIMO_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/FIMO_K562_HepG2_pval_scaled.txt


###----HOMER
# Calculate-Bootstrap-----HOMER---3UTR---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/3UTR

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/HOMER_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/HOMER_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/HOMER_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/HOMER_K562_HepG2_pval_scaled.txt




###----Clustered
# Calculate-Bootstrap-----Clustered---3UTR---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/3UTR

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/Clustered_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/Clustered_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/Clustered_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/Clustered_K562_HepG2_pval_scaled.txt


###----individual
# Calculate-Bootstrap-----individual---3UTR---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/3UTR

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/individual_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/individual_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/individual_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/individual_K562_HepG2_pval_scaled.txt
  
  
###----archetype
# Calculate-Bootstrap-----archetype---3UTR---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/3UTR

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/archetype_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/archetype_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/archetype_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/archetype_K562_HepG2_pval_scaled.txt
  
  
  
  
  
####-----5UTR------------
####-----5UTR------------
####-----5UTR------------
####-----5UTR------------
###----FIMO
# Calculate-Bootstrap-----FIMO---5UTR---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/5UTR

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/FIMO_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/FIMO_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/FIMO_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/FIMO_K562_HepG2_pval_scaled.txt


###----HOMER
# Calculate-Bootstrap-----HOMER---5UTR---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/5UTR

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/HOMER_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/HOMER_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/HOMER_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/HOMER_K562_HepG2_pval_scaled.txt




###----Clustered
# Calculate-Bootstrap-----Clustered---5UTR---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/5UTR

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/Clustered_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/Clustered_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/Clustered_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/Clustered_K562_HepG2_pval_scaled.txt


###----individual
# Calculate-Bootstrap-----individual---5UTR---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/5UTR

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/individual_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/individual_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/individual_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/individual_K562_HepG2_pval_scaled.txt
  
  
  
  
  ###----archetype
# Calculate-Bootstrap-----archetype---5UTR---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/5UTR

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/archetype_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/archetype_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/archetype_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/archetype_K562_HepG2_pval_scaled.txt
  
  
  
  
  
  

####-----CDS------------
####-----CDS------------
####-----CDS------------
####-----CDS------------
###----FIMO
# Calculate-Bootstrap-----FIMO---CDS---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/CDS

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/FIMO_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/FIMO_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/FIMO_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/FIMO_K562_HepG2_pval_scaled.txt


###----HOMER
# Calculate-Bootstrap-----HOMER---CDS---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/CDS

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/HOMER_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/HOMER_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/HOMER_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/HOMER_K562_HepG2_pval_scaled.txt




###----Clustered
# Calculate-Bootstrap-----Clustered---CDS---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/CDS

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/Clustered_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/Clustered_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/Clustered_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/Clustered_K562_HepG2_pval_scaled.txt


###----individual
# Calculate-Bootstrap-----individual---CDS---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/CDS

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/individual_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/individual_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/individual_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/individual_K562_HepG2_pval_scaled.txt
  
  
###----archetype
# Calculate-Bootstrap-----archetype---CDS---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/CDS

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/archetype_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/archetype_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/archetype_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/archetype_K562_HepG2_pval_scaled.txt
  



####-----Introns------------
####-----Introns------------
####-----Introns------------
####-----Introns------------
###----FIMO
# Calculate-Bootstrap-----FIMO---Introns---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/Introns

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/FIMO_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/FIMO_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/FIMO_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/FIMO_K562_HepG2_pval_scaled.txt


###----HOMER
# Calculate-Bootstrap-----HOMER---Introns--K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/Introns

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/HOMER_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/HOMER_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/HOMER_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/HOMER_K562_HepG2_pval_scaled.txt




###----Clustered
# Calculate-Bootstrap-----Clustered---Introns---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/Introns

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/Clustered_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/Clustered_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/Clustered_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/Clustered_K562_HepG2_pval_scaled.txt


###----individual
# Calculate-Bootstrap-----individual---Introns---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/Introns

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/individual_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/individual_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/individual_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/individual_K562_HepG2_pval_scaled.txt
  
  
###----archetype
# Calculate-Bootstrap-----archetype--Introns---K562 and HepG2
TOOL=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2
output_DIR=/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/05_benchmark_ENCODE/z_two_cell_lines/z_jaccard_similarity_nes_K562_HepG2/Introns

Rscript $TOOL/bootstrap_cal.R \
  $output_DIR/archetype_K562_HepG2_scaled_jaccard_similarity.csv \
  $output_DIR/archetype_K562_HepG2_scaled_jaccard_similarity_background.csv \
  Z_Scaled_Jaccard \
  $output_DIR/archetype_K562_HepG2_bootstrap_null_distribution_scaled.pdf \
  $output_DIR/archetype_K562_HepG2_pval_scaled.txt
  
```








