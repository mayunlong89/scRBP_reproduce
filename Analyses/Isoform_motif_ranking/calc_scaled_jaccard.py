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

