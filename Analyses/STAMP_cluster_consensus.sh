
#2024-10-25

#1）将每个cluster motifs 准备为单独的文件：

#这里需要将每个cluster中的motifs分别保存到不同的文件中。

#2）使用STAMP对每个cluster内的motifs进行比对：

#可以将每个cluster中的motifs作为输入文件，使用STAMP进行相似性比对。
#这一步可以帮助找到cluster内的motifs之间的相似性，并生成每个cluster的共识motifs (consensus motifs).

stamp -tf cluster_1.meme -oc cluster_1_output -cc PCC -align SWU

-tf 表示输入 motif 文件。
-oc 是输出文件夹，用于存储比对结果。
-cc 是相似性评分方法，您可以使用 Pearson correlation coefficient（PCC）。
-align SWU 使用 Smith-Waterman 对齐算法。


#3) 生成共识motifs并对比clusters
#对于每个cluster，STAMP可以输出一个共识motif文件；
#可以进一步对不同cluster的共识motifs进行比对，评估cluster之间的相似性；

#use STAMP to compare cluster-based consensus motifs and generate dendrogram plot

stamp -tf cluster_consensus.meme -oc consensus_output -cc PCC -align SWU

#4) Dendrogram
#STAMP 可以自动生成 motifs 之间的相似性树形图，用于可视化 cluster 之间的关系。
#可以从每个 cluster 的共识 motifs 开始生成树形图，显示 cluster 之间的层次关系和相似性。

#使用工具如 iTOL（Interactive Tree of Life）来可视化和注释生成的树。
#选择使用 matplotlib 或 seaborn 等 Python 库来生成个性化的树状图，并标注 cluster 内 motifs 的分布和特征。

###batch running:
#!/bin/bash

for i in {1..184}
do
    stamp -tf cluster_$i.meme -oc cluster_${i}_output -cc PCC -align SWU
done


###Final version for running STAMP:

stamp -tf input.meme -oc output_directory -cc PCC -sd GUMBEL -chp 3 -align SWU



