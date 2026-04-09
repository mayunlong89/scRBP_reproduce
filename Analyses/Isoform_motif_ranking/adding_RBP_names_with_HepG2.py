import pandas as pd

# 读取原始 motif-by-RBP matrix 文件
df = pd.read_csv("motif_by_rbp_matrix_5mingenes.csv")

# 修改列名，第一列是 'Motif' 保留不变，其余加 '_HepG2'
df.columns = [df.columns[0]] + [col + "_HepG2" for col in df.columns[1:]]

# 保存修改后的文件
df.to_csv("motif_by_rbp_matrix_5mingenes_HepG2.csv", index=False)


