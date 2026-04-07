
###------ plot all RBP-gene links correlation vs importances
#---- Plot all links (before filtering) with correlation and importance threshold lines
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
grn_df = pd.read_csv("isoT4k_grn_seed1234_v2_scRBP_isoform_GRNs.tsv", sep="\t") 

plt.figure(figsize=(7, 6))

# 所有 RBP-gene links 的 2D 直方图
h = plt.hist2d(
    grn_df["Correlation"],
    grn_df["Importance"],
    bins=100,
    norm=LogNorm()
)

# colorbar
plt.colorbar(h[3], label="Count")

# 添加三条阈值线
plt.axvline(x=-0.03, color="red", linestyle="--", linewidth=1.5, label="ρ = -0.03")
plt.axvline(x=0.03, color="red", linestyle="--", linewidth=1.5, label="ρ = 0.03")
#plt.axhline(y=0.005, color="blue", linestyle="--", linewidth=1.5, label="Importance = 0.05")

# 图形标签与布局
plt.xlabel("Correlation")
plt.ylabel("Importance")
plt.title("2D Histogram of Correlation vs. Importance (All Links)")
plt.legend(loc="upper right")
plt.tight_layout()

# 保存图像
plt.savefig("02_2Dhist_corr_vs_importance_all_with_thresholds.png", dpi=300)
plt.savefig("02_2Dhist_corr_vs_importance_all_with_thresholds.pdf",
            dpi=300,          # dpi 对矢量 PDF 影响不大，可选
            format="pdf")     # 也可以省略，Matplotlib 会根据后缀自动识别

plt.close()
