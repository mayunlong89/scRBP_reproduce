import pandas as pd

# 1. Assume you already have the following two data frames:
# - corr_df: correlation matrix of RBP x Gene (rows are RBPs, columns are Genes)
# - grn_df: long-format table output from GRNBoost2, with three columns: RBP, Gene, importance

# Example input (please replace with your actual file names)
corr_df = pd.read_parquet("RBP_gene_spearman_corr_200Kcells.parquet")  # Figure 1
grn_df = pd.read_csv("../z_scRBP_grn_network_200Kcells_final_v4_seed42_memory.tsv", sep="\t")  # Figure 2

# 2. Add correlation values to the GRNBoost2 table
# Use row-wise apply to extract the correlation value from corr_df
def get_corr(row):
    rbp, gene = row['RBP'], row['Gene']
    try:
        return corr_df.loc[rbp, gene]
    except KeyError:
        return float('nan')


# Add the Correlation column
grn_df['Correlation'] = grn_df.apply(get_corr, axis=1)

# 3. Optional: remove rows with NaN correlation values
grn_df_filtered = grn_df.dropna(subset=['Correlation'])

# Keep only edges with |Correlation| > 0.03 (i.e., strongly correlated edges)
grn_df_filtered = grn_df_filtered[
    (grn_df_filtered["Correlation"] > 0.03) |
    (grn_df_filtered["Correlation"] < -0.03)
]


# 4. Save results
grn_df_filtered.to_csv("RBP_gene_with_correlation_with_filtered.csv", index=False)
grn_df.to_csv("RBP_gene_with_all_correlation_without_filtered.csv", index=False)



###----plot----for grn_df_filtered
#----- Importance and correlation in filtered links
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
grn_df_filtered = pd.read_csv("z_GRNBoost2_result_200Kcells_correlation_v6_update4_with_correlation.tsv", sep="\t") 

# Create figure
plt.figure(figsize=(7, 6))

# Create 2D histogram
h = plt.hist2d(
    grn_df_filtered["Correlation"],
    grn_df_filtered["Importance"],
    bins=100,
    norm=LogNorm()
)

# Add colorbar
plt.colorbar(h[3], label="Count")

# Set labels
plt.xlabel("Correlation")
plt.ylabel("Importance")
plt.title("2D Histogram of Correlation vs. Importance (Filtered Links)")
plt.tight_layout()

# Save figure
plt.savefig("01_key_2Dhist_corr_vs_importance_filtered.png", dpi=300)
plt.close()



###------ plot all RBP-gene links correlation vs importances
#---- Plot all links (before filtering) with correlation and importance threshold lines
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
grn_df = pd.read_csv("z_GRNBoost2_result_200Kcells_correlation_v6_update4_with_correlation.tsv", sep="\t") 

plt.figure(figsize=(7, 6))

# 2D histogram of all RBP-gene links
h = plt.hist2d(
    grn_df["Correlation"],
    grn_df["Importance"],
    bins=100,
    norm=LogNorm()
)

# colorbar
plt.colorbar(h[3], label="Count")

# Add three threshold lines
plt.axvline(x=-0.03, color="red", linestyle="--", linewidth=1.5, label="ρ = -0.03")
plt.axvline(x=0.03, color="red", linestyle="--", linewidth=1.5, label="ρ = 0.03")
#plt.axhline(y=0.005, color="blue", linestyle="--", linewidth=1.5, label="Importance = 0.05")

# Figure labels and layout
plt.xlabel("Correlation")
plt.ylabel("Importance")
plt.title("2D Histogram of Correlation vs. Importance (All Links)")
plt.legend(loc="upper right")
plt.tight_layout()

# Save figures
plt.savefig("02_2Dhist_corr_vs_importance_all_with_thresholds.png", dpi=300)
plt.savefig("02_2Dhist_corr_vs_importance_all_with_thresholds.pdf",
            dpi=300,          # dpi has little effect on vector PDF, optional
            format="pdf")     # this can also be omitted; Matplotlib will infer it from the suffix

plt.close()









###---------Examine the importance distribution of neutral links:
###---------Examine the importance distribution of neutral links:
###---------Examine the importance distribution of neutral links:
###---------Examine the importance distribution of neutral links:


import matplotlib.pyplot as plt
import seaborn as sns

# Original data: grn_df contains the full unfiltered dataset
# Remove strongly correlated links (|correlation| > 0.03) — the remaining ones are neutral associations
neutral_links = grn_df[
    (grn_df["correlation"] >= -0.03) & (grn_df["correlation"] <= 0.03)
]

# Optional: save to file
neutral_links.to_csv("neutral_RBP_gene_links.csv", index=False)


plt.figure(figsize=(8, 5))

# Main plot: importance distribution (log scale)
sns.histplot(neutral_links["Importance"], bins=100, kde=True, color="gray")

# Add vertical line
plt.axvline(x=0.05, color="red", linestyle="--", linewidth=1.5, label="Importance = 0.05")

# Figure settings
plt.title("Importance Score Distribution of Neutral RBP-Gene Links")
plt.xlabel("Importance")
plt.ylabel("Count")
plt.xscale("log")
plt.legend()
plt.tight_layout()

# Save
plt.savefig("neutral_importance_with_threshold.png", dpi=300)
plt.close()



##-----Importance and correlation corelation
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm  # ✅ Correct import
import numpy as np

# Create figure
plt.figure(figsize=(7, 6))

# Create 2D histogram and capture the returned object
h = plt.hist2d(
    neutral_links["correlation"],
    neutral_links["Importance"],
    bins=100,
    norm=LogNorm()  # ✅ Use imported LogNorm
)

# Add colorbar
plt.colorbar(h[3], label="Count")

# Set figure labels
plt.xlabel("Correlation")
plt.ylabel("Importance")
plt.title("2D Histogram of Correlation vs. Importance (Neutral Links)")
plt.tight_layout()

# Save figure
plt.savefig("2Dhist_corr_vs_importance_neutral_fixed.png", dpi=300)
plt.close()



# Count totals -------- proportion
# Count totals -------- proportion
# Count totals -------- proportion
# Count totals -------- proportion

# 1. Neutral links
total_neutral = len(neutral_links)

# Count the number of links with Importance > 0.005
high_importance = (neutral_links["Importance"] > 0.005).sum()

# Calculate fraction
fraction = high_importance / total_neutral

# Output results
print(f"Total neutral links: {total_neutral}")
print(f"Links with Importance > 0.005: {high_importance}")
print(f"Fraction: {fraction:.4%}")


# 2. Count total number of links (unfiltered)
total_links = len(grn_df)

# Count the number of links with Importance > 0.005
high_importance = (grn_df["Importance"] > 0.005).sum()

# Calculate fraction
fraction = high_importance / total_links

# Output results
print(f"Total links: {total_links}")
print(f"Links with Importance > 0.005: {high_importance}")
print(f"Fraction: {fraction:.4%}")



activating = grn_df_filtered[grn_df_filtered['correlation'] > 0.03]
repressing = grn_df_filtered[grn_df_filtered['correlation'] < -0.03]






##3.1.1 plot volcano plot for activating repressing links
###------ plot all RBP-gene links correlation vs importances
#---- Plot all links (before filtering) with correlation and importance threshold lines
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

plt.figure(figsize=(7, 6))

# 所有 RBP-gene links 的 2D 直方图
h = plt.hist2d(
    grn_df["correlation"],
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


##3.1 boxplot for three group: Repressing links, neutral links, and activating links

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# 1️⃣ 过滤低 importance 边
grn_df_filtered_plot = grn_df[grn_df["Importance"] >= 0.005].copy()

# 2️⃣ 按 correlation 归类
def classify_link(rho):
    if rho < -0.03:
        return "Repressing"
    elif rho > 0.03:
        return "Activating"
    else:
        return "Neutral"

grn_df_filtered_plot["Group"] = grn_df_filtered_plot["correlation"].apply(classify_link)

# 3️⃣ 计算两两比较的 p 值（Mann–Whitney U）
activating_vals = grn_df_filtered_plot.loc[grn_df_filtered_plot["Group"] == "Activating", "Importance"]
repressing_vals  = grn_df_filtered_plot.loc[grn_df_filtered_plot["Group"] == "Repressing",  "Importance"]
neutral_vals     = grn_df_filtered_plot.loc[grn_df_filtered_plot["Group"] == "Neutral",     "Importance"]

p_act_neu = mannwhitneyu(activating_vals, neutral_vals, alternative="two-sided").pvalue
p_rep_neu = mannwhitneyu(repressing_vals,  neutral_vals, alternative="two-sided").pvalue

print(f"Activating vs Neutral p = {p_act_neu:.3e}")
print(f"Repressing  vs Neutral p = {p_rep_neu:.3e}")

# 4️⃣ 绘制 boxplot（自定义顺序）
order = ["Neutral", "Repressing", "Activating"]
palette = {"Neutral":"gray", "Repressing":"steelblue", "Activating":"darkorange"}

plt.figure(figsize=(8, 6))
ax = sns.boxplot(
    data=grn_df_filtered_plot,
    x="Group",
    y="Importance",
    order=order,
    palette=palette
)
plt.yscale("log")
plt.title("Importance Score Distribution by Correlation Group\n(Importance ≥ 0.005)")
plt.xlabel("RBP-Gene Link Type")
plt.ylabel("Importance Score (log scale)")

# 5️⃣ 在图上标注显著性
y_max = grn_df_filtered_plot["Importance"].max()
y_step = y_max * 0.15          # 控制标注高度
# Activating vs Neutral
ax.plot([0, 2], [y_max, y_max], color="black", linewidth=1)
ax.text(1, y_max*1.05, f"p = {p_act_neu:.1e}", ha="center")
# Repressing vs Neutral
ax.plot([0, 1], [y_max - y_step, y_max - y_step], color="black", linewidth=1)
ax.text(0.5, (y_max - y_step)*1.05, f"p = {p_rep_neu:.1e}", ha="center")

plt.tight_layout()
plt.savefig("05_boxplot_importance_by_group_with_stats.png", dpi=300)
plt.savefig("05_boxplot_importance_by_group_with_stats_python.pdf", dpi=300)
plt.close()




