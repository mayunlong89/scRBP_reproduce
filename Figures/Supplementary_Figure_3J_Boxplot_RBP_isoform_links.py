
################################################
#########----> Python version

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import pandas as pd

grn_df = pd.read_csv("isoT4k_grn_seed1234_v2_scRBP_isoform_GRNs.tsv", sep="\t") 

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

grn_df_filtered_plot["Group"] = grn_df_filtered_plot["Correlation"].apply(classify_link)

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


################################################
#########----> R version

# ── 加载所需包 ──────────────────────────────────────────
library(dplyr)
library(ggplot2)

# grn_df 假设已经在环境中；若需读取请用 read.csv()/readr::read_csv()

grn_df <- read.table("isoT4k_grn_seed1234_v2_scRBP_isoform_GRNs.tsv", header=T)

# 1️⃣ 过滤低 importance
grn_df_filtered <- grn_df %>%
  filter(Importance >= 0.005)

# 2️⃣ 分类
grn_df_filtered <- grn_df_filtered %>% 
  mutate(Group = case_when(
    Correlation < -0.03 ~ "Repressing",
    Correlation >  0.03 ~ "Activating",
    TRUE                ~ "Neutral"
  ))

# 3️⃣ 计算 p 值（Wilcoxon / Mann–Whitney U）
act_vals <- grn_df_filtered %>% filter(Group == "Activating") %>% pull(Importance)
rep_vals <- grn_df_filtered %>% filter(Group == "Repressing")  %>% pull(Importance)
neu_vals <- grn_df_filtered %>% filter(Group == "Neutral")     %>% pull(Importance)

#p_act_neu <- wilcox.test(act_vals, neu_vals)$p.value
#p_rep_neu <- wilcox.test(rep_vals, neu_vals)$p.value

# 4️⃣ 绘图
order_levels <- c("Neutral", "Repressing", "Activating")
pal          <- c(Neutral = "gray", Repressing = "steelblue", Activating = "darkorange")

y_max  <- max(grn_df_filtered$Importance)
y_step <- y_max * 0.15              # 控制两条显著性横线的间距

p <- ggplot(grn_df_filtered, aes(x = factor(Group, levels = order_levels),
                                 y = Importance, fill = Group)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.5) +
  scale_fill_manual(values = pal) +
  scale_y_log10() +
  labs(title = "Importance Score Distribution by Correlation Group\n(Importance ≥ 0.005)",
       x = "RBP-Gene Link Type",
       y = "Importance Score (log scale)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") 

  # —— 添加显著性横线与文字 —— 
 # geom_segment(aes(x = 1, xend = 3, y = y_max,       yend = y_max)) +
 # annotate("text", x = 2, y = y_max * 1.05,
 #          label = sprintf("p = %.1e", p_act_neu)) +
 # geom_segment(aes(x = 1, xend = 2, y = y_max - y_step, yend = y_max - y_step)) +
 # annotate("text", x = 1.5, y = (y_max - y_step) * 1.05,
  #         label = sprintf("p = %.1e", p_rep_neu))

# 5️⃣ 保存图像：PNG（改后缀即可存为 PDF）
ggsave("05_boxplot_importance_by_group_with_stats_R.png",
       plot = p, width = 6, height = 5, dpi = 300)








