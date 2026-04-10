#@2026-02-01
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# 1) 输入：细胞类型计数
df_ct <- tribble(
  ~cell_type,              ~n,
  "Glutamatergic neuron",  109285,
  "GABAergic neuron",      35864,
  "Oligodendrocyte",       24613,
  "Radial glia",           16662,
  "Astrocyte",             11985,
  "OPC",                   12105,
  "IPC-EN",                9058,
  "Microglia",             6446,
  "Vascular",              2572,
  "IPC-Glia",              1551,
  "Cajal-Retzius cell",    865
)

# 2) 计算比例 & 排序
df_ct2 <- df_ct %>%
  mutate(
    total = sum(n),
    prop  = n / total,
    pct   = 100 * prop
  ) %>%
  arrange(desc(n)) %>%
  mutate(cell_type = factor(cell_type, levels = rev(cell_type)))  # <- 关键：反转顺序


# 3) 画图：组成柱状图（按占比从高到低）
p <- ggplot(df_ct2, aes(x = cell_type, y = prop)) +
  geom_col(width = 0.75) +
  coord_flip() +
  geom_text(
    aes(label = sprintf("%s (%.1f%%)", comma(n), pct)),
    hjust = -0.05, size = 3.6
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "Cell-type composition",
    x = NULL,
    y = "Proportion of cells"
  ) +
  theme_classic(base_size = 12)

print(p)

# 4)（可选）导出图片
# ggsave("celltype_composition_barplot.pdf", p, width = 8, height = 4.8)
# ggsave("celltype_composition_barplot.png", p, width = 8, height = 4.8, dpi = 300)

# 5)（可选）输出汇总表
df_ct2 %>%
  transmute(cell_type, n, pct = round(pct, 2)) %>%
  print(n = Inf)


suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# df_ct 需要你前面已经定义好：两列 cell_type / n

df_ct2 <- df_ct %>%
  mutate(
    total = sum(n),
    prop  = n / total,
    pct   = 100 * prop
  ) %>%
  arrange(desc(n)) %>%
  mutate(cell_type = factor(cell_type, levels = rev(cell_type)))  # 高的在上

p <- ggplot(df_ct2, aes(x = cell_type, y = prop, fill = prop)) +
  geom_col(width = 0.75) +
  coord_flip() +
  geom_text(
    aes(label = sprintf("%s (%.1f%%)", comma(n), pct)),
    hjust = -0.05, size = 3.6
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.15))
  ) +
  # 关键：按 proportion 做渐变色（橙色系）
  scale_fill_gradient(
    low = "#FBE6D4",
    high = "#E86E4A",
    name = "Proportion",
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    title = "Cell-type composition",
    x = NULL,
    y = "Proportion of cells"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

print(p)



