

#2025-12-24

library(dplyr)
library(ggplot2)

# 0) 读入
cnt <- read.csv("01_8autoimmune_circle_stacking_plot.csv", stringsAsFactors = FALSE)

# 1) 统一 region 顺序（可改）
cnt$region <- factor(cnt$region, levels = c("5UTR","CDS","Introns","3UTR"))

# 2) 计算每个 disease 的 total，并加回到表里
cnt <- cnt %>%
  group_by(disease) %>%
  mutate(total = sum(n)) %>%
  ungroup()

# 3) disease 顺序：按 total 从大到小（或者你自己指定论文顺序）
disease_order <- cnt %>%
  distinct(disease, total) %>%
  arrange(desc(total)) %>%
  pull(disease)

cnt$disease <- factor(cnt$disease, levels = disease_order)


region_cols <- c(
  "5UTR"     = "#E07A5F",  # muted coral
  "CDS"      = "#7FB069",  # olive green
  "Introns"  = "#5DB7B7",  # teal
  "3UTR"     = "#B084F5"   # soft purple
)

# 4) 画圆形堆叠（donut）
hole <- 0.25 * max(cnt$total)  # 中间空洞大小，可调 0.2~0.4

p <- ggplot(cnt, aes(x = disease, y = n, fill = region)) +
  geom_col(width = 0.95, color = "white", linewidth = 0.2) +
  # 甜甜圈洞：给 y 轴留负值
  scale_y_continuous(limits = c(-hole, max(cnt$total) * 1.12)) +
  coord_polar(theta = "x", start = pi/2) +
  # 外圈 total 标注
  geom_text(
    data = cnt %>% distinct(disease, total),
    aes(x = disease, y = total + 0.04 * max(total), label = total),
    inherit.aes = FALSE, size = 3
  ) +
  labs(x = NULL, y = NULL, fill = "Region") +
  theme_void() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 10),
    plot.margin = margin(10, 20, 10, 20)
  )

p
p +
  scale_fill_manual(values = region_cols)
