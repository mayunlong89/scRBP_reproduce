

#2026-1-31

library(dplyr)
library(ggplot2)

# 0) 读入
cnt <- read.csv("01_8brain_diseases_circle_stacking_plot.csv", stringsAsFactors = FALSE)

# 1) 统一 region 顺序（可改）
cnt$region <- factor(cnt$region, levels = c("5'UTR","CDS","Intron","3'UTR"))

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
  "5'UTR"     = "#E07A5F",  # muted coral
  "CDS"      = "#7FB069",  # olive green
  "Intron"  = "#5DB7B7",  # teal
  "3'UTR"     = "#B084F5"   # soft purple
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






###########_--------other version


library(ggplot2)
library(dplyr)


cnt <- read.csv("01_8brain_diseases_circle_stacking_plot.csv")

# 1) disease 顺序（按 total 从大到小更直观；也可以固定你论文的顺序）
disease_order <- cnt %>%
  distinct(disease, total) %>%
  arrange(desc(total)) %>%
  pull(disease)

cnt <- cnt %>%
  mutate(
    disease = factor(disease, levels = disease_order),
    region  = factor(region, levels = c("5'UTR","CDS","Intron","3'UTR"))
  )

# 2) 画圆形堆叠（关键：coord_polar）
hole <- 0.25 * max(cnt$total)   # 中间空洞大小（可调）

p_circ <- ggplot(cnt, aes(x = disease, y = n, fill = region)) +
  geom_col(width = 0.95, color = "white", linewidth = 0.3) +
  # 做“甜甜圈洞”：给 y 轴留负值空间
  scale_y_continuous(limits = c(-hole, max(cnt$total) * 1.10)) +
  coord_polar(theta = "x", start = pi/2) +
  # 标 total（放在外侧）
  geom_text(
    data = cnt %>% distinct(disease, total),
    aes(x = disease, y = total + 0.03 * max(total), label = total),
    inherit.aes = FALSE, size = 3
  ) +
  labs(x = NULL, y = NULL, fill = "Region") +
  theme_void() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 10),
    plot.margin = margin(10, 20, 10, 20)
  )

p_circ





###------version2

library(dplyr)
library(ggplot2)
library(ggforce)

# 把 disease 映射到角度（0~2pi），每个 disease 一个扇区
diseases <- levels(cnt$disease)
k <- length(diseases)

angle_df <- tibble(
  disease = diseases,
  id = seq_along(diseases),
  start = 2*pi*(id-1)/k,
  end   = 2*pi*id/k
)

# 给每个 disease 每个 region 分配“径向堆叠”的内外半径
ring <- cnt %>%
  left_join(angle_df, by="disease") %>%
  arrange(disease, region) %>%
  group_by(disease) %>%
  mutate(r0 = cumsum(lag(n, default = 0)),
         r1 = cumsum(n)) %>%
  ungroup()

# 归一化半径（可选：让最大 total 到 1，图更紧凑）
max_total <- max(ring$total)
ring <- ring %>%
  mutate(r0 = r0 / max_total,
         r1 = r1 / max_total)

p_force <- ggplot(ring) +
  geom_arc_bar(
    aes(x0=0, y0=0, r0=0.35 + r0*0.60, r=0.35 + r1*0.60, start=start, end=end, fill=region),
    color="white", linewidth=0.3
  ) +
  coord_fixed() +
  theme_void() +
  theme(legend.position="right")

p_force

