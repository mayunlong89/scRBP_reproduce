

df <- tibble::tibble(
  Disease = c("MS", "CD", "SLE", "PBC", "T1D", "RA", "UC", "IBD"),
  rho     = c(0.78, 0.68, 0.68, 0.66, 0.65, 0.65, 0.64, 0.61)
)


pal_bar <- c("#FAF3E3", "#FBB394", "#F87C56")

df2 <- df %>%
  mutate(Disease = fct_reorder(Disease, rho))  # 按 rho 排序（小到大，横向图更直观）

p <- ggplot(df2, aes(x = rho, y = Disease, fill = rho)) +
  geom_col(width = 0.85, color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("rho = %.2f", rho)),
            hjust = -0.05, size = 3.6) +
  scale_fill_gradientn(colors = pal_bar, name = expression(rho)) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "TRS reproducibility between datasets\n(170K and 10K cells)",
    x = "rho (TRS correlation)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0, face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

p
