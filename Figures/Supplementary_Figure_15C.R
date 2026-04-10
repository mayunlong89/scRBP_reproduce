
#2025-12-15

library(tidyverse)

##--------IBD
# df 需要包含列：`170K_TRS`, `10K_TRS`
 df <- readr::read_csv("01_IBD_replication.csv", show_col_types = FALSE)

rho <- cor(df$`170K_TRS`, df$`10K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `170K_TRS`, y = `10K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 1.2, alpha = 0.35, color = "red") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.4, linetype = "dashed",color = "#2F6BFF") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    axis.line  = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  labs(
    title = sprintf("IBD: 170K vs 10K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "170K TRS",
    y = "10K TRS"
  )

p_scatter





##--------CD
# df 需要包含列：`170K_TRS`, `10K_TRS`
df <- readr::read_csv("01_CD_replication.csv", show_col_types = FALSE)

rho <- cor(df$`170K_TRS`, df$`10K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `170K_TRS`, y = `10K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 1.2, alpha = 0.35, color = "red") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.4, linetype = "dashed",color = "#2F6BFF") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    axis.line  = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  labs(
    title = sprintf("CD: 170K vs 10K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "170K TRS",
    y = "10K TRS"
  )

p_scatter






##--------MS
# df 需要包含列：`170K_TRS`, `10K_TRS`
df <- readr::read_csv("01_MS_replication.csv", show_col_types = FALSE)

rho <- cor(df$`170K_TRS`, df$`10K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `170K_TRS`, y = `10K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 1.2, alpha = 0.35, color = "red") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.4, linetype = "dashed",color = "#2F6BFF") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    axis.line  = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  labs(
    title = sprintf("MS: 170K vs 10K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "170K TRS",
    y = "10K TRS"
  )

p_scatter





##--------UC
# df 需要包含列：`170K_TRS`, `10K_TRS`
df <- readr::read_csv("01_UC_replication.csv", show_col_types = FALSE)

rho <- cor(df$`170K_TRS`, df$`10K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `170K_TRS`, y = `10K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 1.2, alpha = 0.35, color = "red") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.4, linetype = "dashed",color = "#2F6BFF") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    axis.line  = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  labs(
    title = sprintf("UC: 170K vs 10K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "170K TRS",
    y = "10K TRS"
  )

p_scatter





##--------RA
# df 需要包含列：`170K_TRS`, `10K_TRS`
df <- readr::read_csv("01_RA_replication.csv", show_col_types = FALSE)

rho <- cor(df$`170K_TRS`, df$`10K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `170K_TRS`, y = `10K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 1.2, alpha = 0.35, color = "red") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.4, linetype = "dashed",color = "#2F6BFF") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    axis.line  = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  labs(
    title = sprintf("RA: 170K vs 10K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "170K TRS",
    y = "10K TRS"
  )

p_scatter






##--------SLE
# df 需要包含列：`170K_TRS`, `10K_TRS`
df <- readr::read_csv("01_SLE_replication.csv", show_col_types = FALSE)

rho <- cor(df$`170K_TRS`, df$`10K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `170K_TRS`, y = `10K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 1.2, alpha = 0.35, color = "red") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.4, linetype = "dashed",color = "#2F6BFF") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    axis.line  = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  labs(
    title = sprintf("SLE: 170K vs 10K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "170K TRS",
    y = "10K TRS"
  )

p_scatter




##--------T1D
# df 需要包含列：`170K_TRS`, `10K_TRS`
df <- readr::read_csv("01_T1D_replication.csv", show_col_types = FALSE)

rho <- cor(df$`170K_TRS`, df$`10K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `170K_TRS`, y = `10K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 1.2, alpha = 0.35, color = "red") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.4, linetype = "dashed",color = "#2F6BFF") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    axis.line  = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  labs(
    title = sprintf("T1D: 170K vs 10K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "170K TRS",
    y = "10K TRS"
  )

p_scatter




##--------PBC
# df 需要包含列：`170K_TRS`, `10K_TRS`
df <- readr::read_csv("01_PBC_replication.csv", show_col_types = FALSE)

rho <- cor(df$`170K_TRS`, df$`10K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `170K_TRS`, y = `10K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 1.2, alpha = 0.35, color = "red") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.4, linetype = "dashed",color = "#2F6BFF") +
  coord_cartesian(expand = FALSE) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    axis.line  = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  labs(
    title = sprintf("PBC: 170K vs 10K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "170K TRS",
    y = "10K TRS"
  )

p_scatter

