
#2026-02-02
setwd("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/03_scRBP_trs_fetal_brain_15times/01_Final_scRBP_TRS_fetal_brain_50times/01_8diseases_TRS_correlation_between_discovery_and_validation")
library(tidyverse)

##--------SCZ
# df 需要包含列：`590K_TRS`, `232K_TRS`
 df <- readr::read_csv("01_SCZ_replication.csv", show_col_types = FALSE)
 df <- readr::read_csv("01_SCZ_replication_all.csv", show_col_types = FALSE)
 
rho <- cor(df$`590K_TRS`, df$`232K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `590K_TRS`, y = `232K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 2, alpha = 0.35, color = "#BB8569") +
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
    title = sprintf("SCZ: 590K vs 232K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "590K TRS",
    y = "232K TRS"
  )

p_scatter
rho



##--------ASD
# df 需要包含列：`590K_TRS`, `232K_TRS`
df <- readr::read_csv("01_ASD_replication.csv", show_col_types = FALSE)
df <- readr::read_csv("01_ASD_replication_all.csv", show_col_types = FALSE)

rho <- cor(df$`590K_TRS`, df$`232K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `590K_TRS`, y = `232K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 2, alpha = 0.35, color = "#BB8569") +
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
    title = sprintf("ASD: 590K vs 232K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "590K TRS",
    y = "232K TRS"
  )

p_scatter
rho




##--------BIP
# df 需要包含列：`590K_TRS`, `232K_TRS`
df <- readr::read_csv("01_BIP_replication.csv", show_col_types = FALSE)
df <- readr::read_csv("01_BIP_replication_all.csv", show_col_types = FALSE)

rho <- cor(df$`590K_TRS`, df$`232K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `590K_TRS`, y = `232K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 2, alpha = 0.35, color = "#BB8569") +
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
    title = sprintf("BIP: 590K vs 232K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "590K TRS",
    y = "232K TRS"
  )

p_scatter
rho







##--------ADHD
# df 需要包含列：`590K_TRS`, `232K_TRS`
df <- readr::read_csv("01_ADHD_replication.csv", show_col_types = FALSE)
df <- readr::read_csv("01_ADHD_replication_all.csv", show_col_types = FALSE)

rho <- cor(df$`590K_TRS`, df$`232K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `590K_TRS`, y = `232K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 2, alpha = 0.35, color = "#BB8569") +
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
    title = sprintf("ADHD: 590K vs 232K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "590K TRS",
    y = "232K TRS"
  )

p_scatter
rho



##--------AN
# df 需要包含列：`590K_TRS`, `232K_TRS`
df <- readr::read_csv("01_AN_replication.csv", show_col_types = FALSE)
df <- readr::read_csv("01_AN_replication_all.csv", show_col_types = FALSE)

rho <- cor(df$`590K_TRS`, df$`232K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `590K_TRS`, y = `232K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 2, alpha = 0.35, color = "#BB8569") +
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
    title = sprintf("AN: 590K vs 232K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "590K TRS",
    y = "232K TRS"
  )

p_scatter
rho




##--------OCD
# df 需要包含列：`590K_TRS`, `232K_TRS`
df <- readr::read_csv("01_OCD_replication.csv", show_col_types = FALSE)
df <- readr::read_csv("01_OCD_replication_all.csv", show_col_types = FALSE)

rho <- cor(df$`590K_TRS`, df$`232K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `590K_TRS`, y = `232K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 2, alpha = 0.35, color = "#BB8569") +
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
    title = sprintf("OCD: 590K vs 232K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "590K TRS",
    y = "232K TRS"
  )

p_scatter
rho





##--------TS
# df 需要包含列：`590K_TRS`, `232K_TRS`
df <- readr::read_csv("01_TS_replication.csv", show_col_types = FALSE)
df <- readr::read_csv("01_TS_replication_all.csv", show_col_types = FALSE)

rho <- cor(df$`590K_TRS`, df$`232K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `590K_TRS`, y = `232K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 2, alpha = 0.35, color = "#BB8569") +
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
    title = sprintf("TS: 590K vs 232K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "590K TRS",
    y = "232K TRS"
  )

p_scatter
rho




##--------MDD
# df 需要包含列：`590K_TRS`, `232K_TRS`
df <- readr::read_csv("01_MDD_replication.csv", show_col_types = FALSE)
df <- readr::read_csv("01_MDD_replication_all.csv", show_col_types = FALSE)

rho <- cor(df$`590K_TRS`, df$`232K_TRS`, method="spearman", use="pairwise.complete.obs")

p_scatter <- ggplot(df, aes(x = `590K_TRS`, y = `232K_TRS`)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  geom_point(size = 2, alpha = 0.35, color = "#BB8569") +
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
    title = sprintf("MDD: 590K vs 232K TRS (Spearman \u03c1 = %.2f)", rho),
    x = "590K TRS",
    y = "232K TRS"
  )

p_scatter
rho

