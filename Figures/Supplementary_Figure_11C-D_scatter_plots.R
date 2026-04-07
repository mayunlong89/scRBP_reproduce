#----2025-10-15
setwd("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/02_scRBP_trs_10blood_cell_traits_100times_v2/00_4celltypes_results")

###--------------------------------mix two cell types-----------
#6.5K
data_2celltype <- read.table("Figure_2_monocyte_results_4regions_6.5K_CD14_monocytes.txt",header = TRUE,sep = "\t")
data_2celltype <- read.table("Figure_2_monocyte_results_4regions_6.5K_FCGR3A_monocytes.txt",header = TRUE,sep = "\t")
data_2celltype <- read.table("Figure_2_monocyte_results_4regions_6.5K_CD4_T_cells.txt",header = TRUE,sep = "\t")
data_2celltype <- read.table("Figure_2_monocyte_results_4regions_6.5K_CD8_T_cells.txt",header = TRUE,sep = "\t")

#8.7K
data_2celltype <- read.table("Figure_2_monocyte_results_4regions_8.7K_CD14_monocytes.txt",header = TRUE,sep = "\t")
data_2celltype <- read.table("Figure_2_monocyte_results_4regions_8.7K_FCGR3A_monocytes.txt",header = TRUE,sep = "\t")
data_2celltype <- read.table("Figure_2_monocyte_results_4regions_8.7K_CD4_T_cells.txt",header = TRUE,sep = "\t")
data_2celltype <- read.table("Figure_2_monocyte_results_4regions_8.7K_CD8_T_cells.txt",header = TRUE,sep = "\t")



# Plot using ggplot2
library(ggplot2)
library(ggplot2)
set.seed(42)

eps <- 1e-300
df <- transform(
  data_2celltype,
  x = -log10(pmax(q_RAS, eps)),
  y = -log10(pmax(q_RGS, eps)),
  z = -log10(pmax(q_TRS, eps))   # Also compute -log10(q_TRS)
)

# Also require z to be finite
df <- subset(df, is.finite(x) & is.finite(y) & is.finite(z))

threshold <- -log10(0.05)

# Label as Significant only when all three thresholds are met simultaneously
df$highlight <- ifelse(
  df$x >= threshold & df$y >= threshold & df$z >= threshold,
  "Significant", "Non-significant"
)
# Slight jitter
jit_w <- diff(range(df$x)) * 0.05
jit_h <- diff(range(df$y)) * 0.05
pos_jit <- position_jitter(width = jit_w, height = jit_h)

p <- ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(color = highlight), size = 2, alpha = 0.7, position = pos_jit) +
  stat_ellipse(level = 0.95, color = "lightblue", size = 0.5, type = "norm", na.rm = TRUE) +
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  labs(title = "data_2celltype", x = "RAS", y = "RGS", color = "Group") +
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "gray")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "top")

# Only draw the red ellipse when the subset is large enough and has variance
sig <- subset(df, highlight == "Significant")
if (nrow(sig) >= 3 && sd(sig$x) > 0 && sd(sig$y) > 0) {
  p <- p + stat_ellipse(data = sig, level = 0.95, color = "red", size = 0.5,
                        type = "norm", na.rm = TRUE)
}

p


###_-----final_version_dotplot
library(ggrepel)   # Added for label repelling

# ...(your previous code remains unchanged)...

p <- ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(color = highlight), size = 2, alpha = 0.7, position = pos_jit) +
  stat_ellipse(level = 0.95, color = "lightblue", size = 0.5, type = "norm", na.rm = TRUE) +
  geom_vline(xintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 1) +
  labs(title = "data_2celltype", x = "RAS", y = "RGS", color = "Group") +
  scale_color_manual(values = c("Significant" = "red", "Non-significant" = "gray")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "top")

# Add labels for significant points (in red)
p <- p + geom_text_repel(
  data = sig, aes(label = Regulon),              # The column in your table is named "Regulon"
  color = "red", size = 3,
  position = pos_jit,                            # Same jitter as the scatter points
  max.overlaps = Inf, box.padding = 0.2, point.padding = 0.2,
  min.segment.length = 0
)

# Optionally add a red ellipse for significant points
if (nrow(sig) >= 3 && sd(sig$x) > 0 && sd(sig$y) > 0) {
  p <- p + stat_ellipse(data = sig, level = 0.95, color = "red", size = 0.5,
                        type = "norm", na.rm = TRUE)
}

p
