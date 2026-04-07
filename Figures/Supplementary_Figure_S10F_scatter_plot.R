#----2025-10-13
setwd("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/02_scRBP_trs_10blood_cell_traits_100times_v2")
###--------------------------------mix two cell types-----------
#6.6K
#data_2celltype <- read.table("Figure_2_monocyte_results_4regions_6.6K_monocytes.txt",header = TRUE,sep = "\t")
#data_2celltype <- read.table("Figure_2_monocyte_results_4regions_6.6K_T_cells.txt",header = TRUE,sep = "\t")
#8.8K
data_2celltype <- read.table("Figure_2_monocyte_results_4regions_8.8K_monocytes.txt",header = TRUE,sep = "\t")
data_2celltype <- read.table("Figure_2_monocyte_results_4regions_8.8K_T_cells.txt",header = TRUE,sep = "\t")
# Plot using ggplot2
library(ggplot2)
library(ggplot2)
set.seed(42)
eps <- 1e-300
df <- transform(
  data_2celltype,
  x = -log10(pmax(q_RAS, eps)),
  y = -log10(pmax(q_RGS, eps)),
  z = -log10(pmax(q_TRS, eps))   # also compute -log10(q_TRS)
)
# also require z to be finite
df <- subset(df, is.finite(x) & is.finite(y) & is.finite(z))
threshold <- -log10(0.1)
# labeled as Significant only when all three thresholds are met simultaneously
df$highlight <- ifelse(
  df$x >= threshold & df$y >= threshold & df$z >= threshold,
  "Significant", "Non-significant"
)
# slight jitter
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
# only draw the red ellipse if the subset is large enough and has variance
sig <- subset(df, highlight == "Significant")
if (nrow(sig) >= 3 && sd(sig$x) > 0 && sd(sig$y) > 0) {
  p <- p + stat_ellipse(data = sig, level = 0.95, color = "red", size = 0.5,
                        type = "norm", na.rm = TRUE)
}
p
###_-----final_version_dotplot
library(ggrepel)   # new addition
# ... (keep your preceding code unchanged) ...
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
# add labels to significant points (in red)
p <- p + geom_text_repel(
  data = sig, aes(label = Regulon),              # the column in your table is named "Regulon"
  color = "red", size = 3,
  position = pos_jit,                            # same jitter as the scatter points
  max.overlaps = Inf, box.padding = 0.2, point.padding = 0.2,
  min.segment.length = 0
)
# optionally add a red ellipse around significant points
if (nrow(sig) >= 3 && sd(sig$x) > 0 && sd(sig$y) > 0) {
  p <- p + stat_ellipse(data = sig, level = 0.95, color = "red", size = 0.5,
                        type = "norm", na.rm = TRUE)
}
p
