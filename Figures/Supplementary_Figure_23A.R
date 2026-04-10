##----Boxplots of regulon activity across developmental stages

#regulon_id <- "LIN28B_CDS"
#out_dir <- "LIN28B_CDS_dev_boxplot"



#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# =========================
# Inputs
# =========================
anno_file <- "fb590k_donorStage_v3.mean.pb_anno.csv"
act_file  <- "fb590k_donorStage_v3.mean.activity_pb.devRegulons.csv"


regulon_id <- "LIN28B_CDS"
out_dir <- "LIN28B_CDS_dev_boxplot"



dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 1) Read annotation + activity
# =========================
anno <- fread(anno_file, data.table = FALSE)
stopifnot(all(c("pb_id","stage") %in% colnames(anno)))

act <- fread(act_file, data.table = FALSE, check.names = FALSE)
stopifnot("pb_id" %in% colnames(act))
if (!(regulon_id %in% colnames(act))) {
  stop("Regulon not found in activity file: ", regulon_id)
}

df <- anno %>%
  select(pb_id, stage, donor = any_of("donor"), n_cells = any_of("n_cells")) %>%
  inner_join(act %>% select(pb_id, !!regulon_id), by = "pb_id") %>%
  rename(activity = !!regulon_id) %>%
  mutate(activity = as.numeric(activity))

# =========================
# 2) Stage order (edit to your exact labels if needed)
# =========================
# If stage is numeric 1..9, keep numeric order; if it's string labels, set factor order.
if (is.numeric(df$stage) || all(grepl("^[0-9]+$", df$stage))) {
  df$stage <- as.integer(df$stage)
  df <- df %>% arrange(stage) %>% mutate(stage = factor(stage, levels = sort(unique(stage))))
} else {
  # Example order (replace with your real stage labels if they are strings)
  stage_levels <- unique(df$stage)
  df$stage <- factor(df$stage, levels = stage_levels)
}

# =========================
# 3) Plot: boxplot across stages
# =========================
p <- ggplot(df, aes(x = stage, y = activity, fill = stage)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5) +
  theme_classic(base_size = 11) +
  labs(
    x = "Developmental stage",
    y = "Regulon activity",
    title = paste0(regulon_id, " regulon activity across development")
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(out_dir, paste0(regulon_id, ".boxplot.pdf")),
       p, width = 3.2, height = 2.8, useDingbats = FALSE)
ggsave(file.path(out_dir, paste0(regulon_id, ".boxplot.png")),
       p, width = 3.2, height = 2.8, dpi = 300)

cat("DONE. Outputs in:", out_dir, "\n")

