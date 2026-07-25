#2026-01-28

setwd("~/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/13-scRBP_TRS_test/03_scRBP_trs_fetal_brain_15times/01_Final_scRBP_TRS_fetal_brain_50times")


suppressPackageStartupMessages({
  library(tidyverse)
})

# =========================
# 0) Input
# =========================
infile <- "01_heatmap_8braindiseases.csv"
outdir <- "01_heatmap_8braindiseases_Figure_7_heatmap_stratifiedby_disease_dotplot_celltype_x_regulon"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

df0 <- read.csv(infile, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot("ID" %in% colnames(df0))

# 自动识别 traits
trs_cols <- grep("_TRS$", colnames(df0), value = TRUE)
p_cols   <- grep("_P$",   colnames(df0), value = TRUE)
traits <- sub("_TRS$", "", trs_cols)
keep_traits <- traits[paste0(traits, "_P") %in% p_cols]
stopifnot(length(keep_traits) > 0)

preferred_trait_order <- c("SCZ","ASD","BIP","TS","OCD","MDD","AN","ADHD")
trait_order <- c(preferred_trait_order[preferred_trait_order %in% keep_traits],
                 setdiff(keep_traits, preferred_trait_order))

# =========================
# 1) Parse ID -> CellType + Regulon(=RBP_region)
#    ID 格式：<regulon>_<celltype>（最后一个 "_" 后是 cell type）
# =========================
id_vec <- df0$ID

cell_type <- sub(".*_", "", id_vec) %>% trimws()
regulon   <- sub("_[^_]+$", "", id_vec) %>% trimws()   # 去掉最后一个 _celltype

# 可选：统一 cell type 名称
map_celltype <- c(
  "Excitatory neuron" = "Excitatory Neuron",
  "Inhibitory neuron" = "Inhibitory Neuron",
  "Newborn neuron"    = "Newborn Neuron"
)
cell_type <- dplyr::recode(cell_type, !!!map_celltype, .default = cell_type)

# cell type 排序（按你常用生物顺序）
preferred_celltype_order <- c(
  "Excitatory Neuron","Inhibitory Neuron","Newborn Neuron",
  "IPC","Oligodendrocyte","Astrocyte","Microglia",
  "Pericyte","Endothelial",
  "RG","RG.div","Div","OPC","OPC.div"
)
cell_levels <- c(preferred_celltype_order[preferred_celltype_order %in% unique(cell_type)],
                 setdiff(unique(cell_type), preferred_celltype_order))
cell_type <- factor(cell_type, levels = cell_levels)

# =========================
# 2) helper: significance label
#    你要的符号：#/*/**
#    （你也可以把阈值改成 FDR 或别的）
# =========================
sig_symbol <- function(p) {
  case_when(
    is.na(p)        ~ "",
    p <= 0.001       ~ "**",
    p <= 0.01        ~ "*",
    p <= 0.05        ~ "#",
    TRUE            ~ ""
  )
}

# =========================
# 3) main: per-disease dotplot
# =========================
make_dotplot <- function(trait,
                         p_cut = 0.05,
                         topN_regulons = 60,     # 太多会挤爆；想更全就调大
                         size_by = c("absTRS", "neglog10P", "none")) {
  
  size_by <- match.arg(size_by)
  
  trs_col <- paste0(trait, "_TRS")
  p_col   <- paste0(trait, "_P")
  stopifnot(trs_col %in% colnames(df0), p_col %in% colnames(df0))
  
  df <- tibble(
    Regulon  = regulon,
    CellType = cell_type,
    TRS      = suppressWarnings(as.numeric(df0[[trs_col]])),
    P        = suppressWarnings(as.numeric(df0[[p_col]]))
  ) %>%
    filter(!is.na(Regulon), !is.na(CellType), !is.na(TRS), !is.na(P)) %>%
    mutate(
      TRS_sign = ifelse(TRS >= 0, "Positive TRS", "Negative TRS"),
      sig_lab  = sig_symbol(P),
      neglog10P = -log10(pmax(P, 1e-300)),
      absTRS = abs(TRS)
    )
  
  # 只保留“该病在任一 cell type 显著”的 regulons（更接近你截图：只展示有意义的那批）
  keep_reg <- df %>%
    group_by(Regulon) %>%
    summarise(minP = min(P, na.rm = TRUE), maxAbsTRS = max(absTRS, na.rm = TRUE), .groups = "drop") %>%
    filter(minP < p_cut) %>%
    arrange(minP) %>%
    slice_head(n = topN_regulons) %>%
    pull(Regulon)
  
  if (length(keep_reg) == 0) {
    message("[Skip] ", trait, ": no regulons pass P<", p_cut)
    return(NULL)
  }
  
  df_plot <- df %>%
    filter(Regulon %in% keep_reg) %>%
    mutate(
      Regulon = factor(Regulon, levels = keep_reg)  # x轴按显著性排序（minP从小到大）
    )
  
  # 点大小映射
  if (size_by == "absTRS") {
    df_plot <- df_plot %>% mutate(pt_size = scales::rescale(absTRS, to = c(1.5, 7)))
    size_lab <- "|TRS|"
  } else if (size_by == "neglog10P") {
    df_plot <- df_plot %>% mutate(pt_size = scales::rescale(neglog10P, to = c(1.5, 7)))
    size_lab <- "-log10(P)"
  } else {
    df_plot <- df_plot %>% mutate(pt_size = 4)
    size_lab <- NULL
  }
  
  p <- ggplot(df_plot, aes(x = Regulon, y = CellType)) +
    geom_point(aes(color = TRS_sign, size = pt_size), alpha = 0.85) +
    geom_text(aes(label = sig_lab), color = "black", size = 3) +
    scale_color_manual(values = c("Positive TRS" = "red", "Negative TRS" = "gray70")) +
    guides(size = "none") +
    labs(
      title = paste0(trait, ": Cell type × Regulon TRS"),
      subtitle = paste0("Color: TRS sign (red=positive, gray=negative); labels: #/*/** by P; shown regulons: min(P)<", p_cut),
      x = NULL,
      y = NULL,
      color = "TRS"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.major.y = element_line(color = "gray92")
    )
  
  p
}

# =========================
# 4) Run all diseases
# =========================
for (tr in trait_order) {
  p <- make_dotplot(
    trait = tr,
    p_cut = 0.05,
    topN_regulons = 60,     # 你如果想更像论文图（regulon少），可以 15~30
    size_by = "absTRS"      # 推荐 absTRS；也可改 "neglog10P"
  )
  
  if (is.null(p)) next
  
  png_out <- file.path(outdir, paste0(tr, "_CellType_x_Regulon_dotplot.png"))
  pdf_out <- file.path(outdir, paste0(tr, "_CellType_x_Regulon_dotplot.pdf"))
  
  ggsave(png_out, p, width = 13, height = 6.8, dpi = 300)
  ggsave(pdf_out, p, width = 13, height = 6.8)
  message("[Done] ", tr)
}

message("All done. Output in: ", outdir)


