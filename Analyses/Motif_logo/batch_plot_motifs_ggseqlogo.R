# batch_plot_motifs_ggseqlogo.R

suppressPackageStartupMessages({
  library(universalmotif)
  library(ggseqlogo)
  library(ggplot2)
})

# =========================
# 1. Input / output
# =========================
meme_file <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/11_616RBPs_20746motifs_all_databases.meme"
outdir <- "/mnt/isilon/gandal_lab/mayl/05_RNA_binding_protein/RBP_binding_sites/00_all_motifs_v2/all_20746motifs_logo_ggseqlogo_eps_2"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 2. Read motifs from MEME
# =========================
motifs <- read_meme(meme_file)
cat("Total motifs loaded:", length(motifs), "\n")

# =========================
# 3. Plot each motif
# =========================
for (i in seq_along(motifs)) {
  
  m <- motifs[[i]]
  
  motif_name <- m@name
  if (is.null(motif_name) || motif_name == "") {
    motif_name <- paste0("motif_", i)
  }
  
  # Convert to PPM
  m_ppm <- tryCatch(
    convert_type(m, "PPM"),
    error = function(e) {
      message("Skip ", motif_name, ": cannot convert to PPM")
      return(NULL)
    }
  )
  
  if (is.null(m_ppm)) next
  
  mat <- m_ppm@motif
  
  # Check row names
  if (!all(rownames(mat) %in% c("A", "C", "G", "T", "U"))) {
    message("Skip ", motif_name, ": unexpected rownames")
    next
  }
  
  # RNA style
  if ("T" %in% rownames(mat)) {
    rownames(mat)[rownames(mat) == "T"] <- "U"
    seq_type_use <- "rna"
  } else {
    seq_type_use <- "rna"
  }
  
  p <- ggseqlogo(
    mat,
    method = "bits",
    seq_type = seq_type_use
  ) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 0.8, colour = "black"),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 14, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 18, face = "bold")
    )
  
  ggsave(
    filename = file.path(outdir, paste0(motif_name, ".eps")),
    plot = p,
    width = 6.0,
    height = 3.0,
    units = "in",
    device = cairo_ps
  )
  
  if (i %% 100 == 0) {
    cat("Processed", i, "motifs\n")
  }
}

cat("Done.\n")
