###-----ranking top RBP with significant SNPs
#2025-12-30
###-----ranking top RBP with significant SNPs
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})
#=========================
# 0) Input
#=========================
df <- read.table(text='
RBP Count
DDX3X 0
FUBP3 0
KHSRP 0
RTCB 0
FMR1 0
DDX21 1
RPS11 1
CENPC 1
PRPF8 1
YBX3 1
SRSF1 1
ELAVL1 1
EIF3D 1
DICER1 1
NCBP3 1
DEK 1
SAMD4A 1
EIF4A1 1
VIM 2
NIPBL 2
TRA2A 2
RPS19 2
G3BP1 2
FAM120A 2
RBM47 2
G3BP2 2
TAF15 2
RPS3A 2
YWHAZ 2
NOP58 2
YTHDF3 3
AKAP1 3
SRSF5 3
SRSF2 3
HNRNPU 3
MTDH 3
XRN2 3
METAP2 3
RPS15A 4
IFIT2 4
SERBP1 4
TNRC6C 5
HNRNPC 5
RYBP 5
PABPC1 5
ZRANB2 5
OAS1 5
SRP14 5
TNRC6A 5
QKI 6
ZC3HAV1 6
HNRNPA1 6
RPRD2 6
IGF2BP2 7
CELF1 7
CPEB4 9
DDX6 9
CHD7 10
', header = TRUE, stringsAsFactors = FALSE)
#=========================
# 1) Rank / ordering
#=========================
dfp <- df %>%
  arrange(Count, RBP) %>%
  mutate(RBP = factor(RBP, levels = RBP))
# Color gradient (light to dark warm tones, similar to the example)
pal <- colorRampPalette(c("#FAF3E3", "#FBB394", "#F87C56"))(256)
# Select which points to label: e.g., those with Count >= 6 (adjust threshold as needed)
dfp <- dfp %>% mutate(label = ifelse(Count >= 6, as.character(RBP), NA))
#=========================
# 2) Plot (ranked step + gradient points)
#=========================
p <- ggplot(dfp, aes(x = RBP, y = Count)) +
  geom_step(aes(group = 1), linewidth = 0.9, color = "grey35") +
  geom_point(aes(fill = Count),
             shape = 21, color = "grey20", stroke = 0.25, size = 3.0) +
  scale_fill_gradientn(colours = pal, name = "Count") +
  ggrepel::geom_text_repel(
    aes(label = label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0
  ) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(
    x = NULL,
    y = "Number of traits with ≥1 significant SNP",
    title = "RBP pleiotropy across 10 traits (ranked by Count)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold")
  )
p
#=========================
# 3) Save
#=========================
ggsave("RBP_rank_count_gradient.pdf", p, width = 7, height = 6)
ggsave("RBP_rank_count_gradient.png", p, width = 7, height = 6, dpi = 300)
