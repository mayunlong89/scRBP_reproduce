# --- 2026-02-12
# LIN28B-let-7 target overlap enrichment analysis

# ============================================================
# 1. Load required packages
# ============================================================
library(VennDiagram)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# ============================================================
# 2. Hypergeometric test for target overlap significance
# ============================================================
# Parameter setup
total_genes    <- 20000   # Background gene universe size
let7_targets   <- 1208    # Number of let-7 target genes
lin28b_targets <- 55      # Number of LIN28B regulon target genes
overlap        <- 8       # Observed overlap count

# Hypergeometric test (right-tailed)
p_value <- phyper(q = overlap - 1,
                  m = let7_targets,
                  n = total_genes - let7_targets,
                  k = lin28b_targets,
                  lower.tail = FALSE)

# Expected overlap under the null (random draw)
expected <- (lin28b_targets * let7_targets) / total_genes

# Fold enrichment over expectation
fold_enrichment <- overlap / expected

# Print results
cat("=== Enrichment analysis results ===\n")
cat(sprintf("Observed overlap: %d\n", overlap))
cat(sprintf("Expected overlap: %.2f\n", expected))
cat(sprintf("Fold enrichment: %.2fx\n", fold_enrichment))
cat(sprintf("P-value: %.4f\n", p_value))
cat(sprintf("Significance: %s\n", ifelse(p_value < 0.05, "Significant", "Not significant")))

# ============================================================
# 3. Sensitivity analysis across different background sizes
# ============================================================
background_range <- seq(15000, 25000, 1000)
sensitivity_results <- data.frame(
  background = background_range,
  p_value = sapply(background_range, function(bg) {
    phyper(overlap - 1, let7_targets, bg - let7_targets,
           lin28b_targets, lower.tail = FALSE)
  }),
  fold = sapply(background_range, function(bg) {
    overlap / ((lin28b_targets * let7_targets) / bg)
  })
)

print(sensitivity_results)

# Plot sensitivity analysis
ggplot(sensitivity_results, aes(x = background, y = p_value)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(color = "steelblue", size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Sensitivity Analysis: Background Gene Number",
       x = "Background Gene Number",
       y = "P-value") +
  theme_minimal() +
  annotate("text", x = 20000, y = 0.052, label = "p = 0.05", color = "red")

ggsave("sensitivity_analysis.pdf", width = 6, height = 4)

# ============================================================
# 4. Venn diagram visualization
# ============================================================
venn.plot <- venn.diagram(
  x = list(
    "LIN28B targets\n(n=55)" = 1:55,
    "let-7 targets\n(n=1208)" = c(1:8, 56:1260)  # Simulated index sets
  ),
  filename = NULL,
  fill = c("lightblue", "pink"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20),
  cat.dist = 0.05,
  main = "LIN28B vs let-7 Target Overlap",
  main.cex = 1.5
)

# Save Venn diagram
pdf("venn_diagram.pdf", width = 8, height = 6)
grid.draw(venn.plot)
grid.text(paste0("p = ", format(p_value, digits = 3),
                 "\nFold enrichment = ", format(fold_enrichment, digits = 3)),
          x = 0.5, y = 0.15, gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

# ============================================================
# 5. Enrichment bar plot (observed vs expected)
# ============================================================
enrichment_data <- data.frame(
  Category = c("Observed", "Expected"),
  Count = c(overlap, expected),
  Type = c("Observed", "Expected")
)

p <- ggplot(enrichment_data, aes(x = Category, y = Count, fill = Type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f", Count)), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Observed" = "steelblue",
                               "Expected" = "gray70")) +
  labs(title = "let-7 Target Enrichment in LIN28B Network",
       subtitle = sprintf("p = %.4f, Fold = %.2fx", p_value, fold_enrichment),
       x = "", y = "Number of Genes") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12))

print(p)
ggsave("enrichment_barplot.pdf", width = 6, height = 5)

# ============================================================
# 6. Placeholder for actual gene list analysis (if available)
# ============================================================
# If you have the actual gene lists, replace the simulated data above:
# lin28b_genes  <- c("GENE1", "GENE2", ..., "GENE55")
# let7_genes    <- c("TARGET1", "TARGET2", ..., "TARGET1208")
# overlap_genes <- intersect(lin28b_genes, let7_genes)
# cat("\n=== Overlapping genes ===\n")
# print(overlap_genes)

# ============================================================
# 7. Summary table
# ============================================================
summary_table <- data.frame(
  Metric = c("Total LIN28B targets",
             "Total let-7 targets",
             "Observed overlap",
             "Expected overlap",
             "Fold enrichment",
             "P-value",
             "Significance"),
  Value = c(lin28b_targets,
            let7_targets,
            overlap,
            sprintf("%.2f", expected),
            sprintf("%.2f", fold_enrichment),
            sprintf("%.4f", p_value),
            ifelse(p_value < 0.05, "Yes (p<0.05)", "No"))
)

print(summary_table)
write.csv(summary_table, "enrichment_summary.csv", row.names = FALSE)

# ============================================================
# 8. Pie chart: let-7 dependency among LIN28B targets
# ============================================================
pie_data <- data.frame(
  Category = c("let-7 dependent\n(14.5%)",
               "let-7 independent\n(85.5%)"),
  Count = c(8, 47),
  Percentage = c(14.5, 85.5)
)

ggplot(pie_data, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  scale_fill_manual(values = c("let-7 dependent\n(14.5%)" = "#FF6B6B",
                               "let-7 independent\n(85.5%)" = "#4ECDC4")) +
  labs(title = "LIN28B Targets: let-7 Dependency") +
  theme_void() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_blank())

ggsave("let7_dependency_pie.pdf", width = 6, height = 5)

# ============================================================
# 9. Binomial test with 95% confidence interval (optional)
# ============================================================
binom_test <- binom.test(overlap, lin28b_targets,
                         p = let7_targets / total_genes)
cat("\n=== Binomial test ===\n")
print(binom_test)

# ============================================================
# 10. Save R workspace
# ============================================================
save.image("lin28b_let7_enrichment_analysis.RData")
cat("\nAll analysis results saved.\n")
