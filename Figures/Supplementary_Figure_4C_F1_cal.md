# Cell Population Scalability Benchmarking Pipeline

> **Date:** 2025  
> **Author:** Yunlong Ma

---

## Overview

This pipeline benchmarks **memory consumption** and **runtime** of scRBP methods across increasing cell population sizes (2K → 200K cells). It reads summary statistics from an Excel file, generates line plots and bar plots with error bars, and performs Wilcoxon rank-sum tests to assess statistical significance between methods.

### Input

| File | Sheet | Columns |
|------|-------|---------|
| `01_Summary.xlsx` | `memory` | Cell_count, Method, Memory (MB) |
| `01_Summary.xlsx` | `time` | Cell_count, Method, Time (seconds) |

### Cell count levels

`2K, 5K, 10K, 25K, 50K, 100K, 200K`

---

## Setup

```r
library(tidyverse)
library(ggplot2)
library(readxl)

setwd("/Users/mayunlong/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/01_cell_population_subset_analysis/")

# Ordered factor levels for cell counts
CELL_LEVELS <- c("2K", "5K", "10K", "25K", "50K", "100K", "200K")
```

---

## Step 1 — Helper functions

To avoid repeating the same plotting code 4+ times, we define reusable functions for the line plot, bar plot, and Wilcoxon test.

```r
# ------------------------------------------------------------------
# Function: read and prepare data from Excel
# ------------------------------------------------------------------
read_benchmark_data <- function(file, sheet, value_col, scale_factor = 1) {
  df <- read_excel(file, sheet = sheet)
  df[[value_col]] <- df[[value_col]] * scale_factor
  df$Cell_count <- factor(df$Cell_count, levels = CELL_LEVELS)
  return(df)
}

# ------------------------------------------------------------------
# Function: line plot with error bars
# ------------------------------------------------------------------
plot_line <- function(data, value_col, y_label, title) {
  summary_df <- data %>%
    group_by(Cell_count, Method) %>%
    summarise(mean_val = mean(.data[[value_col]]),
              sd_val   = sd(.data[[value_col]]),
              .groups = "drop")

  ggplot(summary_df, aes(x = Cell_count, y = mean_val, color = Method, group = Method)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val), width = 0.2) +
    labs(title = title, x = "Cell Count", y = y_label) +
    theme_bw() +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# ------------------------------------------------------------------
# Function: grouped bar plot with error bars
# ------------------------------------------------------------------
plot_bar <- function(data, value_col, y_label, title) {
  summary_df <- data %>%
    group_by(Cell_count, Method) %>%
    summarise(mean_val = mean(.data[[value_col]]),
              sd_val   = sd(.data[[value_col]]),
              .groups = "drop")

  ggplot(summary_df, aes(x = Cell_count, y = mean_val, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
    geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                  position = position_dodge(width = 0.9), width = 0.25) +
    labs(title = title, x = "Cell Count", y = y_label) +
    theme_bw() +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# ------------------------------------------------------------------
# Function: Wilcoxon rank-sum test per Cell_count level
# ------------------------------------------------------------------
run_wilcoxon <- function(data, value_col) {
  data %>%
    group_by(Cell_count) %>%
    filter(n_distinct(Method) == 2) %>%
    summarise(
      p_value = wilcox.test(as.formula(paste(value_col, "~ Method")))$p.value,
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_value, method = "fdr"))
}
```

---

## Step 2 — Memory benchmarking

```r
# Read memory data (MB → GB)
mem_data <- read_benchmark_data("01_Summary.xlsx", "memory", "Memory", scale_factor = 1/1024)

# Line plot
plot_line(mem_data, "Memory", "Memory (GB)", "Memory Consumption across Cell Counts")

# Bar plot
plot_bar(mem_data, "Memory", "Memory (GB)", "Memory Consumption across Cell Counts")

# Wilcoxon test
mem_wilcox <- run_wilcoxon(mem_data, "Memory")
print(mem_wilcox)
```

---

## Step 3 — Runtime benchmarking

```r
# Read time data (seconds → hours)
time_data <- read_benchmark_data("01_Summary.xlsx", "time", "Time", scale_factor = 1/3600)

# Line plot
plot_line(time_data, "Time", "Time (Hours)", "Time Consumption across Cell Counts")

# Bar plot
plot_bar(time_data, "Time", "Time (Hours)", "Time Consumption across Cell Counts")

# Wilcoxon test
time_wilcox <- run_wilcoxon(time_data, "Time")
print(time_wilcox)
```

---

## Summary of outputs

| Output | Plot type | Y-axis unit | Description |
|--------|-----------|-------------|-------------|
| Memory line plot | Line + error bars | GB | Memory trend across cell counts per method |
| Memory bar plot | Grouped bar + error bars | GB | Side-by-side memory comparison |
| Time line plot | Line + error bars | Hours | Runtime trend across cell counts per method |
| Time bar plot | Grouped bar + error bars | Hours | Side-by-side runtime comparison |
| Wilcoxon test (memory) | Table | — | Per-cell-count p-values (FDR-corrected) |
| Wilcoxon test (time) | Table | — | Per-cell-count p-values (FDR-corrected) |

## Key simplification

The original script contained **8 near-identical code blocks** (~180 lines) that differed only in the metric (Memory vs. Time), unit conversion factor, and plot type (line vs. bar). This refactored version uses **3 reusable functions** (`plot_line`, `plot_bar`, `run_wilcoxon`), reducing the code to ~60 lines while producing identical outputs.
