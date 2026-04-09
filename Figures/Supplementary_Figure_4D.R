# Load required libraries
library(tidyverse)
library(ggplot2)
library(readxl)

setwd("/Users/mayunlong/Desktop/UPENN/00UPenn-Projects/02-Aimed_projects/05-devBrain-RBP-methods/00-Manuscript-scRBP-2024/01_Data_analysis/01_cell_population_subset_analysis/")

# ============================================================
# 1) Memory consumption — line plot (MB)
# ============================================================

# Read data
memory_data <- read_excel("01_Summary.xlsx", sheet = "memory")

# Convert Cell_count to an ordered factor for proper plot ordering
memory_data$Cell_count <- factor(memory_data$Cell_count, levels = c("2K","5K","10K","25K","50K","100K","200K"))

# Summarize: compute mean and standard deviation
summary_data <- memory_data %>%
  group_by(Cell_count, Method) %>%
  summarise(mean_memory = mean(Memory),
            sd_memory = sd(Memory))

# Line plot with error bars
memory_plot <- ggplot(summary_data, aes(x = Cell_count, y = mean_memory, color = Method, group = Method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_memory - sd_memory, ymax = mean_memory + sd_memory), width = 0.2) +
  labs(title = "Memory Consumption across Cell Counts",
       x = "Cell Count",
       y = "Memory (MB)") +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

memory_plot


# ============================================================
# 2) Time consumption — line plot (seconds)
# ============================================================

time_data <- read_excel("01_Summary.xlsx", sheet = "time")

# Convert Cell_count to an ordered factor
time_data$Cell_count <- factor(time_data$Cell_count, levels = c("2K","5K","10K","25K","50K","100K","200K"))

# Summarize: compute mean and standard deviation
summary_data <- time_data %>%
  group_by(Cell_count, Method) %>%
  summarise(mean_time = mean(Time),
            sd_time = sd(Time))

# Line plot with error bars
time_plot <- ggplot(summary_data, aes(x = Cell_count, y = mean_time, color = Method, group = Method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_time - sd_time, ymax = mean_time + sd_time), width = 0.2) +
  labs(title = "Time Consumption across Cell Counts",
       x = "Cell Count",
       y = "Time (Seconds)") +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

time_plot


# ============================================================
# 3) Time consumption — line plot (converted to hours)
# ============================================================

time_data <- read_excel("01_Summary.xlsx", sheet = "time")

# Convert seconds to hours
time_data <- time_data %>% mutate(Time = Time / 3600)

# Convert Cell_count to an ordered factor
time_data$Cell_count <- factor(time_data$Cell_count, levels = c("2K","5K","10K","25K","50K","100K","200K"))

# Summarize: compute mean and standard deviation
summary_data <- time_data %>%
  group_by(Cell_count, Method) %>%
  summarise(mean_time = mean(Time),
            sd_time = sd(Time))

# Line plot with error bars
time_plot <- ggplot(summary_data, aes(x = Cell_count, y = mean_time, color = Method, group = Method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_time - sd_time, ymax = mean_time + sd_time), width = 0.2) +
  labs(title = "Time Consumption across Cell Counts",
       x = "Cell Count",
       y = "Time (Hours)") +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

time_plot


# ============================================================
# 4) Memory consumption — line plot (converted to GB)
# ============================================================

memory_data <- read_excel("01_Summary.xlsx", sheet = "memory")

# Convert MB to GB
memory_data <- memory_data %>% mutate(Memory = Memory / 1024)

# Convert Cell_count to an ordered factor
memory_data$Cell_count <- factor(memory_data$Cell_count, levels = c("2K","5K","10K","25K","50K","100K","200K"))

# Summarize: compute mean and standard deviation
summary_data <- memory_data %>%
  group_by(Cell_count, Method) %>%
  summarise(mean_memory = mean(Memory),
            sd_memory = sd(Memory))

# Line plot with error bars
memory_plot <- ggplot(summary_data, aes(x = Cell_count, y = mean_memory, color = Method, group = Method)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_memory - sd_memory, ymax = mean_memory + sd_memory), width = 0.2) +
  labs(title = "Memory Consumption across Cell Counts",
       x = "Cell Count",
       y = "Memory (GB)") +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

memory_plot


# ============================================================
# 5) Memory consumption — bar plot (GB) with error bars
# ============================================================

memory_data <- read_excel("01_Summary.xlsx", sheet = "memory")

# Convert MB to GB
memory_data <- memory_data %>%
  mutate(Memory = Memory / 1024)

# Convert Cell_count to an ordered factor
memory_data$Cell_count <- factor(memory_data$Cell_count, levels = c("2K","5K","10K","25K","50K","100K","200K"))

# Summarize: compute mean and standard deviation
summary_data <- memory_data %>%
  group_by(Cell_count, Method) %>%
  summarise(mean_memory = mean(Memory),
            sd_memory = sd(Memory),
            .groups = "drop")

# Grouped bar plot with error bars
memory_barplot <- ggplot(summary_data, aes(x = Cell_count, y = mean_memory, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  geom_errorbar(aes(ymin = mean_memory - sd_memory, ymax = mean_memory + sd_memory),
                position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Memory Consumption across Cell Counts",
       x = "Cell Count",
       y = "Memory (GB)") +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

memory_barplot


# ============================================================
# 6) Time consumption — bar plot (hours) with error bars
# ============================================================

time_data <- read_excel("01_Summary.xlsx", sheet = "time")

# Convert seconds to hours
time_data <- time_data %>%
  mutate(Time = Time / 3600)

# Convert Cell_count to an ordered factor
time_data$Cell_count <- factor(time_data$Cell_count, levels = c("2K","5K","10K","25K","50K","100K","200K"))

# Summarize: compute mean and standard deviation
summary_data <- time_data %>%
  group_by(Cell_count, Method) %>%
  summarise(mean_time = mean(Time),
            sd_time = sd(Time),
            .groups = "drop")

# Grouped bar plot with error bars
time_barplot <- ggplot(summary_data, aes(x = Cell_count, y = mean_time, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  geom_errorbar(aes(ymin = mean_time - sd_time, ymax = mean_time + sd_time),
                position = position_dodge(width = 0.9), width = 0.25) +
  labs(title = "Time Consumption across Cell Counts",
       x = "Cell Count",
       y = "Time (Hours)") +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

time_barplot


# ============================================================
# 7) Statistical test: Wilcoxon test for time (between methods)
# ============================================================

time_data <- read_excel("01_Summary.xlsx", sheet = "time") %>%
  mutate(Time = Time / 3600,
         Cell_count = factor(Cell_count, levels = c("2K", "5K", "10K", "25K", "50K", "100K", "200K")))

# Wilcoxon rank-sum test at each Cell_count level (requires exactly 2 methods)
wilcox_results <- time_data %>%
  group_by(Cell_count) %>%
  filter(n_distinct(Method) == 2) %>%
  summarise(
    p_value = wilcox.test(Time ~ Method)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))  # Optional: FDR correction

print(wilcox_results)


# ============================================================
# 8) Statistical test: Wilcoxon test for memory (between methods)
# ============================================================

memory_data <- read_excel("01_Summary.xlsx", sheet = "memory") %>%
  mutate(Memory = Memory / 1024,  # MB to GB
         Cell_count = factor(Cell_count, levels = c("2K", "5K", "10K", "25K", "50K", "100K", "200K")))

# Wilcoxon rank-sum test at each Cell_count level (requires exactly 2 methods)
memory_wilcox_results <- memory_data %>%
  group_by(Cell_count) %>%
  filter(n_distinct(Method) == 2) %>%
  summarise(
    p_value = wilcox.test(Memory ~ Method)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))  # Optional: FDR correction

print(memory_wilcox_results)
