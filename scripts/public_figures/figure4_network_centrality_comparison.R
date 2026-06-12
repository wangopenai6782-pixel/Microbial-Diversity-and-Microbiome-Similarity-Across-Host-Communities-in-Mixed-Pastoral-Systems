# Author: Wang Linglong
# Public release script derived from: 4\多网络对比.txt
# Purpose: reproduce the corresponding analysis/figure in the manuscript.
#
# Input files should be placed under data/public_figures/.
# Override input/output locations if needed:
#   Sys.setenv(JEV_FIGURE_INPUT_DIR = "/path/to/input")
#   Sys.setenv(JEV_FIGURE_OUTPUT_DIR = "/path/to/output")

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)][1])
script_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path, mustWork = FALSE)) else getwd()
release_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
input_dir <- Sys.getenv("JEV_FIGURE_INPUT_DIR", file.path(release_root, "data", "public_figures"))
output_dir <- Sys.getenv("JEV_FIGURE_OUTPUT_DIR", file.path(release_root, "results", "public_figures"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
public_file <- function(...) file.path(input_dir, ...)
public_output <- function(...) file.path(output_dir, ...)
library(tidyverse)
library(cowplot)
library(showtext)
library(grid)

# === 一键设置全局字体大小函数 ===
set_global_font <- function(size = 14, family = "Arial") {
  theme_set(theme_minimal(base_size = size, base_family = family))
}

# === 设置全局字体 ===
showtext_auto()
font_add("Arial", "arial.ttf")
set_global_font(size = 17, family = "Arial")

# === 数据读取 ===
df <- read.csv("all_microbe_centrality.csv") %>%
  filter(microbe_type %in% c("Protozoa", "Virus", "Bacteria", "Fungi", "Total"))

df$microbe_type <- factor(df$microbe_type, levels = rev(c("Protozoa", "Virus", "Bacteria", "Fungi", "Total")))

# === 颜色方案 ===
colors <- c("#4C72B0", "#55A868", "#C44E52", "#8172B3", "#CCB974")

# === 输出 P 值结果到控制台 ===
p_degree_test <- kruskal.test(degree ~ microbe_type, data = df)
p_betw_test <- kruskal.test(abs(betweenness) ~ microbe_type, data = df)

cat("Degree Kruskal-Wallis test:\n")
print(p_degree_test)
cat("\nBetweenness Kruskal-Wallis test:\n")
print(p_betw_test)

# === Degree 图 ===
df_degree <- df %>% select(microbe_type, degree)

p_degree <- ggplot(df_degree, aes(x = degree, y = microbe_type, fill = microbe_type)) +
  geom_boxplot(alpha = 0.6, color = "gray30", size = 0.9, width = 0.55, outlier.shape = NA) +
  geom_jitter(color = "gray30", alpha = 0.4, size = 2, height = 0.15) +
  scale_fill_manual(values = colors) +
  labs(x = "Degree", y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "gray40", size = 0.8),
    axis.text.x = element_text(size = 12),
    legend.position = "none",
    plot.margin = margin(5,5,5,5)
  )

# === Betweenness 图 ===
df_betw <- df %>% select(microbe_type, betweenness) %>% mutate(betweenness = abs(betweenness))

p_betw <- ggplot(df_betw, aes(x = betweenness, y = microbe_type, fill = microbe_type)) +
  geom_boxplot(alpha = 0.6, color = "gray30", size = 0.9, width = 0.55, outlier.shape = NA) +
  geom_jitter(color = "gray30", alpha = 0.4, size = 2, height = 0.15) +
  scale_x_reverse() +
  scale_fill_manual(values = colors) +
  labs(x = "Betweenness", y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "gray40", size = 0.8),
    axis.text.x = element_text(size = 12),
    legend.position = "none",
    plot.margin = margin(5,5,5,5)
  )

# === 中间标签 ===
network_labels <- levels(df$microbe_type)

middle_labels <- ggplot() +
  geom_text(aes(x = 1, y = seq_along(network_labels)), label = network_labels,
            size = 6, fontface = "bold", family = "Arial") +
  scale_y_reverse(expand = c(0.05,0)) +
  theme_void() +
  coord_cartesian(ylim = c(0.5, length(network_labels) + 0.5)) +
  theme(plot.margin = margin(5,5,5,5))

# === 拼接图形 ===
final_plot <- plot_grid(
  p_betw,
  middle_labels,
  p_degree,
  ncol = 3,
  rel_widths = c(1, 0.3, 1),
  align = "h"
)

# === 保存高分辨率 PDF ===
ggsave(public_output("nature_style_degree_betweenness_mirror_with_jitter.pdf"),
       final_plot, width = 10, height = 6, dpi = 600)


