# Author: Wang Linglong
# Public release script derived from: 3\气泡图.txt
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
library(showtext)
library(paletteer)   # 提供高质量调色板
library(grid)

# === 字体设置 ===
showtext_auto()
font_add("Arial", "arial.ttf")
base_family <- if ("Arial" %in% sysfonts::font_families()) "Arial" else "sans"

# === 数据读取 ===
df <- read_csv(public_file("跨物种.csv")) %>%
  rename(host = `Host.species`,
         virus_family = `Virus.taxo(order:family)`,
         RPM = RPM) %>%
  mutate(RPM = as.numeric(RPM)) %>%
  filter(!is.na(host), !is.na(virus_family), !is.na(RPM)) %>%
  mutate(RPM_log = log10(RPM + 1))

# === 只取前一半病毒科，其他归为 Others ===
family_count <- df %>%
  count(virus_family, sort = TRUE)

top_half <- family_count$virus_family[1:ceiling(nrow(family_count)/2)]

df <- df %>%
  mutate(virus_family_plot = ifelse(virus_family %in% top_half,
                                    virus_family, "Others"))

# === Y轴物种排序（数量多的在上） ===
host_order <- df %>%
  group_by(host) %>%
  summarise(total_RPM = sum(RPM_log)) %>%
  arrange(desc(total_RPM)) %>%
  pull(host)

df$host <- factor(df$host, levels = rev(host_order))

# === 配色方案（离散调色板，不用渐变） ===
virus_levels <- unique(df$virus_family_plot)
n_colors <- length(virus_levels)
family_colors <- paletteer_d("nord::aurora")[1:n_colors]  # 高对比色系
names(family_colors) <- virus_levels
df$virus_family_plot <- factor(df$virus_family_plot, levels = virus_levels)

# === 绘图 ===
p <- ggplot(df, aes(x = virus_family_plot, y = host,
                    size = RPM_log, fill = virus_family_plot)) +
  geom_point(shape = 21, color = "black", alpha = 0.9, stroke = 0.3) +
  scale_size_continuous(range = c(4, 18), name = "log10(RPM + 1)") +
  scale_fill_manual(values = family_colors, name = "Virus family") +
  labs(x = "", y = "") +
  guides(
    fill = guide_legend(
      override.aes = list(size = 6),    # 图例圆圈大小
      ncol = 2                          # 图例分两列，避免太长
    ),
    size = guide_legend(
      override.aes = list(fill = "grey70", shape = 21, color = "black")
    )
  ) +
  theme_classic(base_family = base_family, base_size = 26) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               size = 26, color = "black"),
    axis.text.y = element_text(size = 26, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 22),
    legend.key.size = unit(1.3, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    plot.margin = margin(20, 20, 20, 10),
    panel.grid = element_line(color = "grey90", size = 0.2)
  )

# === 输出高分辨率图 ===
ggsave(public_output("Host_Virus_BubblePlot_NatureStyle.pdf"),
       plot = p, width = 26, height = 16, units = "in", dpi = 600)
ggsave(public_output("Host_Virus_BubblePlot_NatureStyle.png"),
       plot = p, width = 26, height = 16, units = "in", dpi = 600)


