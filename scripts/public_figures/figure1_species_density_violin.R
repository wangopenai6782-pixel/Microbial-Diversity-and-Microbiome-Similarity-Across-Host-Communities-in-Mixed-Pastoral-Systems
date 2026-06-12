# Author: Wang Linglong
# Public release script derived from: 1\小提琴.txt
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
# 📦 安装和加载必要包
packages <- c("ggplot2", "dplyr", "readr", "forcats")
new_packages <- packages[!packages %in% installed.packages()]
if(length(new_packages)) install.packages(new_packages)

library(ggplot2)
library(dplyr)
library(readr)
library(forcats)

# 📂 加载数据
data <- read_csv(public_file("species_data_with_richness.csv"))

# 🧹 数据预处理
data <- data %>%
  mutate(
    habitat = factor(habitat, levels = c("Pasture", "Natural grassland")),
    city = factor(city),
    species = fct_infreq(species)  # 按频率排序
  )

# ❌ 不再筛选 top12，直接保留所有物种（21 个）
plot_data <- data

# 🎨 设置Nature风格配色
habitat_colors <- c("Pasture" = "#D62728", "Natural grassland" = "#1F77B4")  # 红蓝对比
city_colors <- c(
  "Hohhot"   = "#636363",
  "Ordos"    = "#BDBDBD",
  "Baotou"   = "#969696",
  "Xilingol" = "#252525"
)

# ============ 全局字体大小设置 ============
global_font_size <- 26 # 一键控制全局字体大小

# 📈 绘图
p <- ggplot(plot_data, aes(x = habitat, y = species_density, fill = habitat)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9, color = "black", lwd = 0.3) +
  geom_jitter(aes(color = city), size = 0.8, width = 0.15, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.3, 
               fill = "white", color = "black", stroke = 0.3) +
  facet_wrap(~ species, ncol = 3, nrow = 7, scales = "free_y") +   # ✅ 每行 3 个，共 7 行
  scale_fill_manual(values = habitat_colors) +
  scale_color_manual(values = city_colors) +
  theme_minimal(base_size = global_font_size) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "non"
  ) +
  labs(
    y = "Species Density",
    fill = "Habitat",
    color = "City"
  )

# 💾 导出高分辨率图像（PDF 矢量图，更推荐用于出版）
ggsave(public_output("Nature_Style_Species_Violin.pdf"), 
       p, width = 14, height = 18, device = cairo_pdf)   # 调大尺寸以保证 21 图清晰


