# Author: Wang Linglong
# Public release script derived from: 2\微生物堆积图+Z-score 柱状图.txt
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
# 加载必要的包
library(ggplot2)
library(dplyr)
library(readr)

# 读取数据
data <- read_csv(public_file("各文库微生物多样性11.csv"))

# 计算 Z-score 标准化（对每种微生物类别进行标准化）
data_zscore <- data %>%
  group_by(Type) %>%
  mutate(Z_score = (Absolute_abundance - mean(Absolute_abundance, na.rm = TRUE)) / 
                     sd(Absolute_abundance, na.rm = TRUE)) %>%
  ungroup()

# 设定颜色主题（增强对比度）
colors <- c("Viruses" = "#D73027",   # 深红
            "Parasite" = "#4575B4",   # 深蓝
            "Fungi" = "#1A9850",  # 亮绿
            "Bacteria" = "#FDAE61")     # 亮橙

# 绘制 Z-score 柱状图
ggplot(data_zscore, aes(x = Sample, y = Z_score, fill = Type)) +
  geom_col(position = "dodge", width = 0.7, color = "black", alpha = 0.9) +
  scale_fill_manual(values = colors) +
  theme_minimal(base_size = 14) +
  labs(title = "标准化微生物丰度 (Z-score)", 
       x = "样本编号", 
       y = "标准化丰度 (Z-score)", 
       fill = "微生物类别") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),  
        axis.text.y = element_text(size = 12, face = "bold", color = "black"),  
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  
        legend.title = element_text(size = 14, face = "bold"),  
        legend.text = element_text(size = 12)) +
  ylim(-3, 3)

# 保存 Z-score 结果到 CSV 文件
write_csv(data_zscore, public_output("各文库微生物多样性11_Zscore.csv"))

# 提示信息
message("Z-score data saved to: ", public_output("各文库微生物多样性11_Zscore.csv"))












# === 加载必要的 R 包 ===
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggalluvial", quietly = TRUE)) install.packages("ggalluvial")
if (!requireNamespace("alluvial", quietly = TRUE)) install.packages("alluvial")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")

library(ggplot2)
library(ggalluvial)
library(alluvial)
library(patchwork)
library(scales)

# === 一键调整全局字体大小 ===
base_fontsize <- 18  # 🔥 这里修改即可统一字体大小

# === 读取数据 ===
wenkutu <- read.csv(public_file("各文库微生物多样性11_Zscore.csv"))
wenkutu$Z_score <- as.numeric(wenkutu$Z_score)
wenkutu$Sample <- factor(wenkutu$Sample, levels = unique(wenkutu$Sample))

# === 统一配色 ===
my_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F")

# === 上方：Alluvial + 堆叠柱状图 ===
p1 <- ggplot(wenkutu, aes(x = Sample, y = Absolute_abundance, alluvium = Type)) +
  geom_alluvium(aes(fill = Type, colour = Type), alpha = 0.4, width = 0.5) + 
  geom_bar(stat = "identity", aes(fill = Type), width = 0.5) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = base_fontsize - 2, family = "Times New Roman", face = "bold"),
    axis.title = element_text(size = base_fontsize, family = "Times New Roman", face = "bold"),
    legend.title = element_text(size = base_fontsize, family = "Times New Roman", face = "bold"),
    legend.text = element_text(size = base_fontsize - 2, family = "Times New Roman"),
    legend.position = "right",
    plot.margin = margin(5,5,5,5)
  ) +
  labs(y = "Absolute Abundance", x = NULL, fill = "Microbial Type", colour = "Microbial Type")

# === 下方：Z-score 折线 + 散点图 ===
p2 <- ggplot(wenkutu, aes(x = Sample)) +
  geom_line(aes(y = Z_score, group = 1), size = 1, color = "#4C4C4C", linetype = "dashed") +
  geom_point(aes(y = Z_score, fill = Type, shape = Type), size = 3.5, color = "black", stroke = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 0.8) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "red", size = 0.8) +
  scale_fill_manual(values = my_colors) +
  scale_shape_manual(values = c(21,22,23,24,25,3)) +
  scale_y_continuous(
    breaks = seq(-1, 5, 1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = base_fontsize - 2, family = "Times New Roman", face = "bold"),
    axis.text.y = element_text(size = base_fontsize - 2, family = "Times New Roman", face = "bold"),
    axis.title = element_text(size = base_fontsize, family = "Times New Roman", face = "bold"),
    legend.title = element_text(size = base_fontsize, family = "Times New Roman", face = "bold"),
    legend.text = element_text(size = base_fontsize - 2, family = "Times New Roman"),
    legend.position = "right",
    plot.margin = margin(5,5,5,5)
  ) +
  labs(y = "Z-score", x = "Sample", shape = "Microbial Type", fill = "Microbial Type")

# === 上下组合图，比例优化 ===
final_plot <- p1 / p2 + plot_layout(heights = c(2, 1))

# === 输出高清 PNG / PDF ===
ggsave(public_output("alluvial_Zscore_NatureStyle.png"),
       plot = final_plot, dpi = 1200, width = 16, height = 12)

ggsave(public_output("alluvial_Zscore_NatureStyle.pdf"),
       plot = final_plot, dpi = 1200, width = 16, height = 12)



