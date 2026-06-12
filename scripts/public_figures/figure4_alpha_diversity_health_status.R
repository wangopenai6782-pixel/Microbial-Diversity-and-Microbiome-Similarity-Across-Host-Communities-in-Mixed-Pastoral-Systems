# Author: Wang Linglong
# Public release script derived from: 4\FIG5D家畜疾病VS健.txt
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
# === 加载所需包 ===
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)

# === 读取数据 ===
data <- read.csv(public_file("健康VS疾病.csv"), header = TRUE)
colnames(data) <- make.names(colnames(data), unique = TRUE)

# === 定义微生物列 ===
protozoa_names <- colnames(data)[2:29]
virus_names    <- colnames(data)[30:73]
bacteria_names <- colnames(data)[74:114]
fungi_names    <- colnames(data)[115:138]

# === 添加分组信息 ===
data$Group <- factor(data$Group, levels = c("Healthy", "Diseased"))

# === 计算 alpha 多样性 ===
alpha_df <- data.frame(
  Sample = 1:nrow(data),
  Group = data$Group,
  Protozoa  = diversity(data[, protozoa_names], index = "shannon"),
  Virus     = diversity(data[, virus_names], index = "shannon"),
  Bacteria  = diversity(data[, bacteria_names], index = "shannon"),
  Fungi     = diversity(data[, fungi_names], index = "shannon")
)

# === 长格式转换 ===
alpha_long <- alpha_df %>%
  pivot_longer(cols = c(Protozoa, Virus, Bacteria, Fungi),
               names_to = "Microbe_Type", values_to = "Shannon")

# === 计算显著性 p 值 ===
pvals <- alpha_long %>%
  group_by(Microbe_Type) %>%
  summarise(p = t.test(Shannon ~ Group)$p.value) %>%
  mutate(
    signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    )
  )

# === 添加 y 轴用于显著性标注 ===
pvals <- pvals %>%
  rowwise() %>%
  mutate(
    y_pos = max(alpha_long$Shannon[alpha_long$Microbe_Type == Microbe_Type]) + 0.1,
    x_start = 1,
    x_end   = 2
  )

# === 全局字体大小设置 ===
global_font_size <- 16

# === 绘图 ===
p <- ggplot(alpha_long, aes(x = Group, y = Shannon, fill = Group)) +
  # 小提琴图
  geom_violin(trim = FALSE, scale = "width", width = 0.8,
              color = "black", alpha = 0.25) +
  
  # 内嵌箱线图
  geom_boxplot(width = 0.15, position = position_dodge(0.75),
               outlier.shape = NA, color = "black", alpha = 0.6) +
  
  # 散点图
  geom_jitter(shape = 21, size = global_font_size * 0.2, alpha = 0.7,
              width = 0.1, stroke = 0.3, color = "black") +
  
  # Facet 分面图
  facet_wrap(~ Microbe_Type, nrow = 1, scales = "free_y") +
  
  # 连线表示组间差异
  geom_segment(data = pvals, aes(x = x_start, xend = x_end, y = y_pos, yend = y_pos),
               inherit.aes = FALSE, color = "black", size = 0.5) +
  
  # 显著性标注在连线上方
  geom_text(data = pvals, aes(x = 1.5, y = y_pos + 0.03, label = signif),
            inherit.aes = FALSE, size = global_font_size * 0.35, fontface = "bold") +
  
  # 自定义色板（Nature 风格）
  scale_fill_manual(values = c("Healthy" = "#1f77b4", "Diseased" = "#ff7f0e")) +
  
  # 主题美化
  theme_minimal(base_family = "Times New Roman", base_size = global_font_size) +
  theme(
    strip.text = element_text(face = "bold", size = global_font_size * 0.9),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = global_font_size),
    axis.text.x = element_text(face = "bold", size = global_font_size * 0.8),
    axis.text.y = element_text(size = global_font_size * 0.8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", size = 0.3),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = global_font_size * 0.8),
    plot.title = element_text(face = "bold", size = global_font_size * 1.2, hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    y = "Shannon Diversity Index",
    title = "Alpha Diversity Across Microbial Groups"
  ) +
  coord_cartesian(ylim = c(0, max(alpha_long$Shannon) + 0.5))

# === 显示图 ===
print(p)

# === 导出高清PDF/PNG（投稿用） ===
ggsave(public_output("AlphaDiversity_NatureStyle_Line.pdf"), plot = p, width = 12, height = 6, device = "pdf")
ggsave(public_output("AlphaDiversity_NatureStyle_Line.png"), plot = p, width = 12, height = 6, dpi = 800)


