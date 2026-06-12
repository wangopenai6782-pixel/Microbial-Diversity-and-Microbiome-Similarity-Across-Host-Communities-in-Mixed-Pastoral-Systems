# Author: Wang Linglong
# Public release script derived from: 4\血餐.txt
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
# === 安装并加载必要的库 ===
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggalluvial", quietly = TRUE)) install.packages("ggalluvial")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
library(ggplot2)
library(ggalluvial)
library(RColorBrewer)

# === 构造数据 ===
data <- data.frame(
  vector = c(rep("Tick", 6), rep("Mosquito", 10)),
  host_species = c(
    "Human (Homo sapiens)", 
    "Large Flying-Fox (Pteropus vampyrus)", 
    "Brandt's Bat (Myotis brandtii)", 
    "Molossus Bat (Molossus molossus)", 
    "Dromedary Camel (Camelus dromedarius)", 
    "House Mouse (Mus musculus)",
    "Alpaca (Vicugna pacos)", 
    "Cattle (Bos taurus)", 
    "Sheep (Ovis aries)", 
    "Pig (Sus scrofa)", 
    "Horse (Equus)", 
    "Dromedary Camel (Camelus dromedarius)", 
    "House Mouse (Mus musculus)", 
    "Human (Homo sapiens)", 
    "Crow (Corvus)", 
    "Donkey (Equus asinus)"
  )
)

# === 为共享宿主分配相同颜色 ===
unique_hosts <- unique(data$host_species)
color_palette <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(length(unique_hosts)),
  unique_hosts
)

# === 全局字体大小设置 ===
global_font_size <- 20

# === 绘制极致美化的桑基图 ===
p <- ggplot(data, aes(axis1 = vector, axis2 = host_species)) +
  # 流动线条
  geom_alluvium(aes(fill = host_species), width = 0.3, alpha = 0.85,
                knot.pos = 0.45, color = "gray50", size = 0.3) +
  # 节点（方框）
  geom_stratum(width = 0.5, fill = "gray95", color = "black", size = 0.6) +
  # 节点标签
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = global_font_size * 0.35, fontface = "bold", family = "sans", 
            color = "black", lineheight = 0.9) +
  # 颜色映射
  scale_fill_manual(values = color_palette) +
  # 主题风格
  theme_minimal(base_size = global_font_size) +
  labs(
    title = "Blood-Feeding Vectors and Their Hosts",
    subtitle = "Shared blood meal hosts highlight the broad ecological contact of ticks and mosquitoes",
    x = "Vector", y = ""
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.text.x = element_text(face = "bold", size = global_font_size, color = "black"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = global_font_size * 1.5),
    plot.subtitle = element_text(hjust = 0.5, size = global_font_size * 1.1, color = "gray30"),
    plot.margin = margin(10, 15, 10, 15)
  )

# === 保存为高分辨率 PDF（横向展示） ===
output_path <- public_output("Sankey_Plot_Nature_Style.pdf")
ggsave(output_path, plot = p, width = 20, height = 8, units = "in", device = "pdf")

print(paste("PDF 图片已保存至:", output_path))


