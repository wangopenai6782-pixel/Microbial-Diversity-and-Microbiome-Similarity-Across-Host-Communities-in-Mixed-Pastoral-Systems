# Author: Wang Linglong
# Public release script derived from: 5\血餐.txt
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
# 安装并加载必要的库
install.packages("ggplot2")
install.packages("ggalluvial")
install.packages("RColorBrewer")
library(ggplot2)
library(ggalluvial)
library(RColorBrewer)

# 构造数据
data <- data.frame(
  vector = c(rep("Dermacentor silvarum", 6), rep("Aedes sp.", 10)),  # 蜱虫 or 蚊子
  host_species = c("homo", "Pteropus vampyrus", "Myotis brandtii", "Molossus molossus", 
                   "Camelus dromedarius", "Mus musculus",
                   "Vicugna pacos", "Bos taurus", "Ovis aries", "Sus scrofa", 
                   "Equus", "Camelus dromedarius", "Mus musculus", "homo", 
                   "Corvu", "Equus asinus")
)

# 生成足够的颜色
num_colors <- length(unique(data$host_species))  # 计算需要的颜色数量
color_palette <- colorRampPalette(brewer.pal(8, "Set2"))(num_colors)  # 生成足够颜色

# 绘制桑基图
p <- ggplot(data, aes(axis1 = vector, axis2 = host_species)) +
  geom_alluvium(aes(fill = host_species), width = 0.25, alpha = 0.85, knot.pos = 0.4) +  # 平滑流动线
  geom_stratum(width = 0.3, fill = "gray90", color = "black") +  # 让分类更清晰
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4, fontface = "bold") +  # 增大字体
  scale_fill_manual(values = color_palette) +  # 使用高质量配色
  theme_minimal(base_size = 14) +  # 增加基础字体大小
  labs(title = "Blood Meal Hosts of Ticks and Mosquitoes",
       subtitle = "Distribution of Hosts for Different Blood-Feeding Vectors",
       x = "Vector", y = "") +
  theme(
    legend.position = "none", 
    axis.text.y = element_blank(),  # 移除 y 轴刻度
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  # 标题居中
    plot.subtitle = element_text(hjust = 0.5, size = 14),  # 副标题居中
    axis.text.x = element_text(face = "bold")  # 突出 X 轴信息
  )

# 保存高分辨率图片 (1200 DPI)
output_path <- public_output("Sankey_Plot.png")  # macOS/Linux桌面路径
# output_path <- public_output("Sankey_Plot.png")
ggsave(output_path, plot = p, dpi = 1200, width = 10, height = 6, units = "in", device = "png")

# 提示输出路径
print(paste("图片已保存至:", output_path))



