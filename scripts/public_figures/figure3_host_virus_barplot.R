# Author: Wang Linglong
# Public release script derived from: 3\条形图.txt
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

# === 字体设置 ===
showtext_auto()
font_add("Arial", "arial.ttf")
base_family <- if("Arial" %in% sysfonts::font_families()) "Arial" else "sans"

# 一键调整全局字体大小
global_font_size <- 26  

# === 数据读取 ===
df <- read_csv(public_file("跨物种.csv")) %>%
  rename(host = `Host.species`,
         virus_family = `Virus.taxo(order:family)`,
         RPM = RPM) %>%
  filter(!is.na(host), !is.na(virus_family))

# === 每个物种的病毒数量（次数） ===
host_virus_count <- df %>%
  group_by(host) %>%
  summarise(virus_count = n(), .groups = "drop") %>%
  arrange(desc(virus_count))

# === 宿主顺序（数量多的在上面） ===
host_virus_count$host <- factor(
  host_virus_count$host,
  levels = rev(host_virus_count$host)
)

# === 统一颜色设置（Nature 风格蓝灰色） ===
nature_color <- "#4C72B0"

# === 绘图（横向条形图 + Nature 风格） ===
p <- ggplot(host_virus_count, aes(y = host, x = virus_count)) +
  geom_col(fill = nature_color, color = "black", width = 0.7) +
  geom_text(aes(label = virus_count), hjust = -0.2, 
            size = global_font_size * 0.3, family = base_family) +
  labs(x = "Number of viruses", y = "Host species") +
  theme_classic(base_family = base_family, base_size = global_font_size) +
  theme(
    axis.text.x = element_text(size = global_font_size * 0.9, color="black"),
    axis.text.y = element_text(size = global_font_size * 0.9, color="black", face="italic"), # 宿主名斜体
    axis.title = element_text(face="bold", size = global_font_size * 1.1),
    plot.margin = margin(20, 40, 20, 20),   # 右边留空间放标签
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))   # 给右边留点空间放数值

# === 输出高分辨率图 ===
ggsave(public_output("Host_Virus_BarPlot_Horizontal_NatureMono.pdf"), plot=p, width=14, height=10, units="in", dpi=600)
ggsave(public_output("Host_Virus_BarPlot_Horizontal_NatureMono.png"), plot=p, width=14, height=10, units="in", dpi=600)


