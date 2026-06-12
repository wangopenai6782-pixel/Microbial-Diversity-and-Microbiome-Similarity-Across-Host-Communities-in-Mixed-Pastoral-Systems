# Author: Wang Linglong
# Public release script derived from: 3\基因组图.txt
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
# === 加载包 ===
library(tidyverse)
library(ggsci)
library(showtext)
library(grid)

# === 字体设置 ===
showtext_auto()
font_add("Arial", "arial.ttf")
base_family <- "Arial"

# === 数据读取 ===
df <- read_csv(public_file("新病毒.csv"),
               col_types = cols(
                 Virus.species = col_character(),
                 Virus.family = col_character(),
                 Sequence.length = col_double(),
                 N = col_character(),
                 P = col_character(),
                 M = col_character(),
                 G = col_character(),
                 L = col_character(),
                 nsP1a = col_character(),
                 RDRP = col_character(),
                 CP = col_character(),
                 Polyprotein = col_character()
               ))

# === ORF 数据整理 ===
orf_cols <- c("N", "P", "M", "G", "L", "nsP1a", "RDRP", "CP", "Polyprotein")

df_long <- df %>%
  pivot_longer(cols = all_of(orf_cols),
               names_to = "ORF_name", values_to = "range") %>%
  filter(!is.na(range)) %>%
  mutate(range = str_replace_all(range, '"', ''),
         range = str_trim(range),
         start = as.numeric(str_extract(range, "^[0-9]+")),
         end = as.numeric(str_extract(range, "[0-9]+$")),
         direction = ifelse(end < start, "-", "+"))

# === 固定 ORF 顺序 & 色板 ===
orf_levels <- orf_cols[orf_cols %in% unique(df_long$ORF_name)]
df_long$ORF_name <- factor(df_long$ORF_name, levels = orf_levels)
orf_colors <- pal_npg("nrc")(length(orf_levels))
names(orf_colors) <- orf_levels

# === Y 轴顺序（按 CSV 文件顺序） ===
virus_order <- df$Virus.species
df_long$Virus.species <- factor(df_long$Virus.species, levels = virus_order)
df$Virus.species <- factor(df$Virus.species, levels = virus_order)

# === 绘图 ===
p <- ggplot() +
  # 基因组灰线
  geom_segment(data = df, 
               aes(x = 0, xend = Sequence.length,
                   y = as.numeric(Virus.species),
                   yend = as.numeric(Virus.species)),
               color = "gray70", size = 1.5) +
  # ORF 矩形（高度增加）
  geom_rect(data = df_long,
            aes(xmin = pmin(start, end), xmax = pmax(start, end),
                ymin = as.numeric(Virus.species)-0.4,
                ymax = as.numeric(Virus.species)+0.4,
                fill = ORF_name),
            color = "black", size = 0.5) +
  # ORF 标签（矩形内，字体增大）
  geom_text(data = df_long,
            aes(x = (start+end)/2, y = as.numeric(Virus.species), 
                label = ORF_name),
            size = 14, family = base_family, color = "black", fontface = "bold") +
  # 基因组长度标注（右端，字体增大）
  geom_text(data = df,
            aes(x = Sequence.length + 0.02*max(df$Sequence.length), 
                y = as.numeric(Virus.species),
                label = Sequence.length),
            hjust = 0, size = 9, family = base_family, fontface = "bold", color = "black") +
  # Y 轴病毒名称
  scale_y_continuous(breaks = 1:length(virus_order),
                     labels = virus_order, expand = expansion(add = c(0.6,0.6))) +
  # ORF 填充颜色（去掉图例）
  scale_fill_manual(values = orf_colors, guide = "none") +
  # 美化主题
  theme_bw(base_size = 24) +
  labs(x = "Genome position", y = "") +
  theme(
        axis.text.y = element_text(size = 38, face = "bold", family = base_family, color = "black"),
        axis.text.x = element_text(size = 32, family = base_family, color = "black"),
        axis.title = element_text(size = 24, face = "bold", family = base_family),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", size = 0.6),
        plot.margin = margin(20,20,20,20))

# === 输出 PNG & PDF（宽图高分辨率 + 矢量） ===
ggsave(public_output("virus_genome_pubstyle_wide.png"), plot = p, width = 32, height = 16, dpi = 600)
ggsave(public_output("virus_genome_pubstyle_wide.pdf"), plot = p, width = 34, height = 16)

# === 展示图形 ===
print(p)


