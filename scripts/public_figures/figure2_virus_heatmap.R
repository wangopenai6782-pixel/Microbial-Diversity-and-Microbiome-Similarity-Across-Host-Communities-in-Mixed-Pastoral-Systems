# Author: Wang Linglong
# Public release script derived from: 2\病毒热图.txt
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
library(ComplexHeatmap)
library(circlize)
library(cols4all)
library(randomcoloR)
library(grid)

# === 全局字体大小 ===
base_fontsize <- 15

# === 读取病毒数据 ===
dt <- read.csv(public_file("2病毒科占比.csv"), 
               header = TRUE, row.names = 1)
dt <- as.matrix(dt)

# === 热图颜色函数 ===
col_fun <- colorRamp2(c(0,1,2,3,4,5,6),
                      c("#FEFBFB","#F5D4D4","#EDB4B4","#E39292","#D44545","#BF3636","#A82828"))

# === 创建列注释 ===
virus_labels <- c(rep("RNA virus", 35), rep("DNA virus", ncol(dt)-35))
df_col <- data.frame(Viruses = virus_labels)
ha_col <- HeatmapAnnotation(
  df = df_col,
  col = list(Viruses = c("RNA virus"="#F1B543", "DNA virus"="#7552A3")),
  show_legend = TRUE,
  annotation_name_gp = gpar(fontsize = base_fontsize, fontfamily = "serif")
)

# === 行注释英文名（顺序保持一致） ===
species_labels <- c(
  rep("Humans", 5), 
  rep("Cattle", 14), 
  rep("Sheep", 4), 
  rep("Pigs", 2), 
  rep("Rabbits", 4), 
  rep("Chickens", 2), 
  rep("Ticks", 5), 
  rep("Mosquitoes", 2), 
  rep("Rock pigeon", 4), 
  rep("Eurasian tree sparrow", 5), 
  rep("Ground squirrels", 4), 
  rep("Egrets", 3), 
  rep("Great bustard", 2), 
  rep("Eagle owl", 2), 
  rep("Circus cyaneus", 3),
  rep("House mouse", 2)
)
df_row <- data.frame(Animal_taxonomy = species_labels, stringsAsFactors = FALSE)

# === 固定物种颜色 ===
custom_row_colors <- c(
  "Humans" = "#D92B03",
  "Cattle" = "#F1997B",
  "Sheep" = "#F38B2F",
  "Pigs" = "#A5405E",
  "Rabbits" = "#F4C288",
  "Chickens" = "#088C00",
  "Ticks" = "#8B511F",
  "Mosquitoes" = "#DB6C76",
  "Rock pigeon" = "#7552A7",
  "Eurasian tree sparrow" = "#DCA0DD",
  "Ground squirrels" = "#F1B543",
  "Egrets" = "#BF7533",
  "Great bustard" = "#A38277",
  "Eagle owl" = "#592E13",
  "Circus cyaneus" = "#F2CDCF",
  "House mouse" = "#6E86A5"
)

# === 创建行注释 ===
ha_row_taxonomy <- rowAnnotation(
  df = df_row[, "Animal_taxonomy", drop = FALSE],
  col = list(Animal_taxonomy = custom_row_colors),
  show_legend = TRUE,
  annotation_name_gp = gpar(fontsize = base_fontsize, fontfamily = "serif")
)

# === 绘制热图 ===
main_heatmap <- Heatmap(
  dt,
  name = "RPM(log10)",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  border = "black",
  col = col_fun,
  top_annotation = ha_col,
  left_annotation = ha_row_taxonomy,
  rect_gp = gpar(col = "white", lwd = 1),
  row_names_gp = gpar(fontsize = base_fontsize-7, fontfamily = "serif"),
  row_names_side = "right",
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = base_fontsize, fontfamily = "Times", rot = 90, just = "bottom")
)

# === 输出 PDF ===
pdf(public_output("Virus_Heatmap.pdf"), width = 14, height = 10)
draw(main_heatmap,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)
dev.off()


