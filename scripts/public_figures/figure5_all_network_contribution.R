# Author: Wang Linglong
# Public release script derived from: 5\all贡献网络.txt
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
library(igraph)
library(ggraph)

# 1. 读取数据 ----------------------------------------------------------------
data <- read.csv(public_file("all共线网络.csv"), header = TRUE)

# 2. 定义微生物分类（按列位置）-------------------------------------------------
microbe_categories <- list(
  "Protozoa" = 2:29,
  "Virus" = 30:73,
  "Bacteria" = 74:114,
  "Fungi" = 115:138
)

# 3. 构建物种-微生物存在关系矩阵 ------------------------------------------------
edges <- map_dfr(names(microbe_categories), function(type){
  cols <- microbe_categories[[type]]
  
  data %>%
    pivot_longer(
      cols = all_of(cols),
      names_to = "microbe",
      values_to = "abundance"
    ) %>%
    filter(abundance > 0) %>%
    transmute(
      from = species,
      to = microbe,
      weight = 1,
      microbe_type = type
    ) %>%
    distinct()
})

# 4. 创建节点数据 --------------------------------------------------------------
nodes <- bind_rows(
  # 物种节点
  edges %>%
    distinct(name = from) %>%
    mutate(
      type = "Species",
      size = 8,
      color = "#1F77B4"
    ),
  
  # 微生物节点
  edges %>%
    distinct(name = to, microbe_type) %>%
    mutate(
      type = microbe_type,
      size = 5,
      color = case_when(
        microbe_type == "Protozoa" ~ "#FF7F0E",
        microbe_type == "Virus" ~ "#2CA02C",
        microbe_type == "Bacteria" ~ "#D62728",
        microbe_type == "Fungi" ~ "#9467BD"
      )
    )
)

# 5. 构建网络图 ----------------------------------------------------------------
network <- graph_from_data_frame(
  d = edges,
  directed = FALSE,
  vertices = nodes
)
# 6. 可视化优化 ----------------------------------------------------------------
set.seed(123)
p <- ggraph(network, layout = "fr") +
  # 边设置（加粗且更明显）
  geom_edge_link(
    color = "grey40",   # 加深连线颜色
    width = 0.5,        # 加粗连线
    alpha = 0.6,        # 提高透明度
    edge_linetype = 1   # 实线
  ) +
  
  # 节点设置
  geom_node_point(
    aes(color = color, size = size),
    alpha = 0.8
  ) +
  
  # 微生物标签（自动换行）
  geom_node_text(
    aes(label = str_wrap(name, width = 10)),  # 名称自动换行
    data = function(x) filter(x, type != "Species"),  # 筛选微生物节点
    size = 2.8,
    color = "black",
    repel = TRUE,
    max.overlaps = 50,
    segment.color = NA,
    nudge_y = 0.02
  ) +
  
  # 物种标签（加粗显示）
  geom_node_text(
    aes(label = name),
    data = function(x) filter(x, type == "Species"),  # 筛选物种节点
    size = 3.2,
    color = "#1F77B4",
    repel = FALSE,
    fontface = "bold"
  ) +
  
  # 样式设置
  scale_color_identity() +
  scale_size_identity() +
  labs(title = "物种-微生物共生网络") +
  theme_void() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 20)  # 标题下边距
    ),
    plot.margin = margin(2, 2, 2, 2, "cm")  # 画布边距
  )


# 7. 保存高清大图 --------------------------------------------------------------
output_path <- public_output("优化网络图.png")
ggsave(
  output_path,
  plot = p,
  width = 20,        # 加宽画布
  height = 16,       # 加高画布
  dpi = 600,         # 提高分辨率
  bg = "white"       # 白色背景
)

