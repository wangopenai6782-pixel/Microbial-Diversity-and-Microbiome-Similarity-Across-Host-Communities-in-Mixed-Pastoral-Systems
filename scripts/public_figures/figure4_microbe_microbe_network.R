# Author: Wang Linglong
# Public release script derived from: 4\微生物-微生物.txt
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
### Microbial–Microbial Co-occurrence Network (Nature-style, compact layout, font & legend improved)
library(tidyverse)
library(psych)
library(igraph)
library(ggraph)
library(tidygraph)
library(cowplot)
library(showtext)

# === 字体设置 ===
showtext_auto(TRUE)
font_add("Arial", regular = "arial.ttf")

# === 用户可调全局字体大小 ===
GLOBAL_FONT_SIZE <- 24

# === 1. 读取数据 ===
data <- read.csv(public_file("all共线网络.csv"), header = TRUE)
colnames(data) <- make.names(colnames(data), unique = TRUE)

# === 2. 微生物分类 ===
microbe_categories <- list(
  "Protozoa" = 2:29,
  "Virus" = 30:73,
  "Bacteria" = 74:114,
  "Fungi" = 115:138
)
microbe_names <- unlist(map(microbe_categories, ~ colnames(data)[.x]))
microbe_data <- data[, microbe_names]

# === 3. 相关性矩阵（Spearman） ===
cor_result <- corr.test(microbe_data, method = "spearman", adjust = "none")
cor_matrix <- cor_result$r
p_matrix <- cor_result$p

# === 4. 筛选正相关边（R > 0.6 & p < 0.05） ===
threshold <- 0.6
sig_edges <- which((cor_matrix > threshold) & (p_matrix < 0.05), arr.ind = TRUE) %>%
  as.data.frame() %>%
  filter(row < col) %>%
  mutate(
    from = rownames(cor_matrix)[row],
    to   = colnames(cor_matrix)[col],
    weight = cor_matrix[cbind(row, col)]
  ) %>%
  select(from, to, weight)

# === 5. 构建 igraph 对象 ===
g <- graph_from_data_frame(sig_edges, directed = FALSE)

# === 6. 节点属性 ===
node_info <- tibble(name = V(g)$name) %>%
  mutate(
    microbe_type = case_when(
      name %in% colnames(data)[microbe_categories$Protozoa] ~ "Protozoa",
      name %in% colnames(data)[microbe_categories$Virus] ~ "Virus",
      name %in% colnames(data)[microbe_categories$Bacteria] ~ "Bacteria",
      name %in% colnames(data)[microbe_categories$Fungi] ~ "Fungi",
      TRUE ~ "Unknown"
    ),
    degree = degree(g)[name],
    betweenness = betweenness(g)[name]
  )

V(g)$microbe_type <- node_info$microbe_type
V(g)$degree <- node_info$degree
V(g)$betweenness <- node_info$betweenness

# === 7. 颜色定义 ===
microbe_colors <- c(
  "Protozoa" = "#FF7F0E",
  "Virus"    = "#2CA02C",
  "Bacteria" = "#D62728",
  "Fungi"    = "#9467BD",
  "Unknown"  = "#999999"
)

# === 8. 紧凑布局 ===
set.seed(123)
layout_fr <- layout_with_fr(
  g,
  weights = E(g)$weight,
  niter = 1000,
  area = vcount(g)^2 * 0.7,
  repulserad = vcount(g)^1.5
)

# === 9. 绘图 ===
p <- ggraph(g, layout = layout_fr) +
  # 边
  geom_edge_link(aes(width = weight, alpha = weight), color = "grey60") +
  scale_edge_width(range = c(0.2, 1)) +
  scale_edge_alpha(range = c(0.2, 0.8)) +

  # 节点
  geom_node_point(aes(fill = microbe_type, size = degree),
                  shape = 21, color = "white", stroke = 0.3, alpha = 0.9) +

  # 标签
  geom_node_text(aes(label = name),
                 repel = TRUE,
                 size = GLOBAL_FONT_SIZE * 0.15,
                 color = "black",
                 segment.color = "grey80",
                 max.overlaps = Inf) +

  # 颜色和大小映射
  scale_fill_manual(
    values = microbe_colors,
    name = "Microbe Type",
    guide = guide_legend(
      override.aes = list(size = 6) # 图例圆点放大
    )
  ) +
  scale_size_continuous(
    range = c(2, 8),
    name = "Node Degree",
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        fill = "black",
        color = "white",
        stroke = 0.3,
        alpha = 0.9
      )
    )
  ) +

  # 极简主题
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(hjust = 0.5, size = GLOBAL_FONT_SIZE + 2, face = "bold"),
    legend.title = element_text(size = GLOBAL_FONT_SIZE),
    legend.text = element_text(size = GLOBAL_FONT_SIZE * 0.8)
  ) +
  labs(title = "Microbe–Microbe Co-occurrence Network")

# === 10. 保存 PDF ===
output_pdf <- public_file("Microbe_CoOccurrence_Network_LegendLargeDot.pdf")
ggsave(output_pdf, plot = p, width = 12, height = 10, units = "in", dpi = 300)

message("✅ 网络图已保存，Microbe Type 图例圆点加大：Desktop/Microbe_CoOccurrence_Network_LegendLargeDot.pdf")


