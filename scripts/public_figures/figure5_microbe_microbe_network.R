# Author: Wang Linglong
# Public release script derived from: 5\微生物-微生物.txt
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
### Microbial–Microbial Co-occurrence Network (No Species Nodes)
library(tidyverse)
library(psych)
library(igraph)
library(ggraph)
library(tidygraph)
library(cowplot)

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

# === 3. 相关性矩阵（Spearman）===
cor_result <- corr.test(microbe_data, method = "spearman", adjust = "none")
cor_matrix <- cor_result$r
p_matrix <- cor_result$p

# === 4. 筛选正相关边（R > 0.6 & p < 0.05）===
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

# === 5. 构建网络对象 ===
g <- graph_from_data_frame(sig_edges, directed = FALSE)

# === 6. 微生物类型注释 ===
node_info <- tibble(name = V(g)$name) %>%
  mutate(
    microbe_type = case_when(
      name %in% colnames(data)[microbe_categories$Protozoa] ~ "Protozoa",
      name %in% colnames(data)[microbe_categories$Virus] ~ "Virus",
      name %in% colnames(data)[microbe_categories$Bacteria] ~ "Bacteria",
      name %in% colnames(data)[microbe_categories$Fungi] ~ "Fungi",
      TRUE ~ "Unknown"
    ),
    betweenness = betweenness(g)[name],
    degree = degree(g)[name]
  )

V(g)$microbe_type <- node_info$microbe_type
V(g)$betweenness <- node_info$betweenness
V(g)$degree <- node_info$degree

# === 7. 颜色定义 ===
microbe_colors <- c(
  "Protozoa" = "#FF7F0E",
  "Virus"    = "#2CA02C",
  "Bacteria" = "#D62728",
  "Fungi"    = "#9467BD",
  "Unknown"  = "#999999"
)

# === 8. 绘图 ===
set.seed(123)
p <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(edge_width = weight), color = "grey60", alpha = 0.3) +
  geom_node_point(aes(fill = microbe_type, size = betweenness), shape = 21, color = "white", stroke = 0.2) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 2.5,
    color = "black",
    segment.color = "grey80",
    max.overlaps = 100
  ) +
  scale_fill_manual(values = microbe_colors) +
  scale_edge_width(range = c(0.1, 1)) +
  scale_size_continuous(range = c(1, 6)) +
  theme_void(base_family = "Arial") +
  labs(title = "微生物共线网络（Microbe–Microbe Co-occurrence Network）") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  ) +
  guides(
    fill = guide_legend(title = "Microbe Type"),
    size = guide_legend(title = "Betweenness")
  )

# === 9. 保存图像 ===
output_path <- public_output("微生物共线网络.tiff")
ggsave(
  output_path,
  plot = p,
  width = 18,
  height = 16,
  units = "cm",
  dpi = 600,
  bg = "white",
  compression = "lzw"
)

message("✅ 微生物共线网络图已保存：Desktop/微生物共线网络.tiff")


