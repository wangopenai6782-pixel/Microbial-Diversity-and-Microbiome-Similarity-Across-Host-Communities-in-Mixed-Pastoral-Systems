# Author: Wang Linglong
# Public release script derived from: 5\宿主-宿主+微生物.txt
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
# 加载必要库 ----------------------------------------------------------------
library(tidyverse)
library(igraph)
library(ggraph)
library(ggforce)
library(ggrepel)
library(viridis)
library(ggnewscale)

# 固定随机数种子 ---------------------------------------------------------------
set.seed(2024)

# 1. 读取并处理物种相似性数据（“代码1”部分） ----------------------------------
similarity_data <- read_csv(
  public_file("物种相似性结果.csv"),
  col_names = c("from", "to", "similarity")
) %>%
  mutate(similarity = round(similarity, 2))

threshold        <- quantile(similarity_data$similarity, 0.85)
base_edges       <- filter(similarity_data, similarity > threshold)
all_nodes        <- union(similarity_data$from, similarity_data$to)
connected_nodes  <- union(base_edges$from, base_edges$to)
isolated_nodes   <- setdiff(all_nodes, connected_nodes)
supplement_edges <- similarity_data %>%
  filter((from %in% isolated_nodes) | (to %in% isolated_nodes)) %>%
  group_by(from) %>% slice_max(similarity, n = 1) %>%
  bind_rows(
    similarity_data %>%
      filter((from %in% isolated_nodes) | (to %in% isolated_nodes)) %>%
      group_by(to) %>% slice_max(similarity, n = 1) %>%
      rename(from = to, to = from)
  ) %>%
  distinct(from, to, .keep_all = TRUE) %>%
  anti_join(base_edges, by = c("from", "to"))
final_edges <- bind_rows(base_edges, supplement_edges)

net <- graph_from_data_frame(
  d        = final_edges,
  directed = FALSE,
  vertices = tibble(name = all_nodes)
)
V(net)$degree      <- degree(net)
V(net)$betweenness <- betweenness(net)

clusters       <- cluster_louvain(net, weights = E(net)$similarity)
V(net)$cluster <- as.factor(membership(clusters))
custom_layout  <- create_layout(
  graph   = net,
  layout  = "kk",
  kkconst = vcount(net) * 1.2
)

# 2. 读取并处理宿主–微生物共生数据（“代码2”部分） ----------------------------
data2 <- read.csv(
  public_file("all共线网络.csv"),
  header = TRUE
)
microbe_categories <- list(
  "Protozoa" = 2:29,
  "Virus"    = 30:73,
  "Bacteria" = 74:114,
  "Fungi"    = 115:138
)
edges2 <- map_dfr(names(microbe_categories), function(type) {
  cols <- microbe_categories[[type]]
  data2 %>%
    pivot_longer(
      cols       = all_of(cols),
      names_to   = "microbe",
      values_to  = "abundance"
    ) %>%
    filter(abundance > 0) %>%
    transmute(
      from         = species,
      to           = microbe,
      microbe_type = type
    ) %>%
    distinct()
})

# 3. 计算微生物在同一坐标系下的位置 & 构造连线数据 -------------------------------
species_layout <- as_tibble(custom_layout) %>% select(name, x, y)

set.seed(2025)
microbe_layout <- edges2 %>%
  group_by(to) %>%
  summarise(
    x = mean(species_layout$x[match(from, species_layout$name)]),
    y = mean(species_layout$y[match(from, species_layout$name)]),
    .groups = "drop"
  ) %>%
  mutate(
    x    = x + runif(n(), -0.1, 0.1),
    y    = y + runif(n(), -0.1, 0.1),
    name = to
  )

microbe_segments <- edges2 %>%
  left_join(species_layout, by = c("from" = "name")) %>%
  rename(x = x, y = y) %>%
  left_join(microbe_layout, by = c("to" = "name")) %>%
  rename(xend = x.y, yend = y.y, x = x.x, y = y.x)

# 4. 可视化：主次分层 & 美化 ----------------------------------------------------
cluster_plot_final <- ggraph(custom_layout, by = "name") +

  # A. 最底层：宿主–微生物曲线，用微生物类型对应的柔和色，实线
  ggnewscale::new_scale_color() +
  geom_curve(
    data      = microbe_segments,
    aes(x = x, y = y, xend = xend, yend = yend, color = microbe_type),
    curvature = 0.2,
    linewidth  = 0.3,
    linetype   = "solid"
  ) +
  scale_color_manual(
    name   = "Microbe Link",
    values = c(
      Protozoa = alpha("#FF7F0E", 0.6),
      Virus    = alpha("#2CA02C", 0.6),
      Bacteria = alpha("#E73C36", 0.6),
      Fungi    = alpha("#9467BD", 0.6)
    )
  ) +

  # B. 聚类背景
  geom_mark_hull(
    aes(x, y, group = cluster, fill = cluster),
    concavity = 8, expand = unit(8, "mm"), alpha = 0.15,
    show.legend = FALSE
  ) +
  scale_fill_brewer(palette = "Set2") +
  ggnewscale::new_scale_fill() +

  # C. 宿主–宿主主要网络边
  geom_edge_link(
    aes(width = similarity, color = similarity),
    alpha   = 0.9,
    lineend = "round"
  ) +

  # D. 宿主节点 & 标签
  geom_node_point(
    aes(size = degree, fill = betweenness),
    shape  = 21,
    color  = "black",
    stroke = 0.8,
    alpha  = 0.95
  ) +
  geom_text_repel(
    aes(x = x, y = y, label = name),
    size             = 4.5,
    color            = "black",
    fontface         = "bold",
    family           = "Helvetica",
    min.segment.length = 0,
    max.overlaps     = 100,
    box.padding      = 1.5
  ) +

  # E. 主网络图例设置
  scale_edge_color_gradientn(
    name   = "Similarity",
    colors = c("#fee5d9", "#fc9272", "#de2d26", "#a50f15"),
    limits = c(min(final_edges$similarity), max(final_edges$similarity)),
    guide  = guide_colorbar(order = 1,
                            barheight = unit(8, "cm"),
                            barwidth  = unit(0.5, "cm"))
  ) +
  scale_edge_width(
    name  = "Similarity",
    range = c(0.5, 4),
    guide = guide_legend(order = 1,
      override.aes = list(
        linewidth = c(1, 2, 3, 4),
        color     = c("#fee5d9", "#fc9272", "#de2d26", "#a50f15")
      ))
  ) +
  scale_fill_gradientn(
    name   = "Betweenness",
    colors = c("gray", "red"),
    guide  = guide_colorbar(order = 2)
  ) +
  scale_size_continuous(
    name   = "Degree",
    range  = c(4, 12),
    breaks = scales::breaks_pretty(4),
    guide  = guide_legend(order = 3,
      override.aes = list(shape = 21, fill = "grey50"))
  ) +

  labs(title = "Species Clustering Network with Microbes") +
  theme_void() +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.margin     = unit(c(2,2,2,2), "cm"),
    legend.position = "right",
    legend.box      = "vertical",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 10)
  ) +

  # F. 微生物节点（三角 & 不同颜色）
  ggnewscale::new_scale_fill() +
  geom_point(
    data = microbe_layout %>%
      left_join(distinct(edges2, to, microbe_type),
                by = c("name" = "to")),
    aes(x = x, y = y, fill = microbe_type),
    shape  = 24,
    size   = 3,
    color  = "black",
    stroke = 0.7
  ) +
  scale_fill_manual(
    name   = "Microbe Type",
    values = c(
      Protozoa = "#FF7F0E",
      Virus    = "#2CA02C",
      Bacteria = "#E73C36",
      Fungi    = "#9467BD"
    )
  ) +

  # G. 微生物标签
  geom_text_repel(
    data = microbe_layout,
    aes(x = x, y = y, label = name),
    size          =2,
    color         = "black",
    segment.color = NA,
    nudge_y       = 0.03
  )

# 5. 保存图像 ----------------------------------------------------------------
ggsave(
  filename = file.path(Sys.getenv("USERPROFILE"), "Desktop",
                       "Species_with_Microbes_Final.png"),
  plot     = cluster_plot_final,
  width    = 26,
  height   = 22,
  dpi      = 1000,
  bg       = "white"
)


