# Author: Wang Linglong
# Public release script derived from: 1\中国地图内蒙古地图高亮.txt
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
# --- 完整脚本（英文图例 + 统一字体 + 放大字体） ---
library(sf)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(scatterpie)
library(cowplot)
library(ggspatial)
library(Cairo)
library(grid)
library(showtext)

# --- 字体设置（统一为 Arial） ---
showtext_auto(enable = TRUE)
font_add("Arial", regular = "arial.ttf")
base_family <- "Arial"

# --- 读取地图数据 ---
china_map <- st_read(public_file("China_shapefile", "中华人民共和国.shp"), quiet = TRUE)
inner_mongolia <- st_read(public_file("nmg_shapefile", "内蒙古自治区.shp"), quiet = TRUE)

# --- 高亮城市 ---
highlight_cities <- c("呼和浩特市", "鄂尔多斯市", "包头市", "锡林郭勒盟")
highlighted <- inner_mongolia %>% filter(name %in% highlight_cities)

# --- 读取并处理饼图数据 ---
pie_raw <- read_csv(public_file("species_density_by_species.csv"))

pie_data <- pie_raw %>%
  select(location_id, species, total_species_abundance, lon, lat) %>%
  pivot_wider(names_from = species, values_from = total_species_abundance, values_fill = 0) %>%
  distinct()

pie_species_cols <- setdiff(colnames(pie_data), c("location_id", "lon", "lat"))
pie_data$total <- rowSums(pie_data[, pie_species_cols])

# --- 配色 ---
pie_colors <- c(
  "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
  "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
  "#ffff99", "#b15928", "#8dd3c7", "#ffffb3", "#bebada",
  "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
  "#d9d9d9"
)[1:length(pie_species_cols)]

# --- 主图 ---
p_main <- ggplot() +
  geom_sf(data = china_map, fill = "#FFFFFF", color = "#000000", size = 0.8) +
  geom_sf(data = inner_mongolia, fill = "#FCF8E8", color = "#BFB69B", size = 0.5) +
  geom_sf(data = highlighted, fill = "#E6D9B8", color = "#BFB69B", size = 0.7) +
  geom_point(data = pie_data, aes(x = lon, y = lat),
             shape = 16, color = "#D62728", size = 1.2, alpha = 0.9) +
  theme_minimal(base_size = 16, base_family = base_family) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "#FFFFFF", color = NA),
    plot.margin = margin(15, 15, 15, 15)
  )

# --- 饼图大小图例函数 ---
make_nested_circle_legend <- function(pie_data, scale_range = c(0.25, 0.7), base_family="Arial") {
  size_values <- c(min(pie_data$total), median(pie_data$total), max(pie_data$total))
  size_values_int <- round(size_values)
  
  radii_scaled <- scales::rescale(size_values, to = scale_range)
  
  cx <- 0; cy <- 0
  circle_angles <- seq(0, 2*pi, length.out = 200)
  
  circle_paths <- do.call(rbind, lapply(seq_along(radii_scaled), function(i) {
    r <- radii_scaled[i]
    data.frame(
      x = cx + r * cos(circle_angles),
      y = cy + r * sin(circle_angles),
      idx = factor(i),
      radius = r
    )
  }))
  
  total_labels <- length(radii_scaled)
  label_height <- 0.4
  start_y <- cy + max(radii_scaled) + label_height
  y_spacing <- label_height * 2.8
  label_ys <- rev(seq(start_y, start_y - y_spacing*(total_labels - 1), length.out = total_labels))
  
  label_x_const <- cx + max(radii_scaled) + 1.0
  
  df_labels <- data.frame(
    idx = seq_along(radii_scaled),
    radius = radii_scaled,
    label = size_values_int,
    label_x = label_x_const,
    label_y = label_ys,
    stringsAsFactors = FALSE
  )
  
  df_labels <- df_labels %>%
    rowwise() %>%
    mutate(
      dy_target = label_y - cy,
      dy_clamped = ifelse(abs(dy_target) > radius*0.95,
                          sign(dy_target) * radius*0.95,
                          dy_target),
      y_on_circle = cy + dy_clamped,
      x_on_circle = cx + sqrt(pmax(radius^2 - (y_on_circle - cy)^2, 0))
    ) %>%
    ungroup()
  
  ggplot() +
    geom_path(data = circle_paths, aes(x = x, y = y, group = idx),
              color = "black", size = 0.9, alpha = 0.8) +
    geom_segment(data = df_labels,
                 aes(x = x_on_circle, y = y_on_circle,
                     xend = label_x - 0.18, yend = label_y),
                 arrow = arrow(length = unit(0.15, "cm")),
                 color = "black", size = 0.7) +
    geom_text(data = df_labels,
              aes(x = label_x, y = label_y, label = label),
              hjust = 0, vjust = 0.5, size = 5, family = base_family) +
    coord_fixed() +
    xlim(-max(radii_scaled)*1.2, label_x_const + 1.5) +
    ylim(cy - max(radii_scaled)*1.5, start_y + label_height) +
    theme_void() +
    labs(title = "Abundance") +
    theme(
      plot.title = element_text(face = "bold", size = 14, family = base_family)
    )
}

# --- 生成圆圈图例 ---
circle_legend <- make_nested_circle_legend(pie_data, base_family = base_family)

# --- 城市放大图 ---
bbox_cities <- st_bbox(highlighted)
p_zoom <- ggplot() +
  geom_sf(data = highlighted, fill = "#E6D9B8", color = "#BFB69B", size = 1.2) +
  geom_scatterpie(
    data = pie_data,
    aes(x = lon, y = lat, r = scales::rescale(total, to = c(0.25, 0.7))),
    cols = pie_species_cols,
    color = NA, alpha = 0.95
  ) +
  scale_fill_manual(values = pie_colors, name = "Species") +
  coord_sf(xlim = c(bbox_cities$xmin, bbox_cities$xmax),
           ylim = c(bbox_cities$ymin, bbox_cities$ymax)) +
  annotation_scale(location = "bl", width_hint = 0.5,
                   line_width = 0.5, text_cex = 1.2) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, family = base_family),
    legend.text = element_text(size = 12, family = base_family)
  )

legend_species <- cowplot::get_legend(p_zoom + theme(legend.position = "right"))
p_zoom_clean <- p_zoom + theme(legend.position = "none")

# --- 主体拼图 ---
main_grid <- cowplot::plot_grid(
  p_main, p_zoom_clean,
  ncol = 2, rel_widths = c(2.5, 3)
)

# --- 右侧竖直堆叠：species 在上，abundance 在下 ---
right_legends <- cowplot::plot_grid(
  legend_species, circle_legend,
  ncol = 1, rel_heights = c(1, 1.2)
)

# --- 最终合并并保存 PDF ---
final_plot <- cowplot::plot_grid(
  main_grid, right_legends,
  ncol = 2, rel_widths = c(8, 2)
)

ggsave(public_output("final_map_species_abundance_english_font.pdf"),
       final_plot, width = 20, height = 8, device = cairo_pdf)



