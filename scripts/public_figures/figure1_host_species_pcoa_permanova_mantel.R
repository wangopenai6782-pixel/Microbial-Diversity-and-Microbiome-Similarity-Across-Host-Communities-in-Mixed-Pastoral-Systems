# Author: Wang Linglong
# Public release script derived from: 1\新建 文本文档.txt
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
# ============================================================
# PCoA of host species composition based on Bray Curtis distance
# PERMANOVA: habitat and city effects
# Input: species_data_with_richness原始未标准化.csv
# ============================================================

# === 加载必要包 ===
library(tidyverse)
library(vegan)
library(showtext)
library(grid)

# === 字体设置 ===
showtext_auto(enable = TRUE)

# 如果当前系统有 arial.ttf，则使用 Arial；否则使用 sans
try({
  font_add("Arial", regular = "arial.ttf")
}, silent = TRUE)

base_family <- if ("Arial" %in% sysfonts::font_families()) "Arial" else "sans"

# === 一键字体函数 ===
set_global_fontsize <- function(base_size = 18, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "#e0e0e0", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.8),
      axis.title  = element_text(face = "bold", size = rel(1.2), color = "#222222"),
      axis.text   = element_text(size = rel(1.0), color = "#444444"),
      axis.ticks  = element_line(color = "#444444", linewidth = 0.7),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = rel(1.1), color = "#222222"),
      legend.text  = element_text(size = rel(1.0), color = "#222222"),
      legend.background = element_rect(fill = NA, color = NA),
      legend.key = element_rect(fill = NA, color = NA),
      legend.spacing.y = unit(0.3, "cm"),
      legend.box.spacing = unit(0.5, "cm"),
      plot.margin = margin(15, 20, 15, 15)
    )
}

# === 1. 读取数据 ===
data_path <- public_file("species_data_with_richness原始未标准化.csv")
outdir <- public_output("PCoA_host_species_composition_output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

data <- read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)

# === 2. 检查必要字段 ===
required_cols <- c("city", "location_id", "species", "habitat", "species_density")
missing_cols <- setdiff(required_cols, colnames(data))

if (length(missing_cols) > 0) {
  stop(paste("缺少必要字段:", paste(missing_cols, collapse = ", ")))
}

# === 3. 统一 habitat 名称，并创建 SampleID ===
# 注意：SampleID 必须代表一个采样点，而不是每一行
data <- data %>%
  mutate(
    city = as.character(city),
    location_id = as.character(location_id),
    species = as.character(species),
    habitat = as.character(habitat),
    habitat = recode(
      habitat,
      "Farm" = "Pasture",
      "Wild" = "Natural grassland",
      "pasture" = "Pasture",
      "grassland" = "Natural grassland",
      "Grassland" = "Natural grassland"
    ),
    SampleID = paste(city, location_id, habitat, sep = "_")
  )

# === 4. 构建群落矩阵：采样点 × 物种 ===
community_matrix <- data %>%
  group_by(SampleID, species) %>%
  summarise(
    density = sum(species_density, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = species,
    values_from = density,
    values_fill = 0
  ) %>%
  column_to_rownames("SampleID")

# === 5. 元数据 ===
metadata <- data %>%
  select(SampleID, city, habitat) %>%
  distinct() %>%
  column_to_rownames("SampleID")

# 保证 metadata 顺序与 community_matrix 一致
metadata <- metadata[rownames(community_matrix), , drop = FALSE]

# === 6. 去除全零样本 ===
all_zero_rows <- apply(community_matrix, 1, function(x) all(x == 0))

community_matrix_filtered <- community_matrix[!all_zero_rows, , drop = FALSE]
metadata_filtered <- metadata[rownames(community_matrix_filtered), , drop = FALSE]

cat("Removed", sum(all_zero_rows), "samples with zero species.\n")
cat("Remaining samples:", nrow(community_matrix_filtered), "\n")
cat("Number of species:", ncol(community_matrix_filtered), "\n")

# === 7. 设置因子顺序 ===
metadata_filtered <- metadata_filtered %>%
  mutate(
    habitat = factor(habitat, levels = c("Pasture", "Natural grassland")),
    city = factor(city)
  )

# === 8. Bray Curtis 距离矩阵 ===
dist_mat <- vegdist(community_matrix_filtered, method = "bray")

# === 9. PCoA分析 ===
pcoa <- cmdscale(dist_mat, eig = TRUE, k = 2)

pcoa_points <- as.data.frame(pcoa$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points$SampleID <- rownames(pcoa_points)

eig_vals <- pcoa$eig
var_exp <- round(100 * eig_vals[1:2] / sum(eig_vals[eig_vals > 0]), 1)

pcoa_data <- pcoa_points %>%
  left_join(
    metadata_filtered %>% rownames_to_column("SampleID"),
    by = "SampleID"
  )

# === 10. PERMANOVA分析 ===
adonis_result <- adonis2(
  dist_mat ~ habitat + city,
  data = metadata_filtered,
  permutations = 999,
  by = "terms"
)

print(adonis_result)

# 输出 PERMANOVA 结果
permanova_df <- as.data.frame(adonis_result) %>%
  rownames_to_column("Term")

write.csv(
  permanova_df,
  file.path(outdir, "PERMANOVA_results.csv"),
  row.names = FALSE
)

# === 11. 输出 PCoA 坐标 ===
write.csv(
  pcoa_data,
  file.path(outdir, "PCoA_coordinates.csv"),
  row.names = FALSE
)

# === 12. 配色和形状 ===
# 按你原来代码的配色
habitat_colors <- c(
  "Pasture" = "#D62728",
  "Natural grassland" = "#1F77B4"
)

# 城市形状
city_levels <- levels(metadata_filtered$city)
shape_values <- seq(16, 16 + length(city_levels) - 1)
names(shape_values) <- city_levels

# === 13. PCoA绘图 ===
pcoa_plot <- ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2)) +
  geom_point(
    aes(color = habitat, shape = city),
    size = 4.5,
    stroke = 1.2,
    alpha = 0.75
  ) +
  scale_color_manual(values = habitat_colors, drop = FALSE) +
  scale_shape_manual(values = shape_values, drop = FALSE) +
  set_global_fontsize(base_size = 22, base_family = base_family) +
  labs(
    x = paste0("PCoA Axis1 (", var_exp[1], "%)"),
    y = paste0("PCoA Axis2 (", var_exp[2], "%)"),
    color = "Habitat type",
    shape = "City"
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 5, alpha = 1)),
    shape = guide_legend(override.aes = list(size = 5, alpha = 1))
  )

print(pcoa_plot)

# === 14. 保存图片 ===
ggsave(
  filename = file.path(outdir, "PCoA_host_species_composition.pdf"),
  plot = pcoa_plot,
  width = 9,
  height = 7,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(outdir, "PCoA_host_species_composition.png"),
  plot = pcoa_plot,
  width = 9,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)

# === 15. 输出结果摘要 ===
cat("\n=== PCoA variance explained ===\n")
cat("PCoA Axis1:", var_exp[1], "%\n")
cat("PCoA Axis2:", var_exp[2], "%\n")

cat("\n=== PERMANOVA result ===\n")
print(adonis_result)

cat("\nFiles saved to:\n")
cat(outdir, "\n")

