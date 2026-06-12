# Author: Wang Linglong
# Public release script derived from: 1\不同地区物种组成差异.txt
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
# === 加载必要包 ===
library(tidyverse)
library(vegan)
library(showtext)
library(ggsci)
library(grid)
library(ggforce)  # 用于 geom_sina 备选

# === 字体设置 ===
showtext_auto(enable = TRUE)
font_add("Arial", regular = "arial.ttf")
base_family <- if("Arial" %in% sysfonts::font_families()) "Arial" else "sans"

# === 一键字体函数 ===
set_global_fontsize <- function(base_size = 18, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "#e0e0e0", size = 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey70", fill = NA, size = 0.8),
      axis.title  = element_text(face = "bold", size = rel(1.2), color = "#222222"),
      axis.text   = element_text(size = rel(1.0), color = "#444444"),
      axis.ticks  = element_line(color = "#444444", size = 0.7),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = rel(1.1), color = "#222222"),
      legend.text  = element_text(size = rel(1.0), color = "#222222"),
      legend.background = element_rect(fill = NA),
      legend.key = element_rect(fill = NA),
      legend.spacing.y = unit(0.3, "cm"),
      legend.box.spacing = unit(0.5, "cm"),
      plot.margin = margin(15, 20, 15, 15)
    )
}

# === 1. 读取数据 ===
data_path <- public_file("species_data_with_richness.csv")
data <- read.csv(data_path, stringsAsFactors = FALSE)

# === 2. 创建唯一 SampleID ===
data <- data %>% mutate(SampleID = paste(city, habitat, row_number(), sep = "_"))

# === 3. 构建群落矩阵 ===
community_matrix <- data %>%
  group_by(SampleID, species) %>%
  summarise(density = sum(species_density, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = density, values_fill = 0) %>%
  column_to_rownames("SampleID")

# === 4. 元数据 ===
metadata <- data %>% select(SampleID, city, habitat) %>% distinct()

# === 5. 去除全零样本 ===
all_zero_rows <- apply(community_matrix, 1, function(x) all(x == 0))
community_matrix_filtered <- community_matrix[!all_zero_rows, ]
metadata_filtered <- metadata %>% filter(SampleID %in% rownames(community_matrix_filtered))
cat("Removed", sum(all_zero_rows), "samples with zero species.\n")

# === 6. Bray-Curtis 距离矩阵 ===
dist_mat <- vegdist(community_matrix_filtered, method = "bray")

# === 7. PCoA分析 ===
pcoa <- cmdscale(dist_mat, eig = TRUE, k = 2)
pcoa_points <- as.data.frame(pcoa$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points$SampleID <- rownames(pcoa_points)

eig_vals <- pcoa$eig
var_exp <- round(100 * eig_vals[1:2] / sum(eig_vals[eig_vals > 0]), 1)

# === 8. 合并数据 ===
pcoa_data <- pcoa_points %>%
  left_join(metadata_filtered, by = "SampleID") %>%
  mutate(habitat = factor(habitat),
         city = factor(city))

# === 9. PERMANOVA分析 ===
adonis_result <- adonis2(
  community_matrix_filtered ~ city + habitat,
  data = metadata_filtered,
  method = "bray",
  permutations = 999,
  by = "terms"
)
print(adonis_result)

# === 10. PCoA绘图 ===
shape_values <- seq(16, 16 + length(levels(pcoa_data$city)) - 1)

# 指定 habitat 配色
habitat_colors <- c("Pasture" = "#D62728", "Natural grassland" = "#1F77B4")

pcoa_plot <- ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2)) +
  geom_jitter(aes(color = habitat, shape = city),
              width = 0.1, height = 0.1,
              size = 4.5, stroke = 1.2, alpha = 0.6) +
  scale_color_manual(values = habitat_colors) +
  scale_shape_manual(values = shape_values) +
  set_global_fontsize(base_size = 22, base_family = base_family) +   # 一键调整字体
  labs(
    x = paste0("PCoA Axis 1 (", var_exp[1], "% variance)"),
    y = paste0("PCoA Axis 2 (", var_exp[2], "% variance)"),
    color = "Habitat",
    shape = "City"
  )

print(pcoa_plot)

# === 11. 保存PDF ===
ggsave(
  filename = public_output("PCoA_Jittered.pdf"),
  plot = pcoa_plot,
  width = 9,
  height = 7,
  units = "in",
  device = cairo_pdf,
  bg = "transparent"
)











# === 加载必要包 ===
library(tidyverse)
library(vegan)

# === 一键控制函数 ===
set_plot_style <- function(base_size = 22, base_family = "sans") {
  list(
    theme_minimal(base_size = base_size, base_family = base_family) %+replace%
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "#e0e0e0", size = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey70", fill = NA, size = 0.8),
        axis.title  = element_text(face = "bold", size = rel(1.2), color = "#222222"),
        axis.text   = element_text(size = rel(1.0), color = "#444444"),
        axis.ticks  = element_line(color = "#444444", size = 0.7),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = rel(1.1), color = "#222222"),
        legend.text  = element_text(size = rel(1.0), color = "#222222"),
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        legend.spacing.y = unit(0.3, "cm"),
        legend.box.spacing = unit(0.5, "cm"),
        plot.margin = margin(15, 20, 15, 15)
      ),
    # 根据字体大小自动计算图幅
    save_width  = 0.45 * base_size,
    save_height = 0.38 * base_size
  )
}

# === 设置全局参数 ===
plot_cfg <- set_plot_style(base_size = 18)  # 字号加大

# === 1. 读取数据 ===
data_path <- public_file("species_data_with_richness.csv")
data <- read.csv(data_path, stringsAsFactors = FALSE)

# === 2. 创建唯一 SampleID ===
data <- data %>% mutate(SampleID = paste(city, habitat, row_number(), sep = "_"))

# === 3. 构建群落矩阵 ===
community_matrix <- data %>%
  group_by(SampleID, species) %>%
  summarise(density = sum(species_density, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = species, values_from = density, values_fill = 0) %>%
  column_to_rownames("SampleID")

# === 4. 样本元数据 ===
metadata <- data %>% select(SampleID, habitat) %>% distinct()

# === 5. 去除全零样本 ===
all_zero_rows <- apply(community_matrix, 1, function(x) all(x == 0))
community_matrix_filtered <- community_matrix[!all_zero_rows, ]
metadata_filtered <- metadata %>% filter(SampleID %in% rownames(community_matrix_filtered))

# === 6. PERMANOVA 检验 ===
dist_mat <- vegdist(community_matrix_filtered, method = "bray")
adonis_habitat <- adonis2(community_matrix_filtered ~ habitat,
                          data = metadata_filtered,
                          method = "bray",
                          permutations = 999)
p_val <- signif(adonis_habitat$`Pr(>F)`[1], 3)
pval_text <- paste0("PERMANOVA: p = ", p_val)

# === 7. 计算组质心 ===
pcoa <- cmdscale(dist_mat, eig = TRUE, k = nrow(community_matrix_filtered)-1)
scores <- as.data.frame(pcoa$points)
colnames(scores) <- paste0("PCoA", 1:ncol(scores))
scores$SampleID <- rownames(scores)
scores <- scores %>% left_join(metadata_filtered, by = "SampleID")

centroids <- scores %>%
  group_by(habitat) %>%
  summarise(across(starts_with("PCoA"), mean))

# === 8. 计算每个样本到对方组质心的距离 ===
distance_to_other <- scores %>%
  rowwise() %>%
  mutate(
    DistanceToOther = ifelse(
      habitat == "Pasture",
      sqrt(sum((c_across(starts_with("PCoA")) - as.numeric(centroids[centroids$habitat=="Natural grassland",-1]))^2)),
      sqrt(sum((c_across(starts_with("PCoA")) - as.numeric(centroids[centroids$habitat=="Pasture",-1]))^2))
    )
  ) %>%
  ungroup() %>%
  select(SampleID, habitat, DistanceToOther)

# === 9. 箱线图可视化 ===
box_plot <- ggplot(distance_to_other, aes(x = habitat, y = DistanceToOther, fill = habitat)) +
  geom_boxplot(alpha = 0.6, width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 2) +
  scale_fill_manual(values = c("Pasture" = "#D62728", "Natural grassland" = "#1F77B4")) +
  plot_cfg[[1]] +   # 应用字体 & 主题
  labs(
    x = "Habitat",
    y = "Distance to opposite group centroid"
  ) +
  annotate("text",
           x = 1.5,
           y = max(distance_to_other$DistanceToOther) * 1.05,  # 往上移一点
           label = pval_text,
           size = 6, fontface = "bold")

print(box_plot)

# === 10. 保存图像 ===
ggsave(public_output("Pasture_vs_Natural grassland_Between_boxplot.pdf"),
       plot = box_plot,
       width = plot_cfg$save_width,
       height = plot_cfg$save_height,
       units = "in", device = cairo_pdf)



