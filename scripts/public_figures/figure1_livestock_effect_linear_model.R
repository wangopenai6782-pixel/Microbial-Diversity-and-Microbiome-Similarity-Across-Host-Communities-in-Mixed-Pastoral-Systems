# Author: Wang Linglong
# Public release script derived from: 1\家畜影响其他密度新.txt
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
library(ggplot2)
library(dplyr)
library(broom)
library(stringr)

# ============ 清洗物种名 ============
clean_species <- function(x) {
  x <- str_squish(x)
  x <- str_replace_all(x, "\\bspp\\.?\\b", "sp.")
  x
}

# ============ 数据准备 ============
data <- read.csv(public_file("species_data_with_richness原始未标准化.csv"),
                 header = TRUE, stringsAsFactors = FALSE)

# 家畜密度
livestock_data <- data %>%
  filter(species_Order == "Livestock") %>%
  group_by(city, location_id, habitat) %>%
  summarise(livestock_density = sum(species_density, na.rm = TRUE), .groups = "drop")

# 野生动物密度
wildlife_data <- data %>%
  filter(species_Order == "Wildlife") %>%
  group_by(city, location_id, habitat, species) %>%
  summarise(wildlife_density = sum(species_density, na.rm = TRUE), .groups = "drop")

# 合并
joined_data <- left_join(wildlife_data, livestock_data,
                         by = c("city", "location_id", "habitat")) %>%
  mutate(species = clean_species(species))

# ============ 回归结果 (未校正) ============
model_stats <- joined_data %>%
  group_by(species) %>%
  do({
    m <- lm(wildlife_density ~ livestock_density, data = .)
    tibble(
      estimate = coef(m)[["livestock_density"]],
      p.value  = summary(m)$coefficients["livestock_density", "Pr(>|t|)"]
    )
  }) %>%
  ungroup() %>%
  filter(p.value < 0.05)  # 仅保留显著物种

# 选代表性物种（示例取前6个）
top_species <- model_stats %>%
  arrange(p.value) %>%
  slice_head(n = 6) %>%
  pull(species)

plot_data <- joined_data %>% filter(species %in% top_species)

# 左上角标签（仅显示 p 值）
labels_df <- model_stats %>%
  filter(species %in% top_species) %>%
  transmute(
    species,
    x = -Inf, y = Inf,
    label = ifelse(p.value < 0.001,
                   "p < 0.001",
                   paste0("p = ", format.pval(p.value, digits = 2)))
  )

# 坐标范围（统一）
x_lim <- range(joined_data$livestock_density, na.rm = TRUE)
y_lim <- range(joined_data$wildlife_density, na.rm = TRUE)

# ============ 全局字体大小设置 ============
global_font_size <- 26  # 修改这个值即可全局调整字体 

# ============ 绘图 ============
p <- ggplot(plot_data, aes(x = livestock_density, y = wildlife_density)) +
  geom_point(size = 3, alpha = 0.75, color = "black") +
  geom_smooth(method = "lm", se = TRUE,
              linewidth = 1.2, color = "#E15759", fill = "#E15759", alpha = 0.25) +
  facet_wrap(~ species, scales = "fixed", ncol = 3) +
  geom_text(
    data = labels_df,
    aes(x = x, y = y, label = label),
    hjust = -0.05, vjust = 1.2, size = global_font_size / 3.5, color = "black",
    inherit.aes = FALSE
  ) +
  scale_x_continuous(limits = x_lim, expand = expansion(mult = c(0.06, 0.02))) +
  scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0.10, 0.02))) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = global_font_size) +  # 统一全局字体
  labs(
    x = "Livestock density (individuals/km²)",
    y = "Wildlife density (individuals/km²)"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "italic", size = rel(0.85)),  # 相对大小
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = rel(0.85)),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(10, 10, 10, 14)
  )

# ============ 保存为 PDF ============
ggsave(public_output("nature_style_linear_LEFT_noBeta_FDR_globalFont.pdf"),
       plot = p,
       device = cairo_pdf,
       width = 12, height = 8, units = "in")


