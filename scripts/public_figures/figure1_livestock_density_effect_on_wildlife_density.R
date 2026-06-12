# Author: Wang Linglong
# Public release script derived from: 1\家畜密度影响野生动物密度.txt
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
整体

library(dplyr)
library(tidyr)

# 1. 读取数据
data <- read.csv(public_file("species_data_with_richness原始未标准化.csv"), 
                 header = TRUE, stringsAsFactors = FALSE)

# 2. 按采样点和物种状态汇总密度
agg_data <- data %>%
  group_by(city, location_id, habitat, species_Order) %>%
  summarise(total_density = sum(species_density, na.rm = TRUE), .groups = "drop")

# 2.1 过滤掉 species_Order 为空的数据
agg_data <- agg_data %>% filter(species_Order != "")

# 3. 将数据转换为宽格式，每个采样点一行，同时包含Wildlife和Livestock的密度
agg_data_wide <- agg_data %>%
  pivot_wider(names_from = species_Order, values_from = total_density)

# 输出结果查看
print(agg_data_wide)
# 相关性检验：检验Livestock与Wildlife密度之间的相关性
cor_test <- cor.test(agg_data_wide$Livestock, agg_data_wide$Wildlife)
print(cor_test)

# 线性回归分析：检验Wildlife密度是否随着Livestock密度的增加而增加
model <- lm(Wildlife ~ Livestock, data = agg_data_wide)
print(summary(model))




单个野生动物
library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)
library(ggpubr)


# 1. 读取数据
data <- read.csv(public_file("species_data_with_richness原始未标准化.csv"), 
                 header = TRUE, stringsAsFactors = FALSE)

# 2. 计算每个采样点的家畜总密度（仅保留species_Order为"Livestock"的数据）
livestock_data <- data %>%
  filter(species_Order == "Livestock") %>%
  group_by(city, location_id, habitat) %>%
  summarise(livestock_density = sum(species_density, na.rm = TRUE), .groups = "drop")

# 3. 对野生动物数据按物种分别计算密度（仅保留species_Order为"Wildlife"的数据）
wildlife_data <- data %>%
  filter(species_Order == "Wildlife") %>%
  group_by(city, location_id, habitat, species) %>%
  summarise(wildlife_density = sum(species_density, na.rm = TRUE), .groups = "drop")

# 4. 将家畜密度与野生动物数据按采样点合并
joined_data <- left_join(wildlife_data, livestock_data, 
                         by = c("city", "location_id", "habitat"))

# 新增步骤：进行回归分析 -------------------------------------------------
# 对每个野生动物物种运行线性回归
results <- joined_data %>%
  group_by(species) %>%
  do(tidy(lm(wildlife_density ~ livestock_density, data = .))) %>%
  filter(term == "livestock_density") %>%  # 只保留家畜密度的系数
  ungroup() %>%
  mutate(p.adjusted = p.adjust(p.value, method = "fdr")) %>%  # 添加FDR校正
  filter(p.value < 0.05)  # 筛选未校正的显著结果

# 5. 提取显著结果并与原始数据合并
significant_species <- results %>% 
  select(species, estimate, p.value) %>%
  mutate(label = paste0("β = ", round(estimate, 2), 
                       "\np = ", format.pval(p.value, digits = 2)))

# 合并显著物种数据
filtered_data <- joined_data %>%
  semi_join(results, by = "species") %>%
  left_join(significant_species, by = "species")

# 创建更美观的可视化并保存为1200dpi
p <- ggplot(filtered_data, aes(x = livestock_density, y = wildlife_density)) +
  geom_point(color = "#4E79A7", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", color = "#E15759", fill = "#F28E2B", 
              alpha = 0.2, linewidth = 1) +
  geom_text(aes(x = Inf, y = Inf, label = label), 
            hjust = 1.1, vjust = 1.5, size = 3,
            data = distinct(filtered_data, species, .keep_all = TRUE)) +
  facet_wrap(~species, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Significant Effects of Livestock Density on Wildlife Species Density",
    subtitle = "Displaying only species with statistically significant relationships (p < 0.05)",
    x = "Livestock Density (individuals/km²)",
    y = "Wildlife Density (individuals/km²)",
    caption = "Shaded areas represent 95% confidence intervals"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold.italic"),
    plot.caption = element_text(color = "gray50"),
    aspect.ratio = 0.8
  )

# Save the plot as 1200 dpi to your desktop
ggsave(public_output("11111plot.png"), plot = p, dpi = 1200, width = 10, height = 8)


