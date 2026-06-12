# Author: Wang Linglong
# Public release script derived from: 6\AIC热图.txt
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
library(readr)
library(MuMIn)
library(dplyr)
library(broom)
library(tidyr)
library(stringr)

# === Step 1: 读取数据和构建 full model（用于 dredge） ===
data <- read_csv(public_file("最佳拟合模型_修改后.csv"))
full_model <- glm(Virus_Richness ~ ., data = data, family = gaussian())

# === Step 2: 全子集回归（所有变量组合）===
options(na.action = "na.fail")
model_set <- dredge(full_model, trace = FALSE)

# === Step 3: 提取 ΔAIC 最小的前10个模型 ===
model_info <- as.data.frame(model_set)
model_info$Model_ID <- paste0("Model_", seq_len(nrow(model_info)))
model_info$Delta_AIC <- model_info$AIC - min(model_info$AIC)
model_info <- model_info %>% arrange(Delta_AIC)

top_n <- 10
top_info <- model_info %>% slice(1:top_n)
top_models <- get.models(model_set, subset = as.numeric(rownames(top_info)))

# === Step 4: 定义关注变量顺序（与你论文一致）===
target_vars <- c(
  "host_cluster",
  "livestock_Density",
  "host_body_size",
  "host_density",
  "host_habitat",
  "host_health_status",
  "anthelmintics"
)

# === Step 5: 遍历前10个模型，计算解释率 ===
result_list <- list()

for (i in 1:top_n) {
  model <- top_models[[i]]
  model_row <- top_info[i, ]
  
  model_terms <- attr(terms(model), "term.labels")
  coefs <- tidy(model)
  Df_model <- deviance(model)
  Dn_model <- model$null.deviance
  
  # 模型整体解释度
  model_explained <- round((Dn_model - Df_model) / Dn_model * 100, 2)
  
  # 初始化变量解释率
  row <- setNames(rep(NA, length(target_vars)), target_vars)
  
  for (var in target_vars) {
    related_terms <- coefs %>% filter(str_detect(term, paste0("^", var)))
    
    # 如果该变量显著，则计算其解释度
    if (nrow(related_terms) > 0 && any(related_terms$p.value < 0.05, na.rm = TRUE)) {
      reduced_terms <- setdiff(model_terms, grep(paste0("^", var), model_terms, value = TRUE))
      
      if (length(reduced_terms) == 0) next  # 跳过空模型
      reduced_formula <- as.formula(paste("Virus_Richness ~", paste(reduced_terms, collapse = " + ")))
      reduced_model <- glm(reduced_formula, data = data, family = gaussian())
      Di <- deviance(reduced_model)
      
      row[var] <- round((Di - Df_model) / Dn_model * 100, 2)
    }
  }
  
  result_list[[i]] <- c(
    Model_ID = model_row$Model_ID,
    AIC = round(model_row$AIC, 2),
    Delta_AIC = round(model_row$Delta_AIC, 2),
    Explained_Deviance = model_explained,
    row
  )
}

# === Step 6: 输出整理为表格并保存 ===
final_result <- bind_rows(result_list) %>%
  relocate(Model_ID, AIC, Delta_AIC, Explained_Deviance, all_of(target_vars))

# 打印前10行
print(final_result, n = 10)

# 保存为 CSV 文件
write.csv(final_result, public_output("Top10_Models_ExplainedRate_WithTotal.csv"), row.names = FALSE)





library(mgcv)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# 读取数据
data <- read.csv(public_file("最佳拟合模型_修改后gam.csv"))
data$habitat_overlap <- as.factor(data$habitat_overlap)
data$diet_overlap <- as.factor(data$diet_overlap)

# --- Step 1: 所有变量组合（2^3 - 1 = 7种 + full model = 8）
vars <- c("phylogenetic_distance", "habitat_overlap", "diet_overlap")
all_combinations <- unlist(lapply(1:length(vars), function(n) combn(vars, n, simplify = FALSE)), recursive = FALSE)

# --- Step 2: 构建模型 + 计算 AIC 和解释度
model_list <- list()

for (i in seq_along(all_combinations)) {
  term_set <- all_combinations[[i]]
  
  # 构建公式
  formula_terms <- sapply(term_set, function(var) {
    if (var == "phylogenetic_distance") {
      return("s(phylogenetic_distance, k = 5)")
    } else {
      return(var)
    }
  })
  
  fml <- as.formula(paste("shared_protozoa_count ~", paste(formula_terms, collapse = " + ")))
  
  # 拟合模型
  model <- gam(fml, data = data, family = poisson(link = "log"), method = "REML")
  model_name <- paste0("Model_", i)
  
  # 空模型
  null_model <- gam(shared_protozoa_count ~ 1, data = data, family = poisson(link = "log"), method = "REML")
  Dn <- deviance(null_model)
  Df <- deviance(model)
  
  # 总解释度
  total_exp <- round((Dn - Df) / Dn * 100, 2)
  
  # 每个变量的解释率
  var_exp <- setNames(rep(NA, length(vars)), vars)
  for (var in vars) {
    if (var %in% term_set) {
      sub_terms <- setdiff(term_set, var)
      if (length(sub_terms) == 0) next
      sub_formula_terms <- sapply(sub_terms, function(v) {
        if (v == "phylogenetic_distance") "s(phylogenetic_distance, k = 5)" else v
      })
      fml_sub <- as.formula(paste("shared_protozoa_count ~", paste(sub_formula_terms, collapse = " + ")))
      model_sub <- gam(fml_sub, data = data, family = poisson(link = "log"), method = "REML")
      Di <- deviance(model_sub)
      var_exp[var] <- round((Di - Df) / Dn * 100, 2)
    }
  }
  
  model_list[[i]] <- data.frame(
    Model_ID = model_name,
    AIC = round(AIC(model), 2),
    Explained_Deviance = total_exp,
    phylogenetic_distance = var_exp["phylogenetic_distance"],
    habitat_overlap = var_exp["habitat_overlap"],
    diet_overlap = var_exp["diet_overlap"]
  )
}

# --- Step 3: 合并表格，计算 ΔAIC，排序并保留 Top10
result_df <- bind_rows(model_list)
result_df$Delta_AIC <- result_df$AIC - min(result_df$AIC)
result_df <- result_df %>%
  arrange(Delta_AIC) %>%
  mutate(Model_ID = paste0("Model_", row_number())) %>%
  relocate(Model_ID, AIC, Delta_AIC, Explained_Deviance, everything())

# --- Step 4: 导出
top10 <- result_df %>% slice(1:10)
write.csv(top10, public_output("GAM_Top10_Models_WithTotalExplained.csv"), row.names = FALSE)

# --- 查看结果
print(top10)




