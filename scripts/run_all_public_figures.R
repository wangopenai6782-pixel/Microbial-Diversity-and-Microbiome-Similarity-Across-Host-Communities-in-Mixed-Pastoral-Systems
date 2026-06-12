#!/usr/bin/env Rscript
# Author: Wang Linglong
# Run all public R figure scripts in lexical order.
# Usage:
#   Rscript scripts/run_all_public_figures.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_arg <- args[grep(paste0("^", file_arg), args)]
script_path <- if (length(script_arg) > 0) sub(paste0("^", file_arg), "", script_arg[1]) else file.path(getwd(), "scripts", "run_all_public_figures.R")
root <- normalizePath(file.path(dirname(normalizePath(script_path, mustWork = FALSE)), ".."), mustWork = FALSE)
scripts_dir <- file.path(root, "scripts", "public_figures")
output_dir <- Sys.getenv("JEV_FIGURE_OUTPUT_DIR", file.path(root, "results", "public_figures"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

scripts <- list.files(scripts_dir, pattern = "\\.R$", full.names = TRUE)
message("Running ", length(scripts), " R scripts from ", scripts_dir)
for (script in scripts) {
  message("\n==> ", basename(script))
  tryCatch(
    source(script, local = new.env(parent = globalenv()), encoding = "UTF-8"),
    error = function(e) {
      message("FAILED: ", basename(script), "\n", conditionMessage(e))
    }
  )
}
