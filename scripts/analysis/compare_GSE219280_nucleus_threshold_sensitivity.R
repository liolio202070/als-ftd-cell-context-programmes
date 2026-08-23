script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
project_root <- analysis_project_root()
configure_project_library(project_root)
project_lib <- file.path(project_root, "cache/R_library")
.libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages(library(data.table))

thresholds <- c(20L, 30L, 50L)
result_path <- function(threshold) {
  result_dir <- if (threshold == 30L) {
    file.path(project_root, "results/GSE219280_pseudobulk_DE")
  } else {
    file.path(project_root, sprintf("results/GSE219280_pseudobulk_DE_min%d", threshold))
  }
  file.path(result_dir, "all_gene_contrast_results.tsv")
}

paths <- vapply(thresholds, result_path, character(1))
stopifnot(all(file.exists(paths)))

all_results <- rbindlist(lapply(seq_along(thresholds), function(i) {
  x <- fread(paths[[i]])
  x[, threshold := thresholds[[i]]]
  x
}), use.names = TRUE)

primary <- all_results[threshold == 30L, .(
  analysis,
  cell_type,
  contrast,
  gene,
  primary_log2FC = log2FC,
  primary_PValue = PValue,
  primary_FDR = FDR,
  primary_treat_FDR = treat_FDR
)]

comparison <- merge(
  all_results[threshold != 30L, .(
    threshold,
    analysis,
    cell_type,
    contrast,
    gene,
    sensitivity_log2FC = log2FC,
    sensitivity_PValue = PValue,
    sensitivity_FDR = FDR,
    sensitivity_treat_FDR = treat_FDR
  )],
  primary,
  by = c("analysis", "cell_type", "contrast", "gene"),
  all = FALSE
)
comparison[, same_direction := sign(primary_log2FC) == sign(sensitivity_log2FC)]

safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  unname(cor(x[ok], y[ok], method = method))
}

summary <- comparison[, .(
  overlapping_genes = .N,
  pearson_effect = safe_cor(primary_log2FC, sensitivity_log2FC, "pearson"),
  spearman_effect = safe_cor(primary_log2FC, sensitivity_log2FC, "spearman"),
  direction_concordance = mean(same_direction, na.rm = TRUE),
  primary_FDR_0_05 = sum(primary_FDR < 0.05),
  sensitivity_FDR_0_05 = sum(sensitivity_FDR < 0.05),
  jointly_FDR_0_05 = sum(primary_FDR < 0.05 & sensitivity_FDR < 0.05),
  primary_treat_FDR_0_05 = sum(primary_treat_FDR < 0.05),
  sensitivity_treat_FDR_0_05 = sum(sensitivity_treat_FDR < 0.05),
  jointly_treat_FDR_0_05 = sum(primary_treat_FDR < 0.05 & sensitivity_treat_FDR < 0.05)
), by = .(threshold, analysis, cell_type, contrast)]
setorder(summary, threshold, analysis, cell_type, contrast)

out_dir <- file.path(project_root, "results/GSE219280_threshold_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(
  comparison,
  file.path(out_dir, "gene_level_threshold_comparison.tsv"),
  sep = "\t", quote = FALSE, na = ""
)
fwrite(
  summary,
  file.path(out_dir, "threshold_sensitivity_summary.tsv"),
  sep = "\t", quote = FALSE, na = ""
)

print(summary[contrast %chin% c(
  "ALS_vs_Ctrl_mean",
  "FTD_vs_Ctrl_mean",
  "ALS_vs_FTD_mean",
  "Shared_disease_vs_Ctrl",
  "ALS_region_interaction",
  "FTD_region_interaction"
)])
