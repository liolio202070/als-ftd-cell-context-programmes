#!/usr/bin/env Rscript

# Donor-region-level cell composition analysis for GSE219280. This is a
# contextual analysis of captured nuclei proportions, not a direct estimate
# of absolute in-vivo cell abundance.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
project_root <- analysis_project_root()
configure_project_library(project_root)
project_lib <- file.path(project_root, "cache", "R_library")
.libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

out_dir <- file.path(project_root, "results", "GSE219280_donor_composition")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
status_file <- file.path(out_dir, "status.tsv")
if (file.exists(status_file)) file.remove(status_file)

record_status <- function(stage, state, detail, progress) {
  row <- data.table(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stage = stage,
    state = state,
    detail = detail,
    progress = progress
  )
  fwrite(row, status_file, sep = "\t", quote = FALSE,
         append = file.exists(status_file), col.names = !file.exists(status_file))
  cat(sprintf("[%s] %s | %s | %s | %.1f%%\n", row$timestamp, stage, state, detail, progress))
  flush.console()
}

record_status("initialize", "started", "Read clean nucleus metadata", 2)
metadata_path <- file.path(project_root, "data", "metadata", "GSE219280_nuclei_metadata.tsv")
meta <- fread(metadata_path)

sample_columns <- c("sample", "donor_id", "diagnosis", "brain_region", "sex", "age", "sequencing_batch")
sample_meta <- unique(meta[, ..sample_columns])
stopifnot(uniqueN(sample_meta$sample) == nrow(sample_meta))

sample_meta[, diagnosis_short := fcase(
  diagnosis == "C9-ALS", "ALS",
  diagnosis == "C9-FTD", "FTD",
  default = "Ctrl"
)]
sample_meta[, region_short := fifelse(brain_region == "motor cortex", "MCX", "FCX")]
sample_meta[, condition_region := paste(diagnosis_short, region_short, sep = "_")]

cell_types <- sort(unique(meta$annotation_major_cell_type))
counts_long <- meta[, .(n_nuclei = .N), by = .(sample, annotation_major_cell_type)]
complete_grid <- CJ(sample = sample_meta$sample, annotation_major_cell_type = cell_types)
counts_long <- counts_long[complete_grid, on = c("sample", "annotation_major_cell_type")]
counts_long[is.na(n_nuclei), n_nuclei := 0L]
counts_long <- merge(counts_long, sample_meta, by = "sample", all.x = TRUE)
counts_long[, total_nuclei := sum(n_nuclei), by = sample]
counts_long[, proportion := n_nuclei / total_nuclei]
counts_long[, arcsin_sqrt_proportion := asin(sqrt((n_nuclei + 0.5) / (total_nuclei + 1)))]

fwrite(counts_long, file.path(out_dir, "sample_celltype_composition.tsv"), sep = "\t", quote = FALSE)
record_status("aggregate", "complete", sprintf("%d samples; %d major cell types", nrow(sample_meta), length(cell_types)), 25)

wide <- dcast(counts_long, annotation_major_cell_type ~ sample, value.var = "arcsin_sqrt_proportion")
composition_matrix <- as.matrix(wide[, -1])
rownames(composition_matrix) <- wide$annotation_major_cell_type
sample_meta <- sample_meta[match(colnames(composition_matrix), sample)]
stopifnot(identical(colnames(composition_matrix), sample_meta$sample))

required_levels <- c("Ctrl_FCX", "Ctrl_MCX", "ALS_FCX", "ALS_MCX", "FTD_FCX", "FTD_MCX")
if (!all(required_levels %chin% sample_meta$condition_region)) stop("A required diagnosis-region group is missing")
sample_meta[, condition_region := factor(condition_region, levels = required_levels)]
sample_meta[, sex := factor(sex)]
sample_meta[, sequencing_batch := factor(sequencing_batch)]
sample_meta[, age_scaled := as.numeric(scale(age))]

design <- model.matrix(~ 0 + condition_region + sex + age_scaled + sequencing_batch, data = sample_meta)
colnames(design) <- sub("^condition_region", "", colnames(design))
if (qr(design)$rank != ncol(design)) stop("Rank-deficient composition design")

contrast_definitions <- c(
  ALS_vs_Ctrl_mean = "(ALS_FCX + ALS_MCX - Ctrl_FCX - Ctrl_MCX) / 2",
  FTD_vs_Ctrl_mean = "(FTD_FCX + FTD_MCX - Ctrl_FCX - Ctrl_MCX) / 2",
  ALS_vs_FTD_mean = "(ALS_FCX + ALS_MCX - FTD_FCX - FTD_MCX) / 2",
  Shared_disease_vs_Ctrl = "(ALS_FCX + ALS_MCX + FTD_FCX + FTD_MCX) / 4 - (Ctrl_FCX + Ctrl_MCX) / 2",
  ALS_region_interaction = "(ALS_MCX - ALS_FCX) - (Ctrl_MCX - Ctrl_FCX)",
  FTD_region_interaction = "(FTD_MCX - FTD_FCX) - (Ctrl_MCX - Ctrl_FCX)",
  ALS_vs_Ctrl_MCX = "ALS_MCX - Ctrl_MCX",
  ALS_vs_Ctrl_FCX = "ALS_FCX - Ctrl_FCX",
  FTD_vs_Ctrl_MCX = "FTD_MCX - Ctrl_MCX",
  FTD_vs_Ctrl_FCX = "FTD_FCX - Ctrl_FCX"
)
contrast_matrix <- makeContrasts(contrasts = unname(contrast_definitions), levels = design)
colnames(contrast_matrix) <- names(contrast_definitions)

record_status("model", "started", "Fit donor-blocked arcsin-sqrt proportion model", 45)
corfit <- duplicateCorrelation(composition_matrix, design, block = sample_meta$donor_id)
consensus <- corfit$consensus.correlation
if (!is.finite(consensus)) consensus <- 0
fit <- lmFit(composition_matrix, design, block = sample_meta$donor_id, correlation = consensus)
fit <- eBayes(contrasts.fit(fit, contrast_matrix), robust = TRUE)

results <- rbindlist(lapply(colnames(contrast_matrix), function(contrast_name) {
  table <- topTable(fit, coef = contrast_name, number = Inf, sort.by = "none", adjust.method = "BH")
  data.table(
    cell_type = rownames(table),
    contrast = contrast_name,
    transformed_effect = table$logFC,
    AveExpr = table$AveExpr,
    t = table$t,
    PValue = table$P.Value,
    comparison_FDR = table$adj.P.Val,
    analysis_role = fifelse(
      contrast_name %chin% c("ALS_vs_Ctrl_mean", "FTD_vs_Ctrl_mean", "ALS_vs_FTD_mean"),
      "primary_composition_context",
      "secondary_region_context"
    ),
    n_samples = nrow(sample_meta),
    n_donors = uniqueN(sample_meta$donor_id),
    residual_df_nominal = nrow(design) - qr(design)$rank,
    duplicate_correlation = consensus
  )
}), use.names = TRUE)

results[, global_primary_FDR := NA_real_]
results[analysis_role == "primary_composition_context", global_primary_FDR := p.adjust(PValue, method = "BH")]

group_summary <- counts_long[, .(
  n_samples = .N,
  mean_proportion = mean(proportion),
  median_proportion = median(proportion),
  minimum_proportion = min(proportion),
  maximum_proportion = max(proportion),
  mean_nuclei = mean(n_nuclei)
), by = .(annotation_major_cell_type, diagnosis_short, region_short)]

celltype_role <- data.table(
  cell_type = cell_types,
  evidence_role = fifelse(cell_types %chin% c("Astro", "OPC", "Oligo"), "predefined_core_glia", "secondary_composition_context")
)
results <- merge(results, celltype_role, by = "cell_type", all.x = TRUE)
setorder(results, analysis_role, contrast, PValue, cell_type)

fwrite(results, file.path(out_dir, "donor_composition_contrast_results.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(group_summary, file.path(out_dir, "composition_group_summary.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(sample_meta, file.path(out_dir, "composition_model_samples.tsv"), sep = "\t", quote = FALSE, na = "")

primary_summary <- results[analysis_role == "primary_composition_context", .(
  cell_types_tested = .N,
  nominal_P_0_05 = sum(PValue < 0.05),
  comparison_FDR_0_05 = sum(comparison_FDR < 0.05),
  global_primary_FDR_0_05 = sum(global_primary_FDR < 0.05),
  core_glia_global_FDR_0_05 = sum(global_primary_FDR < 0.05 & evidence_role == "predefined_core_glia")
), by = contrast]
fwrite(primary_summary, file.path(out_dir, "composition_primary_summary.tsv"), sep = "\t", quote = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

record_status(
  "pipeline", "complete",
  sprintf("%d primary tests; %d global-FDR composition signals", sum(results$analysis_role == "primary_composition_context"), sum(results$global_primary_FDR < 0.05, na.rm = TRUE)),
  100
)

print(primary_summary)
print(results[analysis_role == "primary_composition_context" & (comparison_FDR < 0.05 | global_primary_FDR < 0.05)])
