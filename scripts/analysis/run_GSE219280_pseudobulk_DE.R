script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
project_root <- analysis_project_root()
configure_project_library(project_root)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(edgeR)
  library(limma)
})

input_rds <- file.path(project_root, "data", "processed", "GSE219280_major_celltype_pseudobulk.rds")
minimum_nuclei <- as.integer(Sys.getenv("MINIMUM_NUCLEI", unset = "30"))
result_dir <- if (minimum_nuclei == 30L) {
  file.path(project_root, "results", "GSE219280_pseudobulk_DE")
} else {
  file.path(project_root, "results", sprintf("GSE219280_pseudobulk_DE_min%d", minimum_nuclei))
}
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

minimum_library_size <- 10000
minimum_fold_change <- 1.2

obj <- readRDS(input_rds)
all_counts <- obj$counts
all_meta <- as.data.table(obj$samples)
stopifnot(identical(colnames(all_counts), all_meta$pseudobulk_id))

all_meta[, diagnosis_short := fifelse(
  diagnosis == "C9-ALS", "ALS",
  fifelse(diagnosis == "C9-FTD", "FTD", "Ctrl")
)]
all_meta[, region_short := fifelse(brain_region == "motor cortex", "MCX", "FCX")]
all_meta[, condition_region := paste(diagnosis_short, region_short, sep = "_")]

core_all_disease <- c("Astro", "OPC", "Oligo")
als_focused <- c("Exc_deep", "Exc_intermediate", "Exc_upper", "Micro")

contrast_all_disease <- c(
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

contrast_als_focused <- c(
  ALS_vs_Ctrl_mean = "(ALS_FCX + ALS_MCX - Ctrl_FCX - Ctrl_MCX) / 2",
  ALS_region_interaction = "(ALS_MCX - ALS_FCX) - (Ctrl_MCX - Ctrl_FCX)",
  ALS_vs_Ctrl_MCX = "ALS_MCX - Ctrl_MCX",
  ALS_vs_Ctrl_FCX = "ALS_FCX - Ctrl_FCX"
)

all_results <- list()
all_qc <- list()
all_inclusion <- list()
result_index <- 0L

fit_cell_type <- function(cell_type, diagnoses, required_levels, contrast_definitions, analysis_name) {
  selected <- which(
    all_meta$annotation_major_cell_type == cell_type &
      all_meta$diagnosis_short %chin% diagnoses &
      all_meta$n_nuclei >= minimum_nuclei &
      all_meta$library_size_cellbender >= minimum_library_size
  )
  meta <- copy(all_meta[selected])
  counts <- all_counts[, selected, drop = FALSE]
  stopifnot(identical(colnames(counts), meta$pseudobulk_id))

  present_levels <- sort(unique(meta$condition_region))
  missing_levels <- setdiff(required_levels, present_levels)
  if (length(missing_levels) > 0L) {
    stop(sprintf("%s/%s missing required groups: %s", analysis_name, cell_type, paste(missing_levels, collapse = ", ")))
  }

  meta[, condition_region := factor(condition_region, levels = required_levels)]
  meta[, sex := factor(sex)]
  meta[, sequencing_batch := factor(sequencing_batch)]
  meta[, age_scaled := as.numeric(scale(age))]

  design <- model.matrix(
    ~ 0 + condition_region + sex + age_scaled + sequencing_batch,
    data = meta
  )
  colnames(design) <- sub("^condition_region", "", colnames(design))
  if (qr(design)$rank != ncol(design)) {
    stop(sprintf("Rank-deficient design for %s/%s", analysis_name, cell_type))
  }

  y <- DGEList(counts = counts)
  keep <- filterByExpr(y, group = meta$condition_region, min.count = 10, min.total.count = 15)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")

  v0 <- voom(y, design, plot = FALSE)
  corfit <- duplicateCorrelation(v0, design, block = meta$donor_id)
  consensus <- corfit$consensus.correlation
  if (!is.finite(consensus)) consensus <- 0
  v <- voom(y, design, plot = FALSE, block = meta$donor_id, correlation = consensus)
  fit <- lmFit(v, design, block = meta$donor_id, correlation = consensus)

  contrast_matrix <- makeContrasts(contrasts = unname(contrast_definitions), levels = design)
  colnames(contrast_matrix) <- names(contrast_definitions)
  fit_contrasts <- contrasts.fit(fit, contrast_matrix)
  fit_ebayes <- eBayes(fit_contrasts, robust = TRUE)
  fit_treat <- treat(fit_contrasts, lfc = log2(minimum_fold_change), robust = TRUE)

  output <- list()
  for (contrast_name in colnames(contrast_matrix)) {
    eb <- topTable(fit_ebayes, coef = contrast_name, number = Inf, sort.by = "none", adjust.method = "BH")
    tr <- topTable(fit_treat, coef = contrast_name, number = Inf, sort.by = "none", adjust.method = "BH")
    stopifnot(identical(rownames(eb), rownames(tr)))
    output[[contrast_name]] <- data.table(
      analysis = analysis_name,
      cell_type = cell_type,
      contrast = contrast_name,
      gene = rownames(eb),
      log2FC = eb$logFC,
      AveExpr = eb$AveExpr,
      t = eb$t,
      PValue = eb$P.Value,
      FDR = eb$adj.P.Val,
      treat_t = tr$t,
      treat_PValue = tr$P.Value,
      treat_FDR = tr$adj.P.Val,
      n_samples = nrow(meta),
      n_donors = uniqueN(meta$donor_id),
      duplicate_correlation = consensus,
      minimum_nuclei = minimum_nuclei
    )
  }

  group_counts <- meta[, .(
    n_samples = .N,
    n_donors = uniqueN(donor_id),
    median_nuclei = as.numeric(median(n_nuclei)),
    minimum_nuclei_observed = as.numeric(min(n_nuclei)),
    median_library_size = as.numeric(median(library_size_cellbender))
  ), by = condition_region]
  group_counts[, `:=`(analysis = analysis_name, cell_type = cell_type)]

  qc <- data.table(
    analysis = analysis_name,
    cell_type = cell_type,
    input_genes = nrow(counts),
    retained_genes = nrow(y),
    samples = nrow(meta),
    donors = uniqueN(meta$donor_id),
    design_columns = ncol(design),
    residual_df_nominal = nrow(design) - qr(design)$rank,
    duplicate_correlation = consensus
  )

  list(results = rbindlist(output), qc = qc, inclusion = group_counts)
}

for (cell_type in core_all_disease) {
  fitted <- fit_cell_type(
    cell_type,
    diagnoses = c("Ctrl", "ALS", "FTD"),
    required_levels = c("Ctrl_FCX", "Ctrl_MCX", "ALS_FCX", "ALS_MCX", "FTD_FCX", "FTD_MCX"),
    contrast_definitions = contrast_all_disease,
    analysis_name = "all_disease_core_glia"
  )
  result_index <- result_index + 1L
  all_results[[result_index]] <- fitted$results
  all_qc[[result_index]] <- fitted$qc
  all_inclusion[[result_index]] <- fitted$inclusion
}

for (cell_type in als_focused) {
  fitted <- fit_cell_type(
    cell_type,
    diagnoses = c("Ctrl", "ALS"),
    required_levels = c("Ctrl_FCX", "Ctrl_MCX", "ALS_FCX", "ALS_MCX"),
    contrast_definitions = contrast_als_focused,
    analysis_name = "ALS_focused"
  )
  result_index <- result_index + 1L
  all_results[[result_index]] <- fitted$results
  all_qc[[result_index]] <- fitted$qc
  all_inclusion[[result_index]] <- fitted$inclusion
}

results <- rbindlist(all_results, use.names = TRUE)
results[, global_FDR := p.adjust(PValue, method = "BH")]
results[, global_treat_FDR := p.adjust(treat_PValue, method = "BH")]
setorder(results, analysis, cell_type, contrast, PValue)

summary <- results[, .(
  tested_genes = .N,
  ebayes_FDR_0_05 = sum(FDR < 0.05),
  treat_FDR_0_05 = sum(treat_FDR < 0.05),
  global_FDR_0_05 = sum(global_FDR < 0.05),
  global_treat_FDR_0_05 = sum(global_treat_FDR < 0.05),
  max_abs_log2FC = max(abs(log2FC))
), by = .(analysis, cell_type, contrast)]
setorder(summary, analysis, cell_type, contrast)

fwrite(results, file.path(result_dir, "all_gene_contrast_results.tsv"), sep = "\t", quote = FALSE)
fwrite(summary, file.path(result_dir, "contrast_summary.tsv"), sep = "\t", quote = FALSE)
fwrite(rbindlist(all_qc), file.path(result_dir, "model_QC.tsv"), sep = "\t", quote = FALSE)
fwrite(rbindlist(all_inclusion), file.path(result_dir, "sample_inclusion_by_group.tsv"), sep = "\t", quote = FALSE)

cat(sprintf("models=%d\n", nrow(rbindlist(all_qc))))
cat(sprintf("result_rows=%d\n", nrow(results)))
print(summary[contrast %chin% c("ALS_vs_Ctrl_mean", "FTD_vs_Ctrl_mean", "ALS_vs_FTD_mean", "Shared_disease_vs_Ctrl")])
