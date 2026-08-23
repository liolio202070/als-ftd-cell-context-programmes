#!/usr/bin/env Rscript

# Leave-one-donor-out stability analysis for the predefined same-cell-type
# mechanism programs. Each omission refits the original donor-blocked
# voom/limma model before recomputing GSEA from gene-level moderated t values.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
project_root <- analysis_project_root()
configure_project_library(project_root)
project_lib <- file.path(project_root, "cache", "R_library")
.libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(edgeR)
  library(limma)
  library(fgsea)
  library(parallel)
})

set.seed(20260811)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

out_dir <- file.path(project_root, "results", "GSE219280_same_celltype_fixed_programs_LOO")
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
  fwrite(
    row, status_file, sep = "\t", quote = FALSE,
    append = file.exists(status_file), col.names = !file.exists(status_file)
  )
  cat(sprintf("[%s] %s | %s | %s | %.1f%%\n", row$timestamp, stage, state, detail, progress))
  flush.console()
}

record_status("initialize", "started", "Load pseudobulk counts and fixed mechanisms", 1)

input_rds <- file.path(project_root, "data", "processed", "GSE219280_major_celltype_pseudobulk.rds")
main_result_path <- file.path(
  project_root, "results", "GSE219280_same_celltype_fixed_programs",
  "fixed_family_fgsea_results.tsv"
)
membership_path <- file.path(
  project_root, "results", "GSE219280_same_celltype_fixed_programs",
  "predefined_gene_set_membership.tsv"
)

obj <- readRDS(input_rds)
all_counts <- obj$counts
all_meta <- as.data.table(obj$samples)
stopifnot(identical(colnames(all_counts), all_meta$pseudobulk_id))

all_meta[, diagnosis_short := fcase(
  diagnosis == "C9-ALS", "ALS",
  diagnosis == "C9-FTD", "FTD",
  default = "Ctrl"
)]
all_meta[, region_short := fifelse(brain_region == "motor cortex", "MCX", "FCX")]
all_meta[, condition_region := paste(diagnosis_short, region_short, sep = "_")]

membership <- unique(fread(membership_path)[, .(family_id, family_label, gene)])
family_pathways <- split(membership$gene, membership$family_id)
family_labels <- unique(membership[, .(set_id = family_id, set_label = family_label)])

main <- fread(main_result_path)[
  threshold == "min30_primary",
  .(
    cell_type, contrast, set_id, set_label,
    full_NES = NES,
    full_pval = pval,
    full_comparison_FDR = comparison_FDR,
    full_global_primary_FDR = global_primary_FDR
  )
]

minimum_nuclei <- 30L
minimum_library_size <- 10000
required_levels <- c("Ctrl_FCX", "Ctrl_MCX", "ALS_FCX", "ALS_MCX", "FTD_FCX", "FTD_MCX")
cell_types <- c("Astro", "OPC", "Oligo")
contrast_definitions <- c(
  ALS_vs_Ctrl_mean = "(ALS_FCX + ALS_MCX - Ctrl_FCX - Ctrl_MCX) / 2",
  FTD_vs_Ctrl_mean = "(FTD_FCX + FTD_MCX - Ctrl_FCX - Ctrl_MCX) / 2",
  ALS_vs_FTD_mean = "(ALS_FCX + ALS_MCX - FTD_FCX - FTD_MCX) / 2"
)

eligible_meta <- all_meta[
  annotation_major_cell_type %chin% cell_types &
    diagnosis_short %chin% c("Ctrl", "ALS", "FTD") &
    n_nuclei >= minimum_nuclei &
    library_size_cellbender >= minimum_library_size
]

tasks <- unique(eligible_meta[, .(
  cell_type = annotation_major_cell_type,
  omitted_donor = donor_id,
  omitted_diagnosis = diagnosis_short
)])
setorder(tasks, cell_type, omitted_diagnosis, omitted_donor)
tasks[, task_id := .I]
fwrite(tasks, file.path(out_dir, "leave_one_donor_out_design.tsv"), sep = "\t", quote = FALSE)

run_fixed_fgsea <- function(ranks) {
  tested <- lapply(family_pathways, function(genes) intersect(genes, names(ranks)))
  tested <- tested[lengths(tested) >= 10L]
  result <- suppressWarnings(fgseaMultilevel(
    pathways = tested,
    stats = ranks,
    minSize = 10L,
    maxSize = 1000L,
    eps = 1e-10,
    scoreType = "std"
  ))
  result <- as.data.table(result)
  result[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
  setnames(result, "pathway", "set_id")
  merge(result, family_labels, by = "set_id", all.x = TRUE)
}

fit_one_omission <- function(task_index) {
  task <- tasks[task_index]
  tryCatch({
    selected <- which(
      all_meta$annotation_major_cell_type == task$cell_type &
        all_meta$diagnosis_short %chin% c("Ctrl", "ALS", "FTD") &
        all_meta$n_nuclei >= minimum_nuclei &
        all_meta$library_size_cellbender >= minimum_library_size &
        all_meta$donor_id != task$omitted_donor
    )
    meta <- copy(all_meta[selected])
    counts <- all_counts[, selected, drop = FALSE]
    stopifnot(identical(colnames(counts), meta$pseudobulk_id))
    if (!all(required_levels %chin% meta$condition_region)) stop("A required diagnosis-region group disappeared")

    meta[, condition_region := factor(condition_region, levels = required_levels)]
    meta[, sex := droplevels(factor(sex))]
    meta[, sequencing_batch := droplevels(factor(sequencing_batch))]
    meta[, age_scaled := as.numeric(scale(age))]
    design <- model.matrix(
      ~ 0 + condition_region + sex + age_scaled + sequencing_batch,
      data = meta
    )
    colnames(design) <- sub("^condition_region", "", colnames(design))
    if (qr(design)$rank != ncol(design)) stop("Rank-deficient design after omission")
    residual_df <- nrow(design) - qr(design)$rank
    if (residual_df < 3L) stop("Fewer than three nominal residual degrees of freedom")

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
    fit <- eBayes(contrasts.fit(fit, contrast_matrix), robust = TRUE)

    results <- rbindlist(lapply(colnames(contrast_matrix), function(contrast_name) {
      gene_table <- topTable(fit, coef = contrast_name, number = Inf, sort.by = "none", adjust.method = "BH")
      ranks <- gene_table$t
      names(ranks) <- rownames(gene_table)
      ranks <- ranks[is.finite(ranks)]
      ranks <- sort(ranks, decreasing = TRUE)
      fg <- run_fixed_fgsea(ranks)
      fg[, `:=`(
        cell_type = task$cell_type,
        contrast = contrast_name,
        omitted_donor = task$omitted_donor,
        omitted_diagnosis = task$omitted_diagnosis,
        ranked_genes = length(ranks),
        n_samples = nrow(meta),
        n_donors = uniqueN(meta$donor_id),
        residual_df_nominal = residual_df,
        duplicate_correlation = consensus
      )]
      fg[, comparison_FDR := p.adjust(pval, method = "BH")]
      fg
    }), use.names = TRUE, fill = TRUE)

    qc <- data.table(
      task_id = task$task_id,
      cell_type = task$cell_type,
      omitted_donor = task$omitted_donor,
      omitted_diagnosis = task$omitted_diagnosis,
      state = "complete",
      error_message = "",
      n_samples = nrow(meta),
      n_donors = uniqueN(meta$donor_id),
      retained_genes = nrow(y),
      residual_df_nominal = residual_df,
      duplicate_correlation = consensus,
      result_rows = nrow(results)
    )
    list(results = results, qc = qc)
  }, error = function(error) {
    list(
      results = data.table(),
      qc = data.table(
        task_id = task$task_id,
        cell_type = task$cell_type,
        omitted_donor = task$omitted_donor,
        omitted_diagnosis = task$omitted_diagnosis,
        state = "error",
        error_message = conditionMessage(error),
        n_samples = NA_integer_, n_donors = NA_integer_, retained_genes = NA_integer_,
        residual_df_nominal = NA_integer_, duplicate_correlation = NA_real_, result_rows = 0L
      )
    )
  })
}

workers <- min(4L, max(1L, detectCores(logical = FALSE)))
batch_size <- workers
batches <- split(seq_len(nrow(tasks)), ceiling(seq_len(nrow(tasks)) / batch_size))
all_results <- list()
all_qc <- list()
result_index <- 0L

record_status(
  "loo_refit",
  "started",
  sprintf("%d cell-type/donor omissions in %d parallel batches", nrow(tasks), length(batches)),
  3
)

for (batch_index in seq_along(batches)) {
  batch_tasks <- batches[[batch_index]]
  batch_results <- mclapply(
    batch_tasks,
    fit_one_omission,
    mc.cores = min(workers, length(batch_tasks)),
    mc.preschedule = FALSE
  )
  for (item in batch_results) {
    result_index <- result_index + 1L
    all_results[[result_index]] <- item$results
    all_qc[[result_index]] <- item$qc
  }
  completed <- min(batch_index * batch_size, nrow(tasks))
  record_status(
    "loo_refit",
    "running",
    sprintf("Completed %d/%d omissions (batch %d/%d)", completed, nrow(tasks), batch_index, length(batches)),
    3 + 87 * completed / nrow(tasks)
  )
}

loo <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
qc <- rbindlist(all_qc, use.names = TRUE, fill = TRUE)
if (any(qc$state != "complete")) {
  fwrite(qc, file.path(out_dir, "leave_one_donor_out_QC.tsv"), sep = "\t", quote = FALSE, na = "")
  stop("One or more leave-one-donor-out fits failed; inspect QC")
}

loo <- merge(
  loo,
  main,
  by = c("cell_type", "contrast", "set_id", "set_label"),
  all.x = TRUE
)
loo[, direction_retained := sign(NES) == sign(full_NES)]
loo[, absolute_delta_NES := abs(NES - full_NES)]
setorder(loo, cell_type, contrast, set_id, omitted_diagnosis, omitted_donor)

fwrite(loo, file.path(out_dir, "leave_one_donor_out_program_results.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(qc, file.path(out_dir, "leave_one_donor_out_QC.tsv"), sep = "\t", quote = FALSE, na = "")

record_status("summarize", "started", "Summarize direction and FDR stability", 93)

stability <- loo[, {
  influence_index <- which.max(absolute_delta_NES)
  .(
    leaveouts_tested = .N,
    full_NES = first(full_NES),
    full_pval = first(full_pval),
    full_comparison_FDR = first(full_comparison_FDR),
    full_global_primary_FDR = first(full_global_primary_FDR),
    direction_retained_fraction = mean(direction_retained),
    nominal_P_0_05_fraction = mean(pval < 0.05),
    comparison_FDR_0_05_fraction = mean(comparison_FDR < 0.05),
    median_LOO_NES = median(NES),
    minimum_LOO_NES = min(NES),
    maximum_LOO_NES = max(NES),
    max_absolute_delta_NES = max(absolute_delta_NES),
    most_influential_donor = omitted_donor[influence_index],
    most_influential_diagnosis = omitted_diagnosis[influence_index]
  )
}, by = .(cell_type, contrast, set_id, set_label)]

comparison_summary <- stability[, .(
  programs = .N,
  minimum_direction_retained_fraction = min(direction_retained_fraction),
  median_direction_retained_fraction = median(direction_retained_fraction),
  full_comparison_FDR_0_05 = sum(full_comparison_FDR < 0.05),
  always_comparison_FDR_0_05 = sum(comparison_FDR_0_05_fraction == 1),
  median_comparison_FDR_0_05_fraction = median(comparison_FDR_0_05_fraction),
  maximum_absolute_delta_NES = max(max_absolute_delta_NES)
), by = .(cell_type, contrast)]

fwrite(stability, file.path(out_dir, "leave_one_donor_out_program_stability.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(comparison_summary, file.path(out_dir, "leave_one_donor_out_comparison_summary.tsv"), sep = "\t", quote = FALSE, na = "")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

record_status(
  "pipeline",
  "complete",
  sprintf(
    "%d/%d fits complete; minimum direction retention %.3f",
    sum(qc$state == "complete"), nrow(qc), min(stability$direction_retained_fraction)
  ),
  100
)

print(comparison_summary)
print(stability[order(direction_retained_fraction, comparison_FDR_0_05_fraction, -max_absolute_delta_NES)])
