#!/usr/bin/env Rscript

# Leave-one-donor-out sensitivity for the independent GSE163122 FTD-GRN
# validation. Refit the donor-blocked model and recompute the seven fixed
# mechanism programs for the region-averaged FTD-versus-control contrast.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
project_root <- analysis_project_root()
configure_project_library(project_root)
project_lib <- file.path(project_root, "cache", "R_library")
.libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(fgsea)
  library(parallel)
})

set.seed(20260811)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

raw_dir <- file.path(project_root, "data", "raw", "GSE163122_FTD_GRN")
main_dir <- file.path(project_root, "results", "GSE163122_FTD_GRN_validation")
out_dir <- file.path(main_dir, "leave_one_donor_out")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
status_file <- file.path(out_dir, "status.tsv")
if (file.exists(status_file)) file.remove(status_file)

record_status <- function(stage, state, detail, progress) {
  row <- data.table(timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                    stage = stage, state = state, detail = detail, progress = progress)
  fwrite(row, status_file, sep = "\t", quote = FALSE,
         append = file.exists(status_file), col.names = !file.exists(status_file))
  cat(sprintf("[%s] %s | %s | %s | %.1f%%\n", row$timestamp, stage, state, detail, progress))
  flush.console()
}

record_status("initialize", "started", "Load independent FTD matrices and metadata", 2)
metadata <- fread(file.path(main_dir, "GSE163122_sample_metadata.tsv"))
main <- fread(file.path(main_dir, "fixed_family_fgsea_results.tsv"))[
  contrast == "FTD_vs_Ctrl_mean",
  .(population, set_id, set_label, full_NES = NES, full_pval = pval,
    full_comparison_FDR = comparison_FDR, full_global_primary_FDR = global_primary_FDR)
]
membership <- unique(fread(file.path(
  project_root, "results", "GSE219280_same_celltype_fixed_programs",
  "predefined_gene_set_membership.tsv"
))[, .(family_id, family_label, gene)])
pathways <- split(membership$gene, membership$family_id)
labels <- unique(membership[, .(set_id = family_id, set_label = family_label)])

required_levels <- c(
  "Ctrl_Frontal", "Ctrl_Temporal", "Ctrl_Occipital",
  "FTD_Frontal", "FTD_Temporal", "FTD_Occipital"
)
contrast_definition <- "(FTD_Frontal + FTD_Temporal + FTD_Occipital - Ctrl_Frontal - Ctrl_Temporal - Ctrl_Occipital) / 3"

read_population <- function(population_value) {
  suffix <- ifelse(population_value == "OLIG2pos", "olig", "neun")
  input <- as.data.table(read.csv(gzfile(file.path(raw_dir, sprintf("GSE163122_countmtx_%s.csv.gz", suffix))), check.names = FALSE))
  counts <- as.matrix(input[, -1L])
  storage.mode(counts) <- "integer"
  rownames(counts) <- input[[1L]]
  meta <- copy(metadata[population == population_value])
  meta <- meta[match(colnames(counts), sample_title)]
  list(counts = counts, meta = meta)
}
objects <- list(NEUNpos = read_population("NEUNpos"), OLIG2pos = read_population("OLIG2pos"))

run_fgsea <- function(ranks) {
  tested <- lapply(pathways, function(genes) intersect(genes, names(ranks)))
  tested <- tested[lengths(tested) >= 10L]
  result <- as.data.table(suppressWarnings(fgseaMultilevel(
    pathways = tested, stats = ranks, minSize = 10L, maxSize = 1000L,
    eps = 1e-10, scoreType = "std"
  )))
  setnames(result, "pathway", "set_id")
  result[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1L))]
  merge(result, labels, by = "set_id", all.x = TRUE)
}

tasks <- rbindlist(lapply(names(objects), function(population_value) {
  unique(objects[[population_value]]$meta[, .(
    population = population_value,
    omitted_donor = donor_id,
    omitted_diagnosis = diagnosis
  )])
}))
setorder(tasks, population, omitted_diagnosis, omitted_donor)
tasks[, task_id := .I]
fwrite(tasks, file.path(out_dir, "LOO_design.tsv"), sep = "\t", quote = FALSE)

fit_task <- function(index) {
  task <- tasks[index]
  tryCatch({
    object <- objects[[task$population]]
    keep_samples <- object$meta$donor_id != task$omitted_donor
    meta <- copy(object$meta[keep_samples])
    counts <- object$counts[, keep_samples, drop = FALSE]
    if (!all(required_levels %chin% meta$condition_region)) stop("Required group disappeared")
    meta[, condition_region := factor(condition_region, levels = required_levels)]
    meta[, sex := droplevels(factor(sex))]
    meta[, age_z := as.numeric(scale(age))]
    meta[, rin_imputed := fifelse(is.na(rin), median(rin, na.rm = TRUE), rin)]
    meta[, rin_z := as.numeric(scale(rin_imputed))]
    design <- model.matrix(~ 0 + condition_region + sex + age_z + rin_z, data = meta)
    colnames(design) <- sub("^condition_region", "", colnames(design))
    if (qr(design)$rank != ncol(design)) stop("Rank-deficient design")

    y <- DGEList(counts = counts)
    keep_genes <- filterByExpr(y, group = meta$condition_region, min.count = 5, min.total.count = 15)
    y <- normLibSizes(y[keep_genes, , keep.lib.sizes = FALSE], method = "TMM")
    v0 <- voom(y, design, plot = FALSE)
    corfit <- duplicateCorrelation(v0, design, block = meta$donor_id)
    consensus <- corfit$consensus.correlation
    if (!is.finite(consensus)) consensus <- 0
    v <- voom(y, design, plot = FALSE, block = meta$donor_id, correlation = consensus)
    fit <- lmFit(v, design, block = meta$donor_id, correlation = consensus)
    contrast <- makeContrasts(contrasts = contrast_definition, levels = design)
    fit <- eBayes(contrasts.fit(fit, contrast), robust = TRUE)
    table <- as.data.table(topTable(fit, number = Inf, sort.by = "none"), keep.rownames = "gene_id")
    symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = table$gene_id, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
    table[, gene := unname(symbols[gene_id])]
    table <- table[!is.na(gene) & nzchar(gene) & is.finite(t)]
    table[, abs_t := abs(t)]
    setorder(table, gene, -AveExpr, -abs_t)
    table <- table[, .SD[1L], by = gene]
    setorder(table, -t, gene)
    ranks <- table$t
    names(ranks) <- table$gene
    result <- run_fgsea(ranks)
    result[, `:=`(
      population = task$population,
      omitted_donor = task$omitted_donor,
      omitted_diagnosis = task$omitted_diagnosis,
      ranked_genes = length(ranks),
      n_samples = nrow(meta), n_donors = uniqueN(meta$donor_id),
      duplicate_correlation = consensus
    )]
    result[, comparison_FDR := p.adjust(pval, method = "BH")]
    qc <- data.table(task_id = task$task_id, population = task$population,
                     omitted_donor = task$omitted_donor, state = "complete",
                     error_message = "", result_rows = nrow(result))
    list(result = result, qc = qc)
  }, error = function(error) list(
    result = data.table(),
    qc = data.table(task_id = task$task_id, population = task$population,
                    omitted_donor = task$omitted_donor, state = "error",
                    error_message = conditionMessage(error), result_rows = 0L)
  ))
}

available_cores <- detectCores(logical = FALSE)
if (is.na(available_cores)) available_cores <- 1L
workers <- min(4L, max(1L, available_cores))
batches <- split(seq_len(nrow(tasks)), ceiling(seq_len(nrow(tasks)) / workers))
items <- list()
item_index <- 0L
record_status("LOO", "started", sprintf("%d fits in %d batches", nrow(tasks), length(batches)), 5)
for (batch_index in seq_along(batches)) {
  results <- mclapply(batches[[batch_index]], fit_task,
                      mc.cores = min(workers, length(batches[[batch_index]])), mc.preschedule = FALSE)
  for (result in results) {
    item_index <- item_index + 1L
    items[[item_index]] <- result
  }
  completed <- min(batch_index * workers, nrow(tasks))
  record_status("LOO", "running", sprintf("Completed %d/%d fits", completed, nrow(tasks)), 5 + 85 * completed / nrow(tasks))
}

qc <- rbindlist(lapply(items, `[[`, "qc"), use.names = TRUE, fill = TRUE)
if (any(qc$state != "complete")) {
  fwrite(qc, file.path(out_dir, "LOO_QC.tsv"), sep = "\t", quote = FALSE, na = "")
  stop("One or more independent FTD LOO fits failed")
}
loo <- rbindlist(lapply(items, `[[`, "result"), use.names = TRUE, fill = TRUE)
loo <- merge(loo, main, by = c("population", "set_id", "set_label"), all.x = TRUE)
loo[, direction_retained := sign(NES) == sign(full_NES)]
loo[, absolute_delta_NES := abs(NES - full_NES)]

stability <- loo[, .(
  leaveouts = .N,
  full_NES = data.table::first(full_NES),
  full_pval = data.table::first(full_pval),
  full_comparison_FDR = data.table::first(full_comparison_FDR),
  full_global_primary_FDR = data.table::first(full_global_primary_FDR),
  direction_retained_fraction = mean(direction_retained),
  comparison_FDR_0_05_fraction = mean(comparison_FDR < 0.05),
  median_LOO_NES = median(NES), min_LOO_NES = min(NES), max_LOO_NES = max(NES),
  max_absolute_delta_NES = max(absolute_delta_NES)
), by = .(population, set_id, set_label)]

fwrite(qc, file.path(out_dir, "LOO_QC.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(loo, file.path(out_dir, "LOO_program_results.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(stability, file.path(out_dir, "LOO_program_stability.tsv"), sep = "\t", quote = FALSE, na = "")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
record_status("pipeline", "complete", sprintf("%d/%d fits complete", sum(qc$state == "complete"), nrow(qc)), 100)
print(stability[order(full_global_primary_FDR, population, set_label)])
