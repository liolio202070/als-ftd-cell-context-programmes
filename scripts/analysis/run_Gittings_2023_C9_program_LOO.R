#!/usr/bin/env Rscript

# Leave-one-donor-out direction stability for the three pre-defined frontal
# contrasts and seven fixed mechanism families in Gittings 2023.  Nuclei are
# never treated as replicates; each iteration removes an entire donor.

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
  library(rhdf5)
  library(fgsea)
})

set.seed(20260811)
out_dir <- file.path(project_root, "results", "Gittings_2023_C9_snRNA")
status_file <- file.path(out_dir, "status.tsv")

record_status <- function(stage, state, detail, progress) {
  row <- data.table(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stage = stage, state = state, detail = detail, progress = progress
  )
  fwrite(row, status_file, sep = "\t", quote = FALSE,
         append = file.exists(status_file), col.names = !file.exists(status_file))
  cat(sprintf("[%s] %s | %s | %s | %.1f%%\n",
              row$timestamp, stage, state, detail, progress))
  flush.console()
}

primary_definitions <- list(
  C9_ALS_vs_Ctrl_Frontal = c(
    ALSPart_ALS_Frontal = 1, ALSPart_Ctrl_Frontal = -1
  ),
  C9_FTD_vs_Ctrl_Frontal = c(
    FTDPart_FTD_Frontal = 1, FTDPart_Ctrl_Frontal = -1
  ),
  C9_ALS_vs_C9_FTD_Frontal_batch_bridged = c(
    ALSPart_ALS_Frontal = 1, ALSPart_Ctrl_Frontal = -1,
    FTDPart_FTD_Frontal = -1, FTDPart_Ctrl_Frontal = 1
  )
)

make_contrasts <- function(definitions, columns) {
  if (!all(vapply(definitions, function(x) all(names(x) %chin% columns), logical(1)))) {
    return(NULL)
  }
  output <- matrix(0, nrow = length(columns), ncol = length(definitions),
                   dimnames = list(columns, names(definitions)))
  for (name in names(definitions)) output[names(definitions[[name]]), name] <- definitions[[name]]
  output
}

membership <- unique(fread(file.path(
  out_dir, "predefined_7_family_gene_membership.tsv"
))[, .(family_id, family_label, gene)])
pathways <- split(membership$gene, membership$family_id)
pathways <- lapply(pathways, unique)
labels <- unique(membership[, .(set_id = family_id, set_label = family_label)])
if (length(pathways) != 7L) stop("Expected seven pre-defined pathways")

run_fgsea <- function(ranks) {
  local <- lapply(pathways, function(x) intersect(x, names(ranks)))
  local <- local[lengths(local) >= 10L]
  output <- as.data.table(suppressWarnings(fgseaMultilevel(
    pathways = local, stats = ranks, minSize = 10L, maxSize = 1000L,
    eps = 1e-10, scoreType = "std", nproc = 1L
  )))
  output[, leadingEdge := vapply(leadingEdge, paste, collapse = ";",
                                  FUN.VALUE = character(1))]
  setnames(output, "pathway", "set_id")
  merge(output, labels, by = "set_id", all.x = TRUE)
}

record_status("program_LOO", "started",
              "Leave one whole donor out; fixed 3x7 frontal family tests", 89)
metadata <- fread(file.path(out_dir, "pseudobulk_sample_metadata.tsv"))
counts_file <- file.path(out_dir, "pseudobulk_raw_counts.h5")
counts <- h5read(counts_file, "counts")
sample_ids <- enc2utf8(as.character(h5read(counts_file, "sample_id")))
gene_ids <- enc2utf8(as.character(h5read(counts_file, "gene_id")))
gene_symbols <- enc2utf8(as.character(h5read(counts_file, "gene_symbol")))
if (nrow(counts) == length(sample_ids) && ncol(counts) == length(gene_ids)) {
  # rhdf5 reverses the dimension order of datasets written by h5py.
  counts <- t(counts)
}
stopifnot(nrow(counts) == length(gene_ids), ncol(counts) == length(sample_ids))
metadata <- metadata[match(sample_ids, sample_id)]
if (anyNA(metadata$sample_id) || !identical(metadata$sample_id, sample_ids)) {
  stop("Pseudobulk sample order mismatch")
}
rownames(counts) <- make.unique(gene_ids)
colnames(counts) <- sample_ids
metadata[, region_model := fifelse(tolower(region) == "frontal", "Frontal", "Occipital")]
metadata[, partition_model := fifelse(grepl("ALS_ALSFTD", partition), "ALSPart", "FTDPart")]
metadata[, diagnosis_model := fcase(
  diagnosis == "Control", "Ctrl", diagnosis == "C9_ALS", "ALS",
  diagnosis == "C9_ALS_FTD", "ALSFTD", diagnosis == "C9_FTD", "FTD",
  default = NA_character_
)]
metadata[, group := paste(partition_model, diagnosis_model, region_model, sep = "_")]
metadata[, sex := factor(sex)]
metadata[, age_z := as.numeric(scale(age))]
metadata <- metadata[n_nuclei >= 30]

full <- fread(file.path(out_dir, "fixed_7_family_fgsea_results.tsv"))[
  analysis_role == "primary_frontal",
  .(cell_type, contrast, set_id, set_label, full_NES = NES,
    full_pval = pval, full_comparison_FDR = comparison_FDR,
    full_global_primary_FDR = global_primary_FDR)
]
model_cell_types <- unique(fread(file.path(out_dir, "program_model_QC.tsv"))$cell_type)
donor_universe <- sort(unique(metadata$donor_id))
expected_rows <- nrow(full)

fit_one <- function(cell_type_value, omitted_donor, definitions) {
  meta <- copy(metadata[cell_type == cell_type_value & donor_id != omitted_donor])
  meta <- meta[match(colnames(counts), sample_id, nomatch = 0L)]
  count_subset <- counts[, meta$sample_id, drop = FALSE]
  meta[, group := factor(group)]
  design <- model.matrix(~ 0 + group + sex + age_z, data = meta)
  colnames(design) <- sub("^group", "", colnames(design))
  if (qr(design)$rank != ncol(design)) {
    return(list(error = "rank_deficient_design"))
  }
  contrast_matrix <- make_contrasts(definitions, colnames(design))
  if (is.null(contrast_matrix)) return(list(error = "primary_group_missing"))
  y <- DGEList(counts = count_subset)
  keep <- filterByExpr(y, design = design, min.count = 5, min.total.count = 15)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y, method = "TMM")
  v0 <- voom(y, design, plot = FALSE)
  correlation_fit <- duplicateCorrelation(v0, design, block = meta$donor_id)
  consensus <- correlation_fit$consensus.correlation
  if (!is.finite(consensus)) consensus <- 0
  v <- voom(y, design, plot = FALSE, block = meta$donor_id, correlation = consensus)
  fit <- lmFit(v, design, block = meta$donor_id, correlation = consensus)
  fit <- eBayes(contrasts.fit(fit, contrast_matrix), robust = TRUE)
  symbol_lookup <- gene_symbols[match(rownames(fit$coefficients), gene_ids)]
  result <- list()
  for (contrast_name in colnames(contrast_matrix)) {
    statistic <- fit$t[, contrast_name]
    rank_table <- data.table(gene = symbol_lookup, statistic = statistic,
                             ave = fit$Amean, abs_t = abs(statistic))
    rank_table <- rank_table[!is.na(gene) & nzchar(gene) & is.finite(statistic)]
    setorder(rank_table, gene, -ave, -abs_t)
    rank_table <- rank_table[, .SD[1L], by = gene]
    setorder(rank_table, -statistic, gene)
    ranks <- rank_table$statistic
    names(ranks) <- rank_table$gene
    fg <- run_fgsea(ranks)
    fg[, `:=`(
      cell_type = cell_type_value, contrast = contrast_name,
      omitted_donor = omitted_donor, refit = TRUE,
      n_samples = nrow(meta), n_donors = uniqueN(meta$donor_id),
      duplicate_correlation = consensus
    )]
    fg[, comparison_FDR := p.adjust(pval, method = "BH")]
    result[[contrast_name]] <- fg
  }
  list(result = rbindlist(result, use.names = TRUE, fill = TRUE), error = "")
}

all_iterations <- list()
qc_iterations <- list()
iteration_index <- 0L

run_cell_for_donor <- function(cell_type_value, omitted) {
  cell_contrasts <- unique(full[cell_type == cell_type_value, contrast])
  cell_definitions <- primary_definitions[names(primary_definitions) %chin% cell_contrasts]
  if (!length(cell_definitions)) {
    return(list(cell_type = cell_type_value, result = NULL,
                error = "no_full_primary_axis", refit = FALSE))
  }
  if (!omitted %chin% metadata[cell_type == cell_type_value, donor_id]) {
    cloned <- copy(full[cell_type == cell_type_value])
    setnames(cloned, c("full_NES", "full_pval", "full_comparison_FDR"),
             c("NES", "pval", "comparison_FDR"))
    cloned[, `:=`(
      omitted_donor = omitted, refit = FALSE,
      n_samples = metadata[cell_type == cell_type_value, .N],
      n_donors = uniqueN(metadata[cell_type == cell_type_value, donor_id]),
      duplicate_correlation = NA_real_, leadingEdge = "", size = NA_integer_,
      ES = NA_real_, log2err = NA_real_, padj = NA_real_
    )]
    cloned[, c("full_global_primary_FDR") := NULL]
    return(list(cell_type = cell_type_value, result = cloned,
                error = "", refit = FALSE))
  }
  fitted <- fit_one(cell_type_value, omitted, cell_definitions)
  list(
    cell_type = cell_type_value,
    result = if (nzchar(fitted$error)) NULL else fitted$result,
    error = fitted$error,
    refit = !nzchar(fitted$error)
  )
}

for (donor_index in seq_along(donor_universe)) {
  omitted <- donor_universe[donor_index]
  workers <- parallel::mclapply(
    model_cell_types, run_cell_for_donor, omitted = omitted,
    mc.cores = min(7L, length(model_cell_types)), mc.preschedule = FALSE,
    mc.set.seed = TRUE
  )
  donor_results <- lapply(workers, `[[`, "result")
  donor_results <- donor_results[!vapply(donor_results, is.null, logical(1))]
  errors <- vapply(workers, function(x) {
    if (nzchar(x$error) && x$error != "no_full_primary_axis")
      paste(x$cell_type, x$error, sep = ":") else ""
  }, character(1))
  errors <- errors[nzchar(errors)]
  donor_table <- rbindlist(donor_results, use.names = TRUE, fill = TRUE)
  complete_iteration <- nrow(donor_table) == expected_rows &&
    !anyDuplicated(donor_table[, .(cell_type, contrast, set_id)])
  donor_table[, LOO_global_primary_FDR := if (complete_iteration) p.adjust(pval, "BH") else NA_real_]
  all_iterations[[omitted]] <- donor_table
  qc_iterations[[omitted]] <- data.table(
    omitted_donor = omitted, rows = nrow(donor_table), expected_rows = expected_rows,
    complete_iteration = complete_iteration, errors = paste(errors, collapse = ";"),
    refitted_cell_types = sum(vapply(workers, `[[`, logical(1), "refit"))
  )
  iteration_index <- iteration_index + 1L
  record_status(
    "program_LOO", "running",
    sprintf("donor %d/%d=%s;rows=%d/%d;complete=%s",
            donor_index, length(donor_universe), omitted, nrow(donor_table),
            expected_rows, complete_iteration),
    89 + 10 * donor_index / length(donor_universe)
  )
}

loo <- rbindlist(all_iterations, use.names = TRUE, fill = TRUE)
loo_qc <- rbindlist(qc_iterations, use.names = TRUE, fill = TRUE)
comparison <- merge(
  loo, full, by = c("cell_type", "contrast", "set_id", "set_label"), all.x = TRUE
)
comparison[, direction_concordant := sign(NES) == sign(full_NES)]
comparison[, globally_significant_same_direction :=
             !is.na(LOO_global_primary_FDR) & LOO_global_primary_FDR < 0.05 &
             direction_concordant]

stability <- comparison[, .(
  full_NES = first(full_NES), full_pval = first(full_pval),
  full_comparison_FDR = first(full_comparison_FDR),
  full_global_primary_FDR = first(full_global_primary_FDR),
  valid_iterations = sum(is.finite(NES)), total_donors = .N,
  direction_concordant_n = sum(direction_concordant, na.rm = TRUE),
  direction_concordant_fraction = mean(direction_concordant, na.rm = TRUE),
  sign_flip_n = sum(!direction_concordant, na.rm = TRUE),
  median_LOO_NES = median(NES, na.rm = TRUE),
  minimum_LOO_NES = min(NES, na.rm = TRUE),
  maximum_LOO_NES = max(NES, na.rm = TRUE),
  global_FDR_lt_0_05_n = sum(LOO_global_primary_FDR < 0.05, na.rm = TRUE),
  global_FDR_lt_0_05_same_direction_n = sum(globally_significant_same_direction, na.rm = TRUE)
), by = .(cell_type, contrast, set_id, set_label)]
stability[, robustness_label := fcase(
  direction_concordant_fraction == 1 & full_global_primary_FDR < 0.05,
  "full_global_significant_and_100pct_direction_stable",
  direction_concordant_fraction >= 0.9 & full_global_primary_FDR < 0.05,
  "full_global_significant_and_ge90pct_direction_stable",
  direction_concordant_fraction == 1,
  "100pct_direction_stable_full_not_global_significant",
  direction_concordant_fraction >= 0.9,
  "ge90pct_direction_stable_full_not_global_significant",
  default = "direction_sensitive"
)]
stability[, inference_note := paste(
  "LOO direction stability is the primary robustness metric;",
  "LOO significance retention is descriptive and not an independent replication"
)]

setorder(comparison, omitted_donor, cell_type, contrast, set_id)
setorder(stability, full_global_primary_FDR, cell_type, contrast, set_id)
fwrite(comparison, file.path(out_dir, "fixed_7_family_LOO_all_iterations.tsv"),
       sep = "\t", quote = FALSE, na = "")
fwrite(stability, file.path(out_dir, "fixed_7_family_LOO_stability_summary.tsv"),
       sep = "\t", quote = FALSE, na = "")
fwrite(loo_qc, file.path(out_dir, "fixed_7_family_LOO_iteration_QC.tsv"),
       sep = "\t", quote = FALSE, na = "")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "program_LOO_sessionInfo.txt"))

record_status(
  "program_LOO", "complete",
  sprintf(
    "donors=%d;axes=%d;complete_iterations=%d;full-significant axes with >=90pct direction stability=%d",
    length(donor_universe), nrow(stability), sum(loo_qc$complete_iteration),
    nrow(stability[full_global_primary_FDR < 0.05 & direction_concordant_fraction >= 0.9])
  ), 100
)

print(loo_qc)
print(stability[full_global_primary_FDR < 0.05])
