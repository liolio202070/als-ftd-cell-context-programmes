#!/usr/bin/env Rscript

# Predefined, donor-level same-cell-type pathway analysis for the C9orf72
# ALS-FTD spectrum. The primary analysis uses the minimum-30-nuclei limma
# results; minimum-20 and minimum-50 are threshold sensitivity analyses.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
project_root <- analysis_project_root()
configure_project_library(project_root)
project_lib <- file.path(project_root, "cache", "R_library")
.libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
  library(fgsea)
})

set.seed(20260811)

out_dir <- file.path(project_root, "results", "GSE219280_same_celltype_fixed_programs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
status_file <- file.path(out_dir, "status.tsv")

record_status <- function(stage, state, detail, progress) {
  row <- data.table(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stage = stage,
    state = state,
    detail = detail,
    progress = progress
  )
  fwrite(
    row,
    status_file,
    sep = "\t",
    quote = FALSE,
    append = file.exists(status_file),
    col.names = !file.exists(status_file)
  )
  cat(sprintf("[%s] %s | %s | %s | %.1f%%\n", row$timestamp, stage, state, detail, progress))
  flush.console()
}

if (file.exists(status_file)) file.remove(status_file)
record_status("initialize", "started", "Load fixed GO definitions and donor-level DE results", 1)

# These 11 terms were already used in the project's leading-edge analysis.
# They are grouped a priori into seven mechanism families to reduce the primary
# testing burden. Term-level results remain secondary and fully traceable.
term_definitions <- data.table(
  go_id = c(
    "GO:0006119", "GO:0042775", "GO:0006457", "GO:0009408",
    "GO:0050821", "GO:0051056", "GO:0007264", "GO:0007416",
    "GO:0034329", "GO:0007409", "GO:0061564"
  ),
  term_label = c(
    "Oxidative phosphorylation",
    "Mitochondrial ATP synthesis coupled electron transport",
    "Protein folding",
    "Heat response",
    "Protein stabilization",
    "Small GTPase regulation",
    "Small GTPase signaling",
    "Synapse assembly",
    "Cell junction assembly",
    "Axonogenesis",
    "Axon development"
  ),
  family_id = c(
    "mitochondrial_energy", "mitochondrial_energy", "protein_folding",
    "heat_response", "protein_stabilization", "small_GTPase",
    "small_GTPase", "synapse_junction", "synapse_junction",
    "axon_development", "axon_development"
  ),
  family_label = c(
    "Mitochondrial energy", "Mitochondrial energy", "Protein folding",
    "Heat response", "Protein stabilization", "Small GTPase",
    "Small GTPase", "Synapse / junction assembly", "Synapse / junction assembly",
    "Axon development", "Axon development"
  )
)

go_membership <- as.data.table(AnnotationDbi::select(
  org.Hs.eg.db,
  keys = term_definitions$go_id,
  keytype = "GO",
  columns = c("SYMBOL", "ONTOLOGY")
))
go_membership <- unique(go_membership[
  ONTOLOGY == "BP" & !is.na(SYMBOL) & nzchar(SYMBOL),
  .(go_id = GO, gene = SYMBOL)
])
go_membership <- merge(go_membership, term_definitions, by = "go_id", all.x = TRUE)

term_pathways <- split(go_membership$gene, go_membership$go_id)
family_pathways <- split(go_membership$gene, go_membership$family_id)
term_pathways <- lapply(term_pathways, unique)
family_pathways <- lapply(family_pathways, unique)

gene_set_manifest <- rbindlist(list(
  term_definitions[, .(
    level = "GO_term",
    set_id = go_id,
    set_label = term_label,
    family_id,
    family_label,
    source_go_ids = go_id,
    annotated_genes = lengths(term_pathways[go_id])
  )],
  unique(term_definitions[, .(family_id, family_label)])[, .(
    level = "mechanism_family",
    set_id = family_id,
    set_label = family_label,
    family_id,
    family_label,
    source_go_ids = vapply(
      family_id,
      function(id) paste(term_definitions[family_id == id, go_id], collapse = ";"),
      character(1)
    ),
    annotated_genes = lengths(family_pathways[family_id])
  )]
), use.names = TRUE)
fwrite(gene_set_manifest, file.path(out_dir, "predefined_gene_set_manifest.tsv"), sep = "\t", quote = FALSE)
fwrite(go_membership, file.path(out_dir, "predefined_gene_set_membership.tsv"), sep = "\t", quote = FALSE)

thresholds <- data.table(
  threshold = c("min30_primary", "min20_sensitivity", "min50_sensitivity"),
  input_dir = c(
    "GSE219280_pseudobulk_DE",
    "GSE219280_pseudobulk_DE_min20",
    "GSE219280_pseudobulk_DE_min50"
  ),
  minimum_nuclei = c(30L, 20L, 50L),
  analysis_role = c("primary", "sensitivity", "sensitivity")
)
cell_types <- c("Astro", "OPC", "Oligo")
contrasts <- c("ALS_vs_Ctrl_mean", "FTD_vs_Ctrl_mean", "ALS_vs_FTD_mean")

make_ranks <- function(dt) {
  x <- copy(dt[
    !is.na(gene) & nzchar(gene) & is.finite(t) & is.finite(log2FC) & is.finite(PValue),
    .(gene, log2FC, t, PValue)
  ])
  x <- x[order(gene, -abs(t), PValue)][, .SD[1L], by = gene]
  setorder(x, -t, gene)
  x[, rank_stat := t + rev(seq_len(.N)) * .Machine$double.eps]
  ranks <- x$rank_stat
  names(ranks) <- x$gene
  ranks
}

run_fixed_fgsea <- function(pathways, ranks, min_size) {
  tested_pathways <- lapply(pathways, function(genes) intersect(genes, names(ranks)))
  tested_pathways <- tested_pathways[lengths(tested_pathways) >= min_size]
  if (!length(tested_pathways)) return(data.table())
  result <- suppressWarnings(fgseaMultilevel(
    pathways = tested_pathways,
    stats = ranks,
    minSize = min_size,
    maxSize = 1000,
    eps = 1e-10,
    scoreType = "std"
  ))
  result <- as.data.table(result)
  if (nrow(result)) {
    result[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
  }
  result
}

family_results <- list()
term_results <- list()
qc_results <- list()
result_index <- 0L
total_comparisons <- nrow(thresholds) * length(cell_types) * length(contrasts)
comparison_index <- 0L

for (threshold_index in seq_len(nrow(thresholds))) {
  threshold_spec <- thresholds[threshold_index]
  input_file <- file.path(
    project_root, "results", threshold_spec$input_dir,
    "all_gene_contrast_results.tsv"
  )
  record_status(
    paste0("load_", threshold_spec$threshold),
    "started",
    sprintf("Read %s", basename(dirname(input_file))),
    3 + 82 * (threshold_index - 1) * length(cell_types) * length(contrasts) / total_comparisons
  )
  de <- fread(input_file)
  de <- de[
    analysis == "all_disease_core_glia" &
      cell_type %chin% cell_types &
      contrast %chin% contrasts
  ]
  expected_keys <- CJ(cell_type = cell_types, contrast = contrasts)
  observed_keys <- unique(de[, .(cell_type, contrast)])
  if (nrow(observed_keys) != nrow(expected_keys)) {
    stop("Missing same-cell-type contrast(s) in ", input_file)
  }

  for (cell_type_value in cell_types) {
    for (contrast_value in contrasts) {
      comparison_index <- comparison_index + 1L
      comparison <- de[cell_type == cell_type_value & contrast == contrast_value]
      ranks <- make_ranks(comparison)
      family <- run_fixed_fgsea(family_pathways, ranks, min_size = 10L)
      term <- run_fixed_fgsea(term_pathways, ranks, min_size = 5L)

      common_columns <- list(
        threshold = threshold_spec$threshold,
        minimum_nuclei = threshold_spec$minimum_nuclei,
        analysis_role = threshold_spec$analysis_role,
        cell_type = cell_type_value,
        contrast = contrast_value,
        ranked_genes = length(ranks)
      )
      if (nrow(family)) {
        for (column_name in names(common_columns)) set(family, j = column_name, value = common_columns[[column_name]])
        setnames(family, "pathway", "set_id")
        family <- merge(
          family,
          unique(term_definitions[, .(set_id = family_id, set_label = family_label)]),
          by = "set_id",
          all.x = TRUE
        )
      }
      if (nrow(term)) {
        for (column_name in names(common_columns)) set(term, j = column_name, value = common_columns[[column_name]])
        setnames(term, "pathway", "set_id")
        term <- merge(
          term,
          term_definitions[, .(set_id = go_id, set_label = term_label, family_id, family_label)],
          by = "set_id",
          all.x = TRUE
        )
      }

      result_index <- result_index + 1L
      family_results[[result_index]] <- family
      term_results[[result_index]] <- term
      qc_results[[result_index]] <- data.table(
        threshold = threshold_spec$threshold,
        minimum_nuclei = threshold_spec$minimum_nuclei,
        analysis_role = threshold_spec$analysis_role,
        cell_type = cell_type_value,
        contrast = contrast_value,
        ranked_genes = length(ranks),
        family_sets_tested = nrow(family),
        term_sets_tested = nrow(term),
        duplicate_gene_rows = nrow(comparison) - uniqueN(comparison$gene)
      )

      progress <- 5 + 82 * comparison_index / total_comparisons
      record_status(
        paste0("gsea_", threshold_spec$threshold),
        "running",
        sprintf("%s | %s (%d/%d)", cell_type_value, contrast_value, comparison_index, total_comparisons),
        progress
      )
    }
  }
  rm(de)
  gc(verbose = FALSE)
}

family <- rbindlist(family_results, use.names = TRUE, fill = TRUE)
term <- rbindlist(term_results, use.names = TRUE, fill = TRUE)
qc <- rbindlist(qc_results, use.names = TRUE, fill = TRUE)

family[, comparison_FDR := p.adjust(pval, method = "BH"), by = .(threshold, cell_type, contrast)]
term[, comparison_FDR := p.adjust(pval, method = "BH"), by = .(threshold, cell_type, contrast)]
family[, global_primary_FDR := NA_real_]
term[, global_primary_FDR := NA_real_]
family[threshold == "min30_primary", global_primary_FDR := p.adjust(pval, method = "BH")]
term[threshold == "min30_primary", global_primary_FDR := p.adjust(pval, method = "BH")]
family[, direction := fifelse(NES > 0, "positive", "negative")]
term[, direction := fifelse(NES > 0, "positive", "negative")]

family_column_order <- c(
  "threshold", "minimum_nuclei", "analysis_role", "cell_type", "contrast",
  "set_id", "set_label", "NES", "pval", "padj", "comparison_FDR",
  "global_primary_FDR", "direction", "size", "ES", "nMoreExtreme",
  "ranked_genes", "leadingEdge"
)
term_column_order <- c(
  "threshold", "minimum_nuclei", "analysis_role", "cell_type", "contrast",
  "set_id", "set_label", "family_id", "family_label", "NES", "pval",
  "padj", "comparison_FDR", "global_primary_FDR", "direction", "size",
  "ES", "nMoreExtreme", "ranked_genes", "leadingEdge"
)
# fgseaMultilevel versions differ in whether nMoreExtreme is returned.
setcolorder(family, c(intersect(family_column_order, names(family)), setdiff(names(family), family_column_order)))
setcolorder(term, c(intersect(term_column_order, names(term)), setdiff(names(term), term_column_order)))

fwrite(family, file.path(out_dir, "fixed_family_fgsea_results.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(term, file.path(out_dir, "fixed_GO_term_fgsea_results.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(qc, file.path(out_dir, "comparison_QC.tsv"), sep = "\t", quote = FALSE, na = "")

record_status("summarize", "started", "Build phenotype and threshold-stability evidence tables", 90)

# One row per cell type and mechanism family. ALS-vs-FTD is a formal direct
# contrast; it is not inferred by comparing two separate disease-control P values.
wide <- dcast(
  family,
  threshold + minimum_nuclei + analysis_role + cell_type + set_id + set_label ~ contrast,
  value.var = c("NES", "pval", "comparison_FDR", "global_primary_FDR")
)
wide[, disease_control_pattern := fcase(
  sign(NES_ALS_vs_Ctrl_mean) == sign(NES_FTD_vs_Ctrl_mean), "shared_direction",
  sign(NES_ALS_vs_Ctrl_mean) != sign(NES_FTD_vs_Ctrl_mean), "opposite_direction",
  default = "undetermined"
)]
wide[, direct_phenotype_direction := fcase(
  NES_ALS_vs_FTD_mean > 0, "ALS_higher",
  NES_ALS_vs_FTD_mean < 0, "FTD_higher",
  default = "no_direction"
)]
wide[, direct_phenotype_FDR_0_05 := comparison_FDR_ALS_vs_FTD_mean < 0.05]
fwrite(wide, file.path(out_dir, "same_celltype_phenotype_evidence.tsv"), sep = "\t", quote = FALSE, na = "")

primary <- family[threshold == "min30_primary", .(
  cell_type, contrast, set_id, set_label,
  primary_NES = NES,
  primary_pval = pval,
  primary_comparison_FDR = comparison_FDR,
  primary_global_FDR = global_primary_FDR
)]
sensitivity <- merge(
  family[threshold != "min30_primary"],
  primary,
  by = c("cell_type", "contrast", "set_id", "set_label"),
  all.x = TRUE
)
sensitivity[, direction_retained := sign(NES) == sign(primary_NES)]
sensitivity[, absolute_NES_difference := abs(NES - primary_NES)]
fwrite(sensitivity, file.path(out_dir, "threshold_sensitivity_by_program.tsv"), sep = "\t", quote = FALSE, na = "")

sensitivity_summary <- sensitivity[, .(
  programs = .N,
  direction_agreement = mean(direction_retained),
  NES_pearson = cor(NES, primary_NES, method = "pearson"),
  NES_spearman = cor(NES, primary_NES, method = "spearman"),
  median_absolute_NES_difference = median(absolute_NES_difference),
  primary_comparison_FDR_0_05 = sum(primary_comparison_FDR < 0.05),
  sensitivity_comparison_FDR_0_05 = sum(comparison_FDR < 0.05)
), by = .(threshold, minimum_nuclei, cell_type, contrast)]
fwrite(sensitivity_summary, file.path(out_dir, "threshold_sensitivity_summary.tsv"), sep = "\t", quote = FALSE, na = "")

primary_summary <- family[threshold == "min30_primary", .(
  family_tests = .N,
  nominal_P_0_05 = sum(pval < 0.05),
  comparison_FDR_0_05 = sum(comparison_FDR < 0.05),
  global_primary_FDR_0_05 = sum(global_primary_FDR < 0.05)
), by = .(cell_type, contrast)]
fwrite(primary_summary, file.path(out_dir, "primary_test_summary.tsv"), sep = "\t", quote = FALSE, na = "")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
record_status(
  "pipeline",
  "complete",
  sprintf(
    "%d primary family tests; %d primary global-FDR signals",
    nrow(family[threshold == "min30_primary"]),
    nrow(family[threshold == "min30_primary" & global_primary_FDR < 0.05])
  ),
  100
)

print(primary_summary)
print(family[
  threshold == "min30_primary" & (comparison_FDR < 0.05 | global_primary_FDR < 0.05),
  .(cell_type, contrast, set_label, NES, pval, comparison_FDR, global_primary_FDR)
][order(global_primary_FDR, comparison_FDR, pval)])
