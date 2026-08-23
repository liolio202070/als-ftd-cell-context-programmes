#!/usr/bin/env Rscript

# Donor-blocked C9-ALS/FTD pseudobulk analysis for the Gittings 2023 cohort.
# Frontal cortex is primary; occipital cortex is a pre-defined regional control.
# The ALS-vs-FTD contrast is a difference of within-partition disease effects,
# because ALS and FTD were distributed across separate CELLxGENE assets.

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
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(fgsea)
})

set.seed(20260811)

out_dir <- file.path(project_root, "results", "Gittings_2023_C9_snRNA")
counts_file <- file.path(out_dir, "pseudobulk_raw_counts.h5")
metadata_file <- file.path(out_dir, "pseudobulk_sample_metadata.tsv")
status_file <- file.path(out_dir, "status.tsv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

candidate_genes <- c(
  "BAG3", "HSP90AA1", "HSP90AB1", "CHORDC1", "ARHGAP35", "PARK7",
  "CCT5", "ARHGAP39", "COX5B", "COX7C", "MAP1B", "MAPT", "BCR",
  "CCT8", "COA3", "DOCK7", "STARD13"
)

family_definition <- data.table(
  go_id = c(
    "GO:0006119", "GO:0042775", "GO:0006457", "GO:0009408", "GO:0050821",
    "GO:0051056", "GO:0007264", "GO:0007416", "GO:0034329", "GO:0007409",
    "GO:0061564"
  ),
  family_id = c(
    "mitochondrial_energy", "mitochondrial_energy", "protein_folding",
    "heat_response", "protein_stabilization", "small_GTPase", "small_GTPase",
    "synapse_junction", "synapse_junction", "axon_development", "axon_development"
  ),
  family_label = c(
    "Mitochondrial energy", "Mitochondrial energy", "Protein folding",
    "Heat response", "Protein stabilization", "Small GTPase",
    "Small GTPase", "Synapse/junction", "Synapse/junction",
    "Axon development", "Axon development"
  )
)

go_map <- as.data.table(AnnotationDbi::select(
  org.Hs.eg.db, keys = unique(family_definition$go_id), keytype = "GOALL",
  columns = c("GOALL", "SYMBOL", "ONTOLOGYALL")
))
go_map <- go_map[ONTOLOGYALL == "BP" & !is.na(SYMBOL) & nzchar(SYMBOL)]
membership <- unique(merge(
  family_definition, go_map[, .(go_id = GOALL, gene = SYMBOL)],
  by = "go_id", allow.cartesian = TRUE
))
membership <- membership[, .(go_id, gene, family_id, family_label)]
setorder(membership, family_id, go_id, gene)
fwrite(family_definition, file.path(out_dir, "predefined_7_family_GO_definition.tsv"),
       sep = "\t", quote = FALSE)
fwrite(membership, file.path(out_dir, "predefined_7_family_gene_membership.tsv"),
       sep = "\t", quote = FALSE)
pathways <- split(membership$gene, membership$family_id)
pathways <- lapply(pathways, unique)
family_labels <- unique(membership[, .(set_id = family_id, set_label = family_label)])
if (length(pathways) != 7L) stop("Expected exactly seven pre-defined mechanism families")

contrast_definitions <- list(
  C9_ALS_vs_Ctrl_Frontal = c(
    ALSPart_ALS_Frontal = 1, ALSPart_Ctrl_Frontal = -1
  ),
  C9_FTD_vs_Ctrl_Frontal = c(
    FTDPart_FTD_Frontal = 1, FTDPart_Ctrl_Frontal = -1
  ),
  C9_ALS_vs_C9_FTD_Frontal_batch_bridged = c(
    ALSPart_ALS_Frontal = 1, ALSPart_Ctrl_Frontal = -1,
    FTDPart_FTD_Frontal = -1, FTDPart_Ctrl_Frontal = 1
  ),
  C9_ALS_FTD_vs_Ctrl_Frontal_transition = c(
    ALSPart_ALSFTD_Frontal = 1, ALSPart_Ctrl_Frontal = -1
  ),
  C9_ALS_vs_Ctrl_Occipital = c(
    ALSPart_ALS_Occipital = 1, ALSPart_Ctrl_Occipital = -1
  ),
  C9_FTD_vs_Ctrl_Occipital = c(
    FTDPart_FTD_Occipital = 1, FTDPart_Ctrl_Occipital = -1
  ),
  C9_ALS_vs_C9_FTD_Occipital_batch_bridged = c(
    ALSPart_ALS_Occipital = 1, ALSPart_Ctrl_Occipital = -1,
    FTDPart_FTD_Occipital = -1, FTDPart_Ctrl_Occipital = 1
  ),
  C9_ALS_FTD_vs_Ctrl_Occipital_transition = c(
    ALSPart_ALSFTD_Occipital = 1, ALSPart_Ctrl_Occipital = -1
  ),
  C9_ALS_Frontal_minus_Occipital_interaction = c(
    ALSPart_ALS_Frontal = 1, ALSPart_Ctrl_Frontal = -1,
    ALSPart_ALS_Occipital = -1, ALSPart_Ctrl_Occipital = 1
  ),
  C9_FTD_Frontal_minus_Occipital_interaction = c(
    FTDPart_FTD_Frontal = 1, FTDPart_Ctrl_Frontal = -1,
    FTDPart_FTD_Occipital = -1, FTDPart_Ctrl_Occipital = 1
  ),
  C9_ALS_vs_C9_FTD_Frontal_minus_Occipital_interaction = c(
    ALSPart_ALS_Frontal = 1, ALSPart_Ctrl_Frontal = -1,
    FTDPart_FTD_Frontal = -1, FTDPart_Ctrl_Frontal = 1,
    ALSPart_ALS_Occipital = -1, ALSPart_Ctrl_Occipital = 1,
    FTDPart_FTD_Occipital = 1, FTDPart_Ctrl_Occipital = -1
  ),
  C9_ALS_FTD_Frontal_minus_Occipital_transition = c(
    ALSPart_ALSFTD_Frontal = 1, ALSPart_Ctrl_Frontal = -1,
    ALSPart_ALSFTD_Occipital = -1, ALSPart_Ctrl_Occipital = 1
  )
)

contrast_role <- function(name) {
  fcase(
    name %chin% c(
      "C9_ALS_vs_Ctrl_Frontal", "C9_FTD_vs_Ctrl_Frontal",
      "C9_ALS_vs_C9_FTD_Frontal_batch_bridged"
    ), "primary_frontal",
    grepl("transition", name), "transition_background",
    grepl("Occipital", name) & !grepl("interaction", name), "regional_control",
    grepl("interaction", name), "regional_interaction",
    default = "other_supportive"
  )
}

make_contrast_matrix <- function(definitions, design_columns) {
  eligible <- vapply(definitions, function(weights) {
    all(names(weights) %chin% design_columns)
  }, logical(1))
  definitions <- definitions[eligible]
  if (!length(definitions)) return(matrix(numeric(0), nrow = length(design_columns)))
  output <- matrix(0, nrow = length(design_columns), ncol = length(definitions),
                   dimnames = list(design_columns, names(definitions)))
  for (name in names(definitions)) {
    weights <- definitions[[name]]
    output[names(weights), name] <- weights
  }
  output
}

run_fgsea <- function(ranks) {
  local_pathways <- lapply(pathways, function(genes) intersect(genes, names(ranks)))
  local_pathways <- local_pathways[lengths(local_pathways) >= 10L]
  result <- as.data.table(suppressWarnings(fgseaMultilevel(
    pathways = local_pathways, stats = ranks, minSize = 10L, maxSize = 1000L,
    eps = 1e-10, scoreType = "std"
  )))
  result[, leadingEdge := vapply(leadingEdge, paste, collapse = ";",
                                 FUN.VALUE = character(1))]
  setnames(result, "pathway", "set_id")
  merge(result, family_labels, by = "set_id", all.x = TRUE)
}

record_status("program_model", "started",
              "Read donor-level raw-count pseudobulk; min30 nuclei", 77)
metadata <- fread(metadata_file)
counts <- h5read(counts_file, "counts")
sample_ids <- h5read(counts_file, "sample_id")
gene_ids <- h5read(counts_file, "gene_id")
gene_symbols <- h5read(counts_file, "gene_symbol")
sample_ids <- enc2utf8(as.character(sample_ids))
gene_ids <- enc2utf8(as.character(gene_ids))
gene_symbols <- enc2utf8(as.character(gene_symbols))
if (nrow(counts) == length(sample_ids) && ncol(counts) == length(gene_ids)) {
  # rhdf5 reverses the dimension order of datasets written by h5py.
  counts <- t(counts)
}
stopifnot(nrow(counts) == length(gene_ids), ncol(counts) == length(sample_ids))
metadata <- metadata[match(sample_ids, sample_id)]
if (anyNA(metadata$sample_id) || !identical(metadata$sample_id, sample_ids)) {
  stop("Pseudobulk metadata/sample order mismatch")
}
rownames(counts) <- make.unique(gene_ids)
colnames(counts) <- sample_ids
metadata[, region_model := fifelse(tolower(region) == "frontal", "Frontal", "Occipital")]
metadata[, partition_model := fifelse(
  grepl("ALS_ALSFTD", partition), "ALSPart", "FTDPart"
)]
metadata[, diagnosis_model := fcase(
  diagnosis == "Control", "Ctrl", diagnosis == "C9_ALS", "ALS",
  diagnosis == "C9_ALS_FTD", "ALSFTD", diagnosis == "C9_FTD", "FTD",
  default = NA_character_
)]
metadata[, group := paste(partition_model, diagnosis_model, region_model, sep = "_")]
metadata[, group := gsub("_", "", group, fixed = TRUE)]
# Restore readable separators used by the pre-defined contrast vectors.
metadata[, group := sub("^ALSPart", "ALSPart_", group)]
metadata[, group := sub("^FTDPart", "FTDPart_", group)]
metadata[, group := sub("(Ctrl|ALSFTD|ALS|FTD)(Frontal|Occipital)$", "\\1_\\2", group)]
metadata[, sex := factor(sex)]
metadata[, age_z := as.numeric(scale(age))]
if (anyNA(metadata$diagnosis_model) || anyNA(metadata$age_z)) {
  stop("Unmapped diagnosis or missing/non-variable age in pseudobulk metadata")
}

metadata_primary <- metadata[n_nuclei >= 30]
sample_audit <- metadata_primary[, .(
  samples = .N, donors = uniqueN(donor_id),
  median_nuclei = as.numeric(median(n_nuclei)),
  minimum_nuclei = as.numeric(min(n_nuclei)),
  maximum_nuclei = as.numeric(max(n_nuclei))
), by = .(cell_type, partition_model, diagnosis, region_model)]
setorder(sample_audit, cell_type, partition_model, diagnosis, region_model)
fwrite(sample_audit, file.path(out_dir, "min30_sample_donor_audit.tsv"),
       sep = "\t", quote = FALSE)

fit_cell_type <- function(cell_type_value) {
  meta <- copy(metadata_primary[cell_type == cell_type_value])
  meta <- meta[match(colnames(counts), sample_id, nomatch = 0L)]
  if (nrow(meta) < 12L || uniqueN(meta$donor_id) < 8L) {
    return(list(skip = data.table(cell_type = cell_type_value,
                                  reason = "insufficient_samples_or_donors")))
  }
  count_subset <- counts[, meta$sample_id, drop = FALSE]
  meta[, group := factor(group)]
  design <- model.matrix(~ 0 + group + sex + age_z, data = meta)
  colnames(design) <- sub("^group", "", colnames(design))
  if (qr(design)$rank != ncol(design)) {
    return(list(skip = data.table(
      cell_type = cell_type_value, reason = "rank_deficient_sex_age_design",
      samples = nrow(meta), donors = uniqueN(meta$donor_id),
      design_columns = ncol(design), design_rank = qr(design)$rank
    )))
  }
  group_donors <- meta[, .(donors = uniqueN(donor_id)), by = group]
  donor_count <- setNames(group_donors$donors, as.character(group_donors$group))
  eligible_definitions <- contrast_definitions[vapply(
    contrast_definitions,
    function(weights) all(names(weights) %chin% names(donor_count)) &&
      all(donor_count[names(weights)] >= 3L),
    logical(1)
  )]
  contrast_matrix <- make_contrast_matrix(eligible_definitions, colnames(design))
  if (!ncol(contrast_matrix)) {
    return(list(skip = data.table(cell_type = cell_type_value,
                                  reason = "no_estimable_predefined_contrast")))
  }

  y <- DGEList(counts = count_subset)
  keep <- filterByExpr(y, design = design, min.count = 5, min.total.count = 15)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y, method = "TMM")
  v0 <- voom(y, design, plot = FALSE)
  correlation_fit <- duplicateCorrelation(v0, design, block = meta$donor_id)
  consensus <- correlation_fit$consensus.correlation
  if (!is.finite(consensus)) consensus <- 0
  v <- voom(y, design, plot = FALSE, block = meta$donor_id,
            correlation = consensus)
  fit <- lmFit(v, design, block = meta$donor_id, correlation = consensus)
  fit <- eBayes(contrasts.fit(fit, contrast_matrix), robust = TRUE)

  retained_gene_id <- rownames(fit$coefficients)
  symbol_lookup <- gene_symbols[match(retained_gene_id, gene_ids)]
  genes_out <- list()
  gsea_out <- list()
  for (contrast_name in colnames(contrast_matrix)) {
    table <- as.data.table(topTable(
      fit, coef = contrast_name, number = Inf, sort.by = "none", adjust.method = "BH"
    ), keep.rownames = "gene_id")
    table[, gene := symbol_lookup]
    setnames(table, c("logFC", "P.Value", "adj.P.Val"),
             c("log2FC", "PValue", "FDR"))
    table[, `:=`(
      cell_type = cell_type_value, contrast = contrast_name,
      analysis_role = contrast_role(contrast_name), n_samples = nrow(meta),
      n_donors = uniqueN(meta$donor_id), duplicate_correlation = consensus
    )]
    genes_out[[contrast_name]] <- table[gene %chin% candidate_genes]

    rank_table <- table[!is.na(gene) & nzchar(gene) & is.finite(t)]
    rank_table[, abs_t := abs(t)]
    setorder(rank_table, gene, -AveExpr, -abs_t)
    rank_table <- rank_table[, .SD[1L], by = gene]
    setorder(rank_table, -t, gene)
    ranks <- rank_table$t
    names(ranks) <- rank_table$gene
    fg <- run_fgsea(ranks)
    fg[, `:=`(
      cell_type = cell_type_value, contrast = contrast_name,
      analysis_role = contrast_role(contrast_name), ranked_genes = length(ranks),
      n_samples = nrow(meta), n_donors = uniqueN(meta$donor_id),
      duplicate_correlation = consensus
    )]
    fg[, comparison_FDR := p.adjust(pval, method = "BH")]
    gsea_out[[contrast_name]] <- fg
  }
  list(
    genes = rbindlist(genes_out, use.names = TRUE, fill = TRUE),
    gsea = rbindlist(gsea_out, use.names = TRUE, fill = TRUE),
    qc = data.table(
      cell_type = cell_type_value, input_samples = nrow(meta),
      donors = uniqueN(meta$donor_id), input_genes = nrow(count_subset),
      retained_genes = nrow(y), design_columns = ncol(design),
      design_rank = qr(design)$rank,
      residual_df_nominal = nrow(design) - qr(design)$rank,
      duplicate_correlation = consensus,
      tested_contrasts = ncol(contrast_matrix),
      tested_primary_contrasts = sum(
        vapply(colnames(contrast_matrix), contrast_role, character(1)) == "primary_frontal"
      ),
      minimum_donors_per_weighted_group = 3L,
      covariates = "group_cell_means+sex+age_z;donor_block"
    ),
    skip = data.table()
  )
}

cell_types <- sort(unique(metadata_primary$cell_type))
fits <- vector("list", length(cell_types))
for (index in seq_along(cell_types)) {
  fits[[index]] <- fit_cell_type(cell_types[index])
  record_status(
    "program_model", "running",
    sprintf("fit %d/%d cell_type=%s", index, length(cell_types), cell_types[index]),
    77 + 10 * index / length(cell_types)
  )
}

gene_results <- rbindlist(lapply(fits, function(x) x$genes), use.names = TRUE, fill = TRUE)
gsea_results <- rbindlist(lapply(fits, function(x) x$gsea), use.names = TRUE, fill = TRUE)
model_qc <- rbindlist(lapply(fits, function(x) x$qc), use.names = TRUE, fill = TRUE)
skipped <- rbindlist(lapply(fits, function(x) x$skip), use.names = TRUE, fill = TRUE)
if (!nrow(gsea_results)) stop("No eligible program tests")
if (anyDuplicated(gsea_results[, .(cell_type, contrast, set_id)])) {
  stop("Duplicate cell-type/contrast/family rows")
}

gsea_results[, global_primary_FDR := NA_real_]
gsea_results[analysis_role == "primary_frontal",
             global_primary_FDR := p.adjust(pval, method = "BH")]
gsea_results[, supportive_global_FDR := NA_real_]
gsea_results[analysis_role %chin% c("regional_control", "regional_interaction"),
             supportive_global_FDR := p.adjust(pval, method = "BH")]
gsea_results[, transition_global_FDR := NA_real_]
gsea_results[analysis_role == "transition_background",
             transition_global_FDR := p.adjust(pval, method = "BH")]

candidate_grid <- CJ(
  cell_type = unique(model_qc$cell_type), contrast = unique(gsea_results$contrast),
  gene = candidate_genes, unique = TRUE
)
candidate_complete <- merge(
  candidate_grid, gene_results, by = c("cell_type", "contrast", "gene"), all.x = TRUE
)
candidate_complete[, tested := !is.na(PValue)]
candidate_complete[is.na(analysis_role), analysis_role := contrast_role(contrast)]
candidate_complete[, non_test_reason := fifelse(
  tested, "", "not_retained_after_predefined_expression_filter"
)]

expected_primary <- 7L * sum(model_qc$tested_primary_contrasts)
observed_primary <- nrow(gsea_results[analysis_role == "primary_frontal"])
test_qc <- data.table(
  primary_cell_types = uniqueN(gsea_results[analysis_role == "primary_frontal", cell_type]),
  expected_primary_family_tests = expected_primary,
  observed_primary_family_tests = observed_primary,
  primary_unique = !anyDuplicated(gsea_results[
    analysis_role == "primary_frontal", .(cell_type, contrast, set_id)
  ]),
  primary_FDR_scope = paste(
    "BH across estimable frontal contrasts x 7 fixed families x eligible cell types;",
    "each weighted diagnosis/partition/region group requires >=3 donors"
  ),
  supportive_FDR_scope = "BH across regional-control and region-interaction tests",
  transition_FDR_scope = "BH only across ALS-FTD transition-background tests",
  ALS_vs_FTD_interpretation = paste(
    "Difference of within-partition disease effects, anchored by partition-specific controls;",
    "not an unadjusted direct phenotype comparison and not independent of asset batch"
  )
)
if (observed_primary != expected_primary) {
  stop(sprintf("Primary test count mismatch: %d != %d", observed_primary, expected_primary))
}

setorder(gsea_results, analysis_role, cell_type, contrast, pval)
setorder(candidate_complete, cell_type, contrast, PValue, gene, na.last = TRUE)
fwrite(gsea_results, file.path(out_dir, "fixed_7_family_fgsea_results.tsv"),
       sep = "\t", quote = FALSE, na = "")
fwrite(candidate_complete, file.path(out_dir, "predefined_17_gene_results.tsv"),
       sep = "\t", quote = FALSE, na = "")
fwrite(model_qc, file.path(out_dir, "program_model_QC.tsv"),
       sep = "\t", quote = FALSE, na = "")
fwrite(skipped, file.path(out_dir, "program_model_skipped_celltypes.tsv"),
       sep = "\t", quote = FALSE, na = "")
fwrite(test_qc, file.path(out_dir, "program_test_FDR_QC.tsv"),
       sep = "\t", quote = FALSE, na = "")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "program_model_sessionInfo.txt"))

record_status(
  "program_model", "complete",
  sprintf(
    "cell_types=%d;primary_tests=%d;global_FDR_sig=%d;batch-bridged ALS-vs-FTD labelled",
    nrow(model_qc), observed_primary,
    sum(gsea_results$global_primary_FDR < 0.05, na.rm = TRUE)
  ), 88
)

print(model_qc)
print(gsea_results[analysis_role == "primary_frontal", .(
  cell_type, contrast, set_label, NES, pval, comparison_FDR, global_primary_FDR
)][order(global_primary_FDR, pval)])
