#!/usr/bin/env Rscript

# Independent FTD-GRN/TDP-43 human cortex validation using GEO GSE163122.
# The processed matrices are FANS-sorted NEUN+ neuronal and OLIG2+ lineage
# nuclei from frontal, temporal and occipital cortex (4 controls, 5 FTD-GRN).
# This validates cross-FTD-subtype generalizability, not C9 specificity.

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
})

set.seed(20260811)

raw_dir <- file.path(project_root, "data", "raw", "GSE163122_FTD_GRN")
out_dir <- file.path(project_root, "results", "GSE163122_FTD_GRN_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
status_file <- file.path(out_dir, "status.tsv")
if (file.exists(status_file)) file.remove(status_file)

record_status <- function(stage, state, detail, progress) {
  row <- data.table(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stage = stage, state = state, detail = detail, progress = progress
  )
  fwrite(row, status_file, sep = "\t", quote = FALSE,
         append = file.exists(status_file), col.names = !file.exists(status_file))
  cat(sprintf("[%s] %s | %s | %s | %.1f%%\n", row$timestamp, stage, state, detail, progress))
  flush.console()
}

xml_unescape <- function(x) {
  x <- gsub("&gt;", ">", x, fixed = TRUE)
  x <- gsub("&lt;", "<", x, fixed = TRUE)
  x <- gsub("&amp;", "&", x, fixed = TRUE)
  trimws(x)
}

tag_value <- function(lines, tag) {
  hit <- grep(sprintf("<%s(?: [^>]*)?>", tag), lines, perl = TRUE)
  if (!length(hit)) return(NA_character_)
  index <- hit[1L]
  same_line <- sub(sprintf("^.*<%s(?: [^>]*)?>(.*?)</%s>.*$", tag, tag), "\\1", lines[index], perl = TRUE)
  if (!identical(same_line, lines[index])) return(xml_unescape(same_line))
  if (index < length(lines)) return(xml_unescape(lines[index + 1L]))
  NA_character_
}

parse_miniml <- function(path) {
  lines <- readLines(path, warn = FALSE)
  starts <- grep('^  <Sample iid="', lines)
  stops <- grep('^  </Sample>', lines)
  stopifnot(length(starts) == length(stops))
  rows <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    block <- lines[starts[i]:stops[i]]
    fields <- list()
    char_hits <- grep('^      <Characteristics tag="', block)
    if (length(char_hits)) {
      for (index in char_hits) {
        field <- sub('^.*tag="([^"]+)".*$', "\\1", block[index])
        fields[[field]] <- if (index < length(block)) xml_unescape(block[index + 1L]) else NA_character_
      }
    }
    title <- tag_value(block, "Title")
    rows[[i]] <- data.table(
      gsm = sub('^.*iid="([^"]+)".*$', "\\1", block[1L]),
      sample_title = title,
      donor_field = fields[["donor"]] %||% NA_character_,
      sex = fields[["Sex"]] %||% NA_character_,
      age = suppressWarnings(as.numeric(fields[["age"]] %||% NA_character_)),
      brain_region = fields[["brain region"]] %||% NA_character_,
      nuclei_population = fields[["nuclei population"]] %||% NA_character_,
      rin = suppressWarnings(as.numeric(fields[["rin tissue block"]] %||% NA_character_))
    )
  }
  rbindlist(rows, fill = TRUE)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

record_status("metadata", "started", "Parse GEO MINiML and audit donor labels", 3)
metadata <- parse_miniml(file.path(raw_dir, "GSE163122_family.xml"))
metadata <- metadata[grepl("_(neun|olig)$", sample_title)]
metadata[, donor_title := sub("^([CP][0-9]+).*$", "\\1", sample_title)]
metadata[, donor_label_discordant := donor_field != donor_title]
# GEO has a known internal donor-label inconsistency for P4D; the sample title,
# age/sex and paired population identify it as donor P4. Use the title-derived
# donor consistently and retain the discordance flag in QC.
metadata[, donor_id := donor_title]
metadata[, diagnosis := fifelse(grepl("^P", donor_id), "FTD_GRN", "Control")]
metadata[, region := fcase(
  brain_region == "Frontal Cortex", "Frontal",
  brain_region == "Temporal Cortex", "Temporal",
  brain_region == "Occipital Cortex", "Occipital",
  default = NA_character_
)]
metadata[, population := fifelse(grepl("_olig$", sample_title), "OLIG2pos", "NEUNpos")]
metadata[, condition_region := paste(fifelse(diagnosis == "FTD_GRN", "FTD", "Ctrl"), region, sep = "_")]

metadata_qc <- metadata[, .(
  samples = .N,
  donors = uniqueN(donor_id),
  controls = uniqueN(donor_id[diagnosis == "Control"]),
  FTD_GRN = uniqueN(donor_id[diagnosis == "FTD_GRN"]),
  regions = uniqueN(region),
  missing_RIN = sum(is.na(rin)),
  donor_label_discordances = sum(donor_label_discordant)
), by = population]
fwrite(metadata, file.path(out_dir, "GSE163122_sample_metadata.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(metadata_qc, file.path(out_dir, "GSE163122_metadata_QC.tsv"), sep = "\t", quote = FALSE, na = "")

family_membership <- unique(fread(file.path(
  project_root, "results", "GSE219280_same_celltype_fixed_programs",
  "predefined_gene_set_membership.tsv"
))[, .(family_id, family_label, gene)])
family_pathways <- split(family_membership$gene, family_membership$family_id)
family_labels <- unique(family_membership[, .(set_id = family_id, set_label = family_label)])

required_levels <- c(
  "Ctrl_Frontal", "Ctrl_Temporal", "Ctrl_Occipital",
  "FTD_Frontal", "FTD_Temporal", "FTD_Occipital"
)
contrast_definitions <- c(
  FTD_vs_Ctrl_mean = "(FTD_Frontal + FTD_Temporal + FTD_Occipital - Ctrl_Frontal - Ctrl_Temporal - Ctrl_Occipital) / 3",
  FTD_vs_Ctrl_Frontal = "FTD_Frontal - Ctrl_Frontal",
  FTD_vs_Ctrl_Temporal = "FTD_Temporal - Ctrl_Temporal",
  FTD_vs_Ctrl_Occipital = "FTD_Occipital - Ctrl_Occipital",
  Frontal_interaction_vs_Occipital = "(FTD_Frontal - Ctrl_Frontal) - (FTD_Occipital - Ctrl_Occipital)",
  Temporal_interaction_vs_Occipital = "(FTD_Temporal - Ctrl_Temporal) - (FTD_Occipital - Ctrl_Occipital)"
)

run_fixed_fgsea <- function(ranks) {
  pathways <- lapply(family_pathways, function(genes) intersect(genes, names(ranks)))
  pathways <- pathways[lengths(pathways) >= 10L]
  result <- as.data.table(suppressWarnings(fgseaMultilevel(
    pathways = pathways, stats = ranks, minSize = 10L, maxSize = 1000L,
    eps = 1e-10, scoreType = "std"
  )))
  result[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
  setnames(result, "pathway", "set_id")
  merge(result, family_labels, by = "set_id", all.x = TRUE)
}

fit_population <- function(population_value) {
  matrix_file <- file.path(raw_dir, sprintf(
    "GSE163122_countmtx_%s.csv.gz",
    ifelse(population_value == "OLIG2pos", "olig", "neun")
  ))
  # Use base gzfile support to avoid an unnecessary R.utils dependency.
  input <- as.data.table(read.csv(gzfile(matrix_file), check.names = FALSE))
  gene_id <- input[[1L]]
  counts <- as.matrix(input[, -1L])
  storage.mode(counts) <- "integer"
  rownames(counts) <- gene_id
  rm(input)

  meta <- copy(metadata[population == population_value])
  meta <- meta[match(colnames(counts), sample_title)]
  if (anyNA(meta$sample_title)) stop("Count matrix column missing from metadata: ", population_value)
  if (!all(required_levels %chin% meta$condition_region)) stop("Required group missing: ", population_value)

  meta[, condition_region := factor(condition_region, levels = required_levels)]
  meta[, sex := factor(sex)]
  meta[, age_z := as.numeric(scale(age))]
  meta[, rin_imputed := fifelse(is.na(rin), median(rin, na.rm = TRUE), rin)]
  meta[, rin_z := as.numeric(scale(rin_imputed))]
  design <- model.matrix(~ 0 + condition_region + sex + age_z + rin_z, data = meta)
  colnames(design) <- sub("^condition_region", "", colnames(design))
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient design: ", population_value)

  y <- DGEList(counts = counts)
  keep <- filterByExpr(y, group = meta$condition_region, min.count = 5, min.total.count = 15)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y, method = "TMM")
  v0 <- voom(y, design, plot = FALSE)
  corfit <- duplicateCorrelation(v0, design, block = meta$donor_id)
  consensus <- corfit$consensus.correlation
  if (!is.finite(consensus)) consensus <- 0
  v <- voom(y, design, plot = FALSE, block = meta$donor_id, correlation = consensus)
  fit <- lmFit(v, design, block = meta$donor_id, correlation = consensus)
  contrast_matrix <- makeContrasts(contrasts = unname(contrast_definitions), levels = design)
  colnames(contrast_matrix) <- names(contrast_definitions)
  fit <- eBayes(contrasts.fit(fit, contrast_matrix), robust = TRUE)

  ensembl_ids <- rownames(fit$coefficients)
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = ensembl_ids, keytype = "ENSEMBL",
    column = "SYMBOL", multiVals = "first"
  )

  gene_results <- list()
  gsea_results <- list()
  for (contrast_name in colnames(contrast_matrix)) {
    table <- as.data.table(topTable(
      fit, coef = contrast_name, number = Inf, sort.by = "none", adjust.method = "BH"
    ), keep.rownames = "gene_id")
    table[, gene := unname(symbols[gene_id])]
    table[, `:=`(
      population = population_value,
      contrast = contrast_name,
      n_samples = nrow(meta),
      n_donors = uniqueN(meta$donor_id),
      residual_df_nominal = nrow(design) - qr(design)$rank,
      duplicate_correlation = consensus
    )]
    setnames(table, c("logFC", "P.Value", "adj.P.Val"), c("log2FC", "PValue", "FDR"))
    gene_results[[contrast_name]] <- table

    ranks_table <- table[!is.na(gene) & nzchar(gene) & is.finite(t)]
    ranks_table[, abs_t_for_dedup := abs(t)]
    setorder(ranks_table, gene, -AveExpr, -abs_t_for_dedup)
    ranks_table <- ranks_table[, .SD[1L], by = gene]
    setorder(ranks_table, -t, gene)
    ranks <- ranks_table$t
    names(ranks) <- ranks_table$gene
    fg <- run_fixed_fgsea(ranks)
    fg[, `:=`(
      population = population_value,
      contrast = contrast_name,
      analysis_role = fcase(
        contrast_name == "FTD_vs_Ctrl_mean", "primary_external_FTD",
        contrast_name == "FTD_vs_Ctrl_Frontal", "predefined_frontal_support",
        default = "regional_context"
      ),
      ranked_genes = length(ranks),
      n_samples = nrow(meta),
      n_donors = uniqueN(meta$donor_id),
      duplicate_correlation = consensus
    )]
    fg[, comparison_FDR := p.adjust(pval, method = "BH")]
    gsea_results[[contrast_name]] <- fg
  }
  list(
    genes = rbindlist(gene_results, use.names = TRUE, fill = TRUE),
    gsea = rbindlist(gsea_results, use.names = TRUE, fill = TRUE),
    qc = data.table(
      population = population_value,
      input_genes = nrow(counts), retained_genes = nrow(y),
      mapped_symbols = sum(!is.na(symbols)), samples = nrow(meta),
      donors = uniqueN(meta$donor_id), controls = uniqueN(meta[diagnosis == "Control", donor_id]),
      FTD_GRN = uniqueN(meta[diagnosis == "FTD_GRN", donor_id]),
      design_columns = ncol(design), residual_df_nominal = nrow(design) - qr(design)$rank,
      duplicate_correlation = consensus
    )
  )
}

record_status("model", "started", "Fit donor-blocked NEUN+ and OLIG2+ models", 20)
fits <- lapply(c("NEUNpos", "OLIG2pos"), fit_population)
genes <- rbindlist(lapply(fits, `[[`, "genes"), use.names = TRUE, fill = TRUE)
gsea <- rbindlist(lapply(fits, `[[`, "gsea"), use.names = TRUE, fill = TRUE)
qc <- rbindlist(lapply(fits, `[[`, "qc"), use.names = TRUE, fill = TRUE)

gsea[, global_primary_FDR := NA_real_]
gsea[analysis_role == "primary_external_FTD", global_primary_FDR := p.adjust(pval, method = "BH")]

source_ftd_long <- fread(file.path(
  project_root, "results", "GSE219280_same_celltype_fixed_programs",
  "fixed_family_fgsea_results.tsv"
))[
  threshold == "min30_primary" & contrast == "FTD_vs_Ctrl_mean" & cell_type %chin% c("OPC", "Oligo"),
  .(set_id, set_label, cell_type, NES, global_primary_FDR)
]
source_ftd <- dcast(
  source_ftd_long, set_id + set_label ~ cell_type,
  value.var = c("NES", "global_primary_FDR")
)
source_ftd[, source_mean_NES := rowMeans(cbind(NES_OPC, NES_Oligo), na.rm = TRUE)]
source_ftd[, source_any_global_significant :=
             (!is.na(global_primary_FDR_OPC) & global_primary_FDR_OPC < 0.05) |
             (!is.na(global_primary_FDR_Oligo) & global_primary_FDR_Oligo < 0.05)]
external_olig <- gsea[
  population == "OLIG2pos" & contrast == "FTD_vs_Ctrl_mean",
  .(set_id, set_label, external_NES = NES, external_pval = pval,
    external_comparison_FDR = comparison_FDR, external_global_primary_FDR = global_primary_FDR)
]
concordance <- merge(source_ftd, external_olig, by = c("set_id", "set_label"), all = TRUE)
concordance[, direction_concordant_with_source_mean := sign(source_mean_NES) == sign(external_NES)]
concordance[, strict_source_axis_concordant :=
              (!is.na(global_primary_FDR_OPC) & global_primary_FDR_OPC < 0.05 & sign(NES_OPC) == sign(external_NES)) |
              (!is.na(global_primary_FDR_Oligo) & global_primary_FDR_Oligo < 0.05 & sign(NES_Oligo) == sign(external_NES))]
concordance[, evidence_label := fcase(
  external_global_primary_FDR < 0.05 & strict_source_axis_concordant,
  "strict_directional_cross_subtype_support",
  external_global_primary_FDR < 0.05 & direction_concordant_with_source_mean,
  "external_significant_direction_matches_source_mean_but_source_axis_not_global_significant",
  external_global_primary_FDR < 0.05,
  "external_significant_discordance",
  direction_concordant_with_source_mean,
  "direction_only_external_not_significant",
  default = "discordant_external_not_significant"
)]
concordance[, interpretation := fcase(
  evidence_label == "strict_directional_cross_subtype_support",
  "Directionally concordant cross-genetic-subtype support; not C9-specific replication",
  evidence_label == "external_significant_discordance",
  "Significant external signal is directionally discordant; supports context dependence, not replication",
  grepl("external_not_significant", evidence_label),
  "External result is not significant; do not label as replication",
  default = "External signal does not meet strict source-axis replication criteria"
)]

setorder(genes, population, contrast, PValue)
setorder(gsea, population, contrast, pval)
fwrite(genes, file.path(out_dir, "all_gene_contrast_results.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(gsea, file.path(out_dir, "fixed_family_fgsea_results.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(qc, file.path(out_dir, "model_QC.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(concordance, file.path(out_dir, "OLIG2pos_GSE219280_concordance.tsv"), sep = "\t", quote = FALSE, na = "")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

record_status(
  "pipeline", "complete",
  sprintf(
    "2 populations; %d/14 external primary program signals at global FDR < 0.05",
    sum(gsea$global_primary_FDR < 0.05, na.rm = TRUE)
  ),
  100
)

print(metadata_qc)
print(qc)
print(gsea[analysis_role == "primary_external_FTD", .(
  population, set_label, NES, pval, comparison_FDR, global_primary_FDR
)][order(global_primary_FDR, pval)])
print(concordance)
