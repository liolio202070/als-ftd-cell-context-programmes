#!/usr/bin/env Rscript

# Frozen, donor-level external validation in Ruf et al. 2026 same-nucleus
# motor-cortex RNA+ATAC.  No data-driven redefinition of the seven programmes
# or the 17 candidate genes is permitted in this script.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
root <- analysis_project_root()
configure_project_library(root)
.libPaths(c(file.path(root, "cache/R_library"), .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(edgeR)
  library(limma)
  library(fgsea)
  library(ArchR)
})

set.seed(20260813)
out <- file.path(root, "results/Ruf_2026_frozen_validation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

message_log <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
              paste0(..., collapse = "")))
  flush.console()
}

safe_bh <- function(p) {
  ans <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  ans[ok] <- p.adjust(p[ok], method = "BH")
  ans
}

families <- fread(file.path(
  root, "results/Gittings_2023_C9_snRNA/predefined_7_family_gene_membership.tsv"
))
families <- unique(families[, .(family_id, family_label, gene)])
pathways <- split(families$gene, families$family_id)
pathways <- lapply(pathways, unique)
family_labels <- unique(families[, .(family_id, family_label)])
stopifnot(length(pathways) == 7L)

shortlist <- fread(file.path(
  root, "results/leading_edge_evidence_v2/leading_edge_driver_shortlist.tsv"
))
candidate_genes <- unique(shortlist$gene)
stopifnot(length(candidate_genes) == 17L)

meta <- fread(file.path(out, "pseudobulk_sample_metadata.tsv"))

load_counts <- function(modality, cell_type) {
  path <- file.path(out, sprintf("%s_pseudobulk_counts_%s.rds", modality, cell_type))
  value <- readRDS(path)
  cell_type_value <- cell_type
  expected <- meta[cell_type == cell_type_value, sample_id]
  stopifnot(identical(colnames(value), expected))
  value
}

meta[, diagnosis := factor(diagnosis, levels = c("HC", "ALS", "ALS_FTD"))]
meta[, case_type := factor(case_type,
                           levels = c("HC", "ALS", "ALS_FTD", "C9_ALS_FTD"))]
meta[, sex := factor(sex)]
meta[, log10_RNA_UMI := log10(pmax(median_RNA_UMI, 1))]
meta[, log10_ATAC_fragments := log10(pmax(median_ATAC_fragments, 1))]

cell_type_map <- c(
  Astrocytes = "Astrocyte", Exc_Neurons = "Excitatory neuron",
  Inh_Neurons = "Inhibitory neuron", Microglia = "Microglia",
  OPC = "OPC", Oligodendrocytes = "Oligodendrocyte"
)
meta[, cell_type_label := unname(cell_type_map[cell_type])]
stopifnot(!anyNA(meta$cell_type_label))

primary_contrasts <- list(
  ALS_vs_Control = c(diagnosisALS = 1, diagnosisHC = -1),
  ALSFTD_vs_Control = c(diagnosisALS_FTD = 1, diagnosisHC = -1),
  ALS_vs_ALSFTD = c(diagnosisALS = 1, diagnosisALS_FTD = -1)
)
c9_contrast <- list(
  C9_ALSFTD_vs_Control_sensitivity = c(case_typeC9_ALS_FTD = 1, case_typeHC = -1)
)

make_design <- function(m, group_var, technical_variable = NULL, age = FALSE) {
  m <- copy(m)
  m[, analysis_group := droplevels(get(group_var))]
  rhs <- c("0 + analysis_group")
  if (uniqueN(m$sex) > 1L) rhs <- c(rhs, "sex")
  if (!is.null(technical_variable)) {
    if (is.finite(sd(m[[technical_variable]], na.rm = TRUE)) &&
        sd(m[[technical_variable]], na.rm = TRUE) > 0)
      rhs <- c(rhs, sprintf("scale(%s)", technical_variable))
  }
  if (age && sum(is.finite(m$age)) == nrow(m) && sd(m$age) > 0)
    rhs <- c(rhs, "scale(age)")
  design <- model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), m)
  colnames(design) <- sub("^analysis_group", paste0(group_var), colnames(design))
  if (qr(design)$rank != ncol(design)) stop("Non-full-rank design")
  design
}

make_contrast <- function(weights, design_cols) {
  if (!all(names(weights) %chin% design_cols)) return(NULL)
  ans <- rep(0, length(design_cols)); names(ans) <- design_cols
  ans[names(weights)] <- weights
  matrix(ans, ncol = 1L, dimnames = list(design_cols, NULL))
}

fit_count_model <- function(counts, m, contrast_defs, group_var = "diagnosis",
                            technical_variable = NULL, age = FALSE,
                            filter = TRUE) {
  design <- make_design(m, group_var, technical_variable = technical_variable,
                        age = age)
  y <- DGEList(counts = counts)
  keep <- if (filter) filterByExpr(y, design = design) else rowSums(y$counts) >= 5
  if (sum(keep) < 20L) stop("Too few retained features")
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  results <- list()
  for (nm in names(contrast_defs)) {
    cm <- make_contrast(contrast_defs[[nm]], colnames(design))
    if (is.null(cm)) next
    qlf <- glmQLFTest(fit, contrast = cm[, 1])
    tab <- as.data.table(topTags(qlf, n = Inf, sort.by = "none")$table,
                        keep.rownames = "feature")
    tab[, `:=`(contrast = nm, retained_features = nrow(tab))]
    results[[nm]] <- tab
  }
  list(tables = results, y = y, design = design)
}

run_fgsea <- function(tab) {
  ranks <- sign(tab$logFC) * sqrt(pmax(tab$F, 0))
  names(ranks) <- tab$feature
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  local <- lapply(pathways, function(x) intersect(x, names(ranks)))
  local <- local[lengths(local) >= 10L]
  ans <- as.data.table(suppressWarnings(fgseaMultilevel(
    pathways = local, stats = ranks, minSize = 10L, maxSize = 2000L,
    eps = 1e-10, scoreType = "std"
  )))
  if (!nrow(ans)) return(ans)
  ans[, leadingEdge := vapply(leadingEdge, paste, collapse = ";",
                               FUN.VALUE = character(1))]
  setnames(ans, "pathway", "family_id")
  merge(ans, family_labels, by = "family_id", all.x = TRUE)
}

run_rna_models <- function(technical = FALSE, age = FALSE, suffix = "primary") {
  family_out <- list(); gene_out <- list(); model_qc <- list(); k <- 0L
  for (ct in names(cell_type_map)) {
    rna_ct <- load_counts("RNA", ct)
    m_all <- meta[cell_type == ct]
    idx <- which(m_all$n_nuclei >= 30)
    m <- copy(m_all[idx])
    if (age) {
      ok <- is.finite(m$age)
      m <- m[ok]; idx <- idx[ok]
    }
    group_n <- m[, .(n = uniqueN(donor_id)), by = diagnosis]
    if (!all(c("HC", "ALS", "ALS_FTD") %chin% as.character(group_n$diagnosis)) ||
        min(group_n$n) < 3L) next
    fit <- fit_count_model(
      rna_ct[, idx, drop = FALSE], m, primary_contrasts,
      group_var = "diagnosis",
      technical_variable = if (technical) "log10_RNA_UMI" else NULL,
      age = age
    )
    k <- k + 1L
    model_qc[[k]] <- data.table(
      analysis = suffix, cell_type = ct, samples = nrow(m), donors = uniqueN(m$donor_id),
      design_rank = qr(fit$design)$rank, design_columns = ncol(fit$design),
      retained_genes = nrow(fit$y)
    )
    for (nm in names(fit$tables)) {
      tab <- fit$tables[[nm]]
      fg <- run_fgsea(tab)
      if (nrow(fg)) fg[, `:=`(cell_type = ct, cell_type_label = cell_type_map[[ct]],
                               contrast = nm, analysis = suffix,
                               n_donors = uniqueN(m$donor_id))]
      family_out[[length(family_out) + 1L]] <- fg
      cand <- tab[feature %chin% candidate_genes]
      cand[, `:=`(gene = feature, cell_type = ct,
                  cell_type_label = cell_type_map[[ct]], analysis = suffix,
                  n_donors = uniqueN(m$donor_id))]
      gene_out[[length(gene_out) + 1L]] <- cand
    }
  }
  list(family = rbindlist(family_out, fill = TRUE),
       gene = rbindlist(gene_out, fill = TRUE), qc = rbindlist(model_qc, fill = TRUE))
}

resume_atac <- identical(Sys.getenv("RUF_RESUME_ATAC"), "1") ||
               identical(Sys.getenv("RUF_RESUME_FINAL"), "1")
resume_final <- identical(Sys.getenv("RUF_RESUME_FINAL"), "1")
if (!resume_atac) {
message_log("Running frozen RNA programme and candidate-gene models")
rna_primary <- run_rna_models(FALSE, FALSE, "primary_sex_adjusted")
rna_technical <- run_rna_models(TRUE, FALSE, "technical_depth_adjusted")
rna_age <- run_rna_models(FALSE, TRUE, "age_complete_case")

rna_family <- rna_primary$family
rna_family[, comparison_FDR := safe_bh(pval), by = .(cell_type, contrast)]
rna_family[, global_fixed_family_FDR := safe_bh(pval)]
tech_family <- rna_technical$family[, .(
  cell_type, contrast, family_id, technical_NES = NES,
  technical_p = pval, technical_padj = padj
)]
age_family <- rna_age$family[, .(
  cell_type, contrast, family_id, age_complete_case_NES = NES,
  age_complete_case_p = pval, age_complete_case_padj = padj
)]
rna_family <- merge(rna_family, tech_family,
                    by = c("cell_type", "contrast", "family_id"), all.x = TRUE)
rna_family <- merge(rna_family, age_family,
                    by = c("cell_type", "contrast", "family_id"), all.x = TRUE)
rna_family[, technical_direction_stable :=
             !is.na(technical_NES) & sign(NES) == sign(technical_NES)]
rna_family[, age_direction_stable :=
             is.na(age_complete_case_NES) | sign(NES) == sign(age_complete_case_NES)]
setorder(rna_family, global_fixed_family_FDR, cell_type, contrast, family_id)
fwrite(rna_family, file.path(out, "RNA_7_family_results.tsv"), sep = "\t")

rna_gene <- rna_primary$gene
rna_gene[, within_17_FDR := safe_bh(PValue), by = .(cell_type, contrast)]
rna_gene[, global_17_axis_FDR := safe_bh(PValue)]
rna_gene[, ci_low := logFC - qnorm(.975) * sqrt(pmax(0, logCPM * 0 + 1/F))]
# Exact QL standard errors are not exported by topTags; CI columns above are
# removed rather than exposing a pseudo-CI.
rna_gene[, c("ci_low") := NULL]
tech_gene <- rna_technical$gene[, .(
  cell_type, contrast, gene, technical_logFC = logFC,
  technical_PValue = PValue, technical_full_FDR = FDR
)]
age_gene <- rna_age$gene[, .(
  cell_type, contrast, gene, age_complete_case_logFC = logFC,
  age_complete_case_PValue = PValue, age_complete_case_full_FDR = FDR
)]
rna_gene <- merge(rna_gene, tech_gene,
                  by = c("cell_type", "contrast", "gene"), all.x = TRUE)
rna_gene <- merge(rna_gene, age_gene,
                  by = c("cell_type", "contrast", "gene"), all.x = TRUE)
rna_gene[, technical_direction_stable :=
           !is.na(technical_logFC) & sign(logFC) == sign(technical_logFC)]
rna_gene[, age_direction_stable :=
           is.na(age_complete_case_logFC) | sign(logFC) == sign(age_complete_case_logFC)]
setorder(rna_gene, FDR, cell_type, contrast, gene)
fwrite(rna_gene, file.path(out, "RNA_17_gene_results.tsv"), sep = "\t")
fwrite(rbindlist(list(rna_primary$qc, rna_technical$qc, rna_age$qc), fill = TRUE),
       file.path(out, "RNA_model_QC.tsv"), sep = "\t")

# Low-power C9 ALS-FTD sensitivity is deliberately separate from the main
# all-ALS-FTD estimand.
c9_family <- list(); c9_gene <- list()
for (ct in names(cell_type_map)) {
  rna_ct <- load_counts("RNA", ct)
  m_all <- meta[cell_type == ct]
  idx <- which(m_all$n_nuclei >= 30 &
                 as.character(m_all$case_type) %chin% c("HC", "C9_ALS_FTD"))
  m <- droplevels(copy(m_all[idx]))
  if (uniqueN(m$donor_id[m$case_type == "C9_ALS_FTD"]) < 3L) next
  fit <- fit_count_model(rna_ct[, idx, drop = FALSE], m, c9_contrast,
                         group_var = "case_type")
  tab <- fit$tables[[1L]]
  fg <- run_fgsea(tab)
  if (nrow(fg)) fg[, `:=`(cell_type = ct, contrast = names(c9_contrast),
                           n_C9_ALSFTD = uniqueN(m$donor_id[m$case_type == "C9_ALS_FTD"]),
                           n_HC = uniqueN(m$donor_id[m$case_type == "HC"]))]
  c9_family[[length(c9_family) + 1L]] <- fg
  cg <- tab[feature %chin% candidate_genes]
  cg[, `:=`(gene = feature, cell_type = ct, contrast = names(c9_contrast))]
  c9_gene[[length(c9_gene) + 1L]] <- cg
}
c9_family <- rbindlist(c9_family, fill = TRUE)
if (nrow(c9_family)) c9_family[, global_sensitivity_FDR := safe_bh(pval)]
c9_gene <- rbindlist(c9_gene, fill = TRUE)
if (nrow(c9_gene)) c9_gene[, global_17_sensitivity_FDR := safe_bh(PValue)]
fwrite(c9_family, file.path(out, "RNA_C9_ALSFTD_sensitivity_7_family.tsv"), sep = "\t")
fwrite(c9_gene, file.path(out, "RNA_C9_ALSFTD_sensitivity_17_gene.tsv"), sep = "\t")

# Full donor-deletion re-fit for the six cell types.  The full-count model and
# FGSEA are recomputed after each deletion; this is not a score-only shortcut.
message_log("Running donor-deletion RNA refits")
loo_long <- list(); loo_k <- 0L
for (ct in names(cell_type_map)) {
  rna_ct <- load_counts("RNA", ct)
  m_all <- meta[cell_type == ct]
  base_idx <- which(m_all$n_nuclei >= 30)
  donors <- unique(m_all$donor_id[base_idx])
  for (drop_donor in donors) {
    idx <- base_idx[m_all$donor_id[base_idx] != drop_donor]
    m <- droplevels(copy(m_all[idx]))
    attempt <- try(fit_count_model(rna_ct[, idx, drop = FALSE], m, primary_contrasts,
                                   group_var = "diagnosis"), silent = TRUE)
    if (inherits(attempt, "try-error")) next
    for (nm in names(attempt$tables)) {
      tab <- attempt$tables[[nm]]
      fg <- run_fgsea(tab)
      if (nrow(fg)) {
        loo_k <- loo_k + 1L
        fg[, `:=`(cell_type = ct, contrast = nm, omitted_donor = drop_donor,
                  result_type = "family", feature = family_id, effect = NES)]
        loo_long[[loo_k]] <- fg[, .(cell_type, contrast, omitted_donor,
                                     result_type, feature, effect)]
      }
      cg <- tab[feature %chin% candidate_genes]
      if (nrow(cg)) {
        loo_k <- loo_k + 1L
        cg[, `:=`(cell_type = ct, omitted_donor = drop_donor,
                  result_type = "candidate_gene", effect = logFC)]
        loo_long[[loo_k]] <- cg[, .(cell_type, contrast, omitted_donor,
                                     result_type, feature, effect)]
      }
    }
  }
  message_log("LOO completed for ", ct)
}
loo_long <- rbindlist(loo_long, fill = TRUE)
full_family_effect <- rna_family[, .(cell_type, contrast, result_type = "family",
                                     feature = family_id, full_effect = NES)]
full_gene_effect <- rna_gene[, .(cell_type, contrast, result_type = "candidate_gene",
                                 feature = gene, full_effect = logFC)]
full_effects <- rbindlist(list(full_family_effect, full_gene_effect))
loo_long <- merge(loo_long, full_effects,
                  by = c("cell_type", "contrast", "result_type", "feature"), all.x = TRUE)
loo_summary <- loo_long[, .(
  valid_deletions = .N,
  direction_retained = sum(sign(effect) == sign(full_effect), na.rm = TRUE),
  direction_fraction = mean(sign(effect) == sign(full_effect), na.rm = TRUE),
  min_effect = min(effect, na.rm = TRUE), max_effect = max(effect, na.rm = TRUE)
), by = .(cell_type, contrast, result_type, feature, full_effect)]
loo_summary[, stable_90pct := direction_fraction >= 0.90]
fwrite(loo_summary, file.path(out, "LOO_stability_summary.tsv"), sep = "\t")
} else {
  message_log("Resuming after completed RNA and LOO; loading frozen outputs")
  rna_family <- fread(file.path(out, "RNA_7_family_results.tsv"))
  rna_gene <- fread(file.path(out, "RNA_17_gene_results.tsv"))
  loo_summary <- fread(file.path(out, "LOO_stability_summary.tsv"))
  stopifnot(nrow(rna_family) > 0L, nrow(rna_gene) > 0L, nrow(loo_summary) > 0L)
}

# Map public peaks to a frozen ±2 kb promoter window around the 17 genes.
if (!resume_final) {
message_log("Mapping candidate promoters and fitting peak-wide ATAC models")
addArchRGenome("hg38")
ga <- as.data.table(as.data.frame(getGeneAnnotation()$genes))
ga <- ga[symbol %chin% candidate_genes]
ga[, tss := fifelse(strand == "+", start, end)]
atac_features <- fread(file.path(root,
  "data/raw/Ruf_2026_multiome/Multiome_Dataset_ATAC_features.tsv"),
  header = FALSE)[[1]]
peak_dt <- data.table(peak = atac_features)
peak_dt[, c("chr", "start", "end") := tstrsplit(peak, "-", fixed = TRUE,
                                                  type.convert = TRUE)]
promoter_pairs <- rbindlist(lapply(seq_len(nrow(ga)), function(i) {
  hit <- peak_dt[chr == ga$seqnames[i] & end >= ga$tss[i] - 2000L &
                   start <= ga$tss[i] + 2000L]
  if (!nrow(hit)) return(NULL)
  hit[, .(gene = ga$symbol[i], peak, chr, start, end,
           tss = ga$tss[i], strand = as.character(ga$strand[i]))]
}))
promoter_pairs <- unique(promoter_pairs, by = c("gene", "peak"))
fwrite(promoter_pairs, file.path(out, "frozen_17_gene_promoter_peak_map.tsv"), sep = "\t")

fit_atac_low_memory <- function(logcpm, m, technical_variable = NULL) {
  design <- make_design(m, "diagnosis", technical_variable = technical_variable)
  cm_list <- lapply(primary_contrasts, make_contrast, design_cols = colnames(design))
  cm_list <- cm_list[!vapply(cm_list, is.null, logical(1))]
  cm <- do.call(cbind, lapply(cm_list, function(x) x[, 1]))
  colnames(cm) <- names(cm_list)
  fit <- lmFit(logcpm, design)
  fit <- contrasts.fit(fit, cm)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  output <- list()
  target_rows <- which(rownames(logcpm) %chin% promoter_pairs$peak)
  for (j in seq_len(ncol(cm))) {
    all_p <- fit$p.value[, j]
    all_fdr <- safe_bh(all_p)
    output[[colnames(cm)[j]]] <- data.table(
      feature = rownames(logcpm)[target_rows],
      logFC = fit$coefficients[target_rows, j],
      AveExpr = fit$Amean[target_rows],
      t = fit$t[target_rows, j],
      PValue = all_p[target_rows],
      FDR = all_fdr[target_rows],
      retained_features = nrow(logcpm)
    )
  }
  rm(fit); gc()
  output
}

atac_primary <- list(); atac_technical <- list(); atac_logcpm_by_cell <- list()
for (ct in names(cell_type_map)) {
  message_log("ATAC low-memory full-peak model started for ", ct)
  atac_ct <- load_counts("ATAC", ct)
  m_all <- meta[cell_type == ct]
  idx <- which(m_all$n_nuclei >= 30)
  m <- copy(m_all[idx])
  design_filter <- make_design(m, "diagnosis")
  y <- DGEList(atac_ct[, idx, drop = FALSE])
  keep <- filterByExpr(y, design = design_filter)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y)
  atac_logcpm <- cpm(y, log = TRUE, prior.count = 2)
  rm(y, atac_ct); gc()
  primary_tables <- fit_atac_low_memory(atac_logcpm, m)
  technical_tables <- fit_atac_low_memory(
    atac_logcpm, m, technical_variable = "log10_ATAC_fragments"
  )
  target_present <- intersect(promoter_pairs$peak, rownames(atac_logcpm))
  atac_logcpm_by_cell[[ct]] <- atac_logcpm[target_present, , drop = FALSE]
  for (nm in names(primary_tables)) {
    tab <- primary_tables[[nm]]
    tab <- merge(tab, promoter_pairs, by.x = "feature", by.y = "peak",
                 all.x = TRUE, allow.cartesian = TRUE)
    tab[, `:=`(peak = feature, cell_type = ct, contrast = nm,
               n_donors = uniqueN(m$donor_id))]
    atac_primary[[length(atac_primary) + 1L]] <- tab
    tt <- technical_tables[[nm]][, .(
      peak = feature, technical_logFC = logFC,
      technical_PValue = PValue, technical_peak_wide_FDR = FDR
    )]
    tt[, `:=`(cell_type = ct, contrast = nm)]
    atac_technical[[length(atac_technical) + 1L]] <- tt
  }
  rm(atac_logcpm, primary_tables, technical_tables); gc()
  message_log("ATAC low-memory full-peak model completed for ", ct)
}
atac_result <- rbindlist(atac_primary, fill = TRUE)
atac_result[, target_promoter_FDR := safe_bh(PValue), by = .(cell_type, contrast)]
atac_result[, global_target_promoter_FDR := safe_bh(PValue)]
atac_technical <- rbindlist(atac_technical, fill = TRUE)
atac_result <- merge(atac_result, atac_technical,
                     by = c("peak", "cell_type", "contrast"), all.x = TRUE)
atac_result[, technical_direction_stable :=
              !is.na(technical_logFC) & sign(logFC) == sign(technical_logFC)]
setorder(atac_result, FDR, cell_type, contrast, gene, peak)
fwrite(atac_result, file.path(out, "ATAC_targeted_results.tsv"), sep = "\t")

# Same-donor residual RNA--promoter accessibility association.  This is an
# association test, not peak-to-gene causality.
association <- list()
# Only frozen genes and their promoter peaks are materialized as dense logCPM
# matrices, avoiding an unnecessary ~200,000 x sample dense ATAC object.
for (ct in names(cell_type_map)) {
  rna_ct <- load_counts("RNA", ct)
  m_all <- meta[cell_type == ct]
  idx <- which(m_all$n_nuclei >= 30)
  m <- copy(m_all[idx])
  rna_assoc_features <- intersect(candidate_genes, rownames(rna_ct))
  rna_logcpm <- cpm(DGEList(rna_ct[rna_assoc_features, idx, drop = FALSE]),
                    log = TRUE, prior.count = 2)
  cov_design <- model.matrix(~ diagnosis + sex + scale(log10_RNA_UMI) +
                               scale(log10_ATAC_fragments), m)
  for (i in seq_len(nrow(promoter_pairs))) {
    gene <- promoter_pairs$gene[i]; peak <- promoter_pairs$peak[i]
    atac_local <- atac_logcpm_by_cell[[ct]]
    if (!gene %chin% rownames(rna_logcpm) || !peak %chin% rownames(atac_local)) next
    x <- as.numeric(rna_logcpm[gene, ]); y <- as.numeric(atac_local[peak, ])
    rx <- residuals(lm.fit(cov_design, x)); ry <- residuals(lm.fit(cov_design, y))
    test <- suppressWarnings(cor.test(rx, ry, method = "pearson"))
    association[[length(association) + 1L]] <- data.table(
      cell_type = ct, gene = gene, peak = peak, n_donors = length(rx),
      residual_correlation = unname(test$estimate), association_p = test$p.value
    )
  }
}
association <- rbindlist(association, fill = TRUE)
association[, association_FDR := safe_bh(association_p), by = cell_type]
fwrite(association, file.path(out, "RNA_ATAC_residual_associations.tsv"), sep = "\t")

joint <- merge(
  rna_gene[, .(cell_type, contrast, gene, RNA_logFC = logFC,
                RNA_full_FDR = FDR, RNA_within17_FDR = within_17_FDR,
                RNA_technical_stable = technical_direction_stable)],
  atac_result[, .(cell_type, contrast, gene, peak, ATAC_logFC = logFC,
                   ATAC_peak_wide_FDR = FDR,
                   ATAC_target_promoter_FDR = target_promoter_FDR,
                   ATAC_technical_stable = technical_direction_stable)],
  by = c("cell_type", "contrast", "gene"), allow.cartesian = TRUE
)
joint <- merge(joint, association,
               by = c("cell_type", "gene", "peak"), all.x = TRUE)
gene_loo <- loo_summary[result_type == "candidate_gene",
                        .(cell_type, contrast, gene = feature,
                          RNA_LOO_direction_fraction = direction_fraction,
                          RNA_LOO_stable_90pct = stable_90pct)]
joint <- merge(joint, gene_loo, by = c("cell_type", "contrast", "gene"), all.x = TRUE)
joint[, concordant_direction := sign(RNA_logFC) == sign(ATAC_logFC)]
joint[, complete_joint_chain :=
        RNA_full_FDR < 0.05 & ATAC_peak_wide_FDR < 0.05 &
        association_FDR < 0.05 & concordant_direction &
        RNA_technical_stable & ATAC_technical_stable &
        RNA_LOO_stable_90pct]
joint[, interpretation := fifelse(
  complete_joint_chain,
  "strict same-nucleus donor-level RNA-promoter joint support",
  "incomplete joint evidence; no regulatory or causal claim"
)]
setorder(joint, -complete_joint_chain, RNA_full_FDR, ATAC_peak_wide_FDR,
         association_FDR)
fwrite(joint, file.path(out, "RNA_ATAC_joint_results.tsv"), sep = "\t")
} else {
  message_log("Resuming final integration from completed ATAC and joint outputs")
  atac_result <- fread(file.path(out, "ATAC_targeted_results.tsv"))
  joint <- fread(file.path(out, "RNA_ATAC_joint_results.tsv"))
  promoter_pairs <- fread(file.path(out, "frozen_17_gene_promoter_peak_map.tsv"))
  gene_loo <- loo_summary[result_type == "candidate_gene",
                          .(cell_type, contrast, gene = feature,
                            RNA_LOO_direction_fraction = direction_fraction,
                            RNA_LOO_stable_90pct = stable_90pct)]
  stopifnot(nrow(atac_result) > 0L, nrow(joint) > 0L, nrow(promoter_pairs) > 0L)
}

strict_family <- merge(
  rna_family,
  loo_summary[result_type == "family",
              .(cell_type, contrast, family_id = feature,
                LOO_direction_fraction = direction_fraction,
                LOO_stable_90pct = stable_90pct)],
  by = c("cell_type", "contrast", "family_id"), all.x = TRUE
)
strict_family[, strict_external_programme_replication :=
                global_fixed_family_FDR < 0.05 & LOO_stable_90pct &
                technical_direction_stable & age_direction_stable]

# Test the previously frozen ALS-vs-control programme directions from the two
# earlier cohorts.  This is the external-replication decision table; Ruf effect
# signs are not used to redefine the expected direction.
prior_discovery <- fread(file.path(
  root, "results/submission_cross_cohort_evidence_v1/exact_estimand_program_replication.tsv"
))
prior_discovery <- prior_discovery[estimand == "ALS_vs_Control",
  .(contrast = estimand, cell_type, family_id,
    GSE219280_NES = effect_discovery,
    GSE219280_global_FDR = adjusted_p_discovery,
    Gittings_NES = effect_replication,
    Gittings_global_FDR = adjusted_p_replication,
    prior_exact_direction_concordant = same_direction,
    prior_both_global_FDR = both_global_FDR_0_05,
    prior_replication_class = replication_class)]
setnames(prior_discovery, "cell_type", "cell_type_label")
ruf_validation <- strict_family[contrast == "ALS_vs_Control",
  .(cell_type, cell_type_label, contrast, family_id, Ruf_NES = NES,
    Ruf_global_FDR = global_fixed_family_FDR,
    Ruf_LOO_direction_fraction = LOO_direction_fraction,
    Ruf_technical_stable = technical_direction_stable,
    Ruf_age_stable = age_direction_stable,
    Ruf_strict_internal_gate = strict_external_programme_replication)]
external_replication <- merge(
  prior_discovery, ruf_validation,
  by = c("cell_type_label", "contrast", "family_id"), all = TRUE
)
external_replication[, Ruf_vs_GSE219280_direction_concordant :=
                       !is.na(Ruf_NES) & !is.na(GSE219280_NES) &
                       sign(Ruf_NES) == sign(GSE219280_NES)]
external_replication[, Ruf_vs_Gittings_direction_concordant :=
                       !is.na(Ruf_NES) & !is.na(Gittings_NES) &
                       sign(Ruf_NES) == sign(Gittings_NES)]
external_replication[, strict_Ruf_replication_of_discovery_direction :=
                       Ruf_strict_internal_gate &
                       Ruf_vs_GSE219280_direction_concordant]
external_replication[, strict_three_cohort_same_direction :=
                       Ruf_strict_internal_gate & prior_both_global_FDR &
                       prior_exact_direction_concordant &
                       Ruf_vs_GSE219280_direction_concordant &
                       Ruf_vs_Gittings_direction_concordant]
external_replication[, interpretation := fcase(
  strict_three_cohort_same_direction,
  "strict same-direction three-cohort programme replication",
  strict_Ruf_replication_of_discovery_direction,
  "Ruf validates the discovery direction; earlier external cohort may be heterogeneous",
  Ruf_strict_internal_gate,
  "Ruf-significant and internally robust but direction does not replicate discovery",
  default = "not strictly replicated in Ruf"
)]
fwrite(external_replication,
       file.path(out, "three_cohort_ALS_programme_replication.tsv"), sep = "\t")
fwrite(strict_family, file.path(out, "RNA_7_family_strict_replication.tsv"), sep = "\t")

qc <- data.table(
  metric = c(
    "donors", "nuclei", "cell_types", "frozen_families", "frozen_genes",
    "primary_family_axes", "strict_programme_axes",
    "strict_three_cohort_same_direction_axes", "candidate_gene_axes",
    "strict_single_gene_replications", "promoter_peaks_tested",
    "peak_wide_FDR_target_hits", "complete_RNA_ATAC_joint_chains"
  ),
  value = c(
    uniqueN(meta$donor_id), sum(meta$n_nuclei), uniqueN(meta$cell_type),
    length(pathways), length(candidate_genes), nrow(rna_family),
    sum(strict_family$strict_external_programme_replication, na.rm = TRUE),
    sum(external_replication$strict_three_cohort_same_direction, na.rm = TRUE),
    nrow(rna_gene),
    sum(rna_gene$FDR < 0.05 & rna_gene$technical_direction_stable &
          rna_gene$age_direction_stable &
          paste(rna_gene$cell_type, rna_gene$contrast, rna_gene$gene) %chin%
          paste(gene_loo$cell_type, gene_loo$contrast, gene_loo$gene)[gene_loo$RNA_LOO_stable_90pct],
        na.rm = TRUE),
    uniqueN(promoter_pairs$peak), sum(atac_result$FDR < 0.05, na.rm = TRUE),
    sum(joint$complete_joint_chain, na.rm = TRUE)
  )
)
fwrite(qc, file.path(out, "final_QC.tsv"), sep = "\t")

claim_limits <- data.table(
  topic = c("independence", "C9 phenotype", "clinical covariates", "RNA_ATAC",
            "candidate genes"),
  permitted_claim = c(
    "external public cohort; no public donor-ID overlap detected",
    "C9 ALS-FTD versus control is a seven-donor sensitivity analysis",
    "age/PMI are complete-case sensitivities because nine matrix IDs lack an exact public clinical-table match",
    "same-nucleus donor-level association only when all predeclared gates pass",
    "17 genes were discovery-derived and frozen before Ruf external testing"
  ),
  prohibited_claim = c(
    "proven donor independence without a cross-study private crosswalk",
    "independent C9 ALS versus C9 FTD phenotypic-switch replication",
    "imputed or guessed clinical-ID matches",
    "ATAC-to-RNA causality or TF activity",
    "outcome-independent discovery or de novo biomarker selection"
  )
)
fwrite(claim_limits, file.path(out, "claim_boundaries.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "inferential_sessionInfo.txt"))
message_log("Ruf 2026 frozen validation complete")
