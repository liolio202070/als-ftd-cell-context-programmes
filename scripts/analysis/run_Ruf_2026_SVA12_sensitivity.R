#!/usr/bin/env Rscript

# Reproduce the source study's 12-SV adjustment as a frozen sensitivity layer.
# The seven programmes and 17 genes remain exactly those fixed before Ruf data.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
root <- analysis_project_root()
configure_project_library(root)
.libPaths(c(file.path(root, "cache/R_library"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(fgsea)
  library(sva)
})
set.seed(20260814)

out <- file.path(root, "results/Ruf_2026_frozen_validation")
meta <- fread(file.path(out, "pseudobulk_sample_metadata.tsv"))
families <- unique(fread(file.path(
  root, "results/Gittings_2023_C9_snRNA/predefined_7_family_gene_membership.tsv"
))[, .(family_id, family_label, gene)])
pathways <- lapply(split(families$gene, families$family_id), unique)
labels <- unique(families[, .(family_id, family_label)])
genes <- unique(fread(file.path(
  root, "results/leading_edge_evidence_v2/leading_edge_driver_shortlist.tsv"
))$gene)

contrasts <- list(
  ALS_vs_Control = c(diagnosisALS = 1, diagnosisHC = -1),
  ALSFTD_vs_Control = c(diagnosisALS_FTD = 1, diagnosisHC = -1),
  ALS_vs_ALSFTD = c(diagnosisALS = 1, diagnosisALS_FTD = -1)
)

safe_bh <- function(p) {
  z <- rep(NA_real_, length(p)); ok <- is.finite(p); z[ok] <- p.adjust(p[ok], "BH"); z
}

run_fgsea <- function(effect, feature) {
  names(effect) <- feature
  effect <- sort(effect[is.finite(effect)], decreasing = TRUE)
  local <- lapply(pathways, function(x) intersect(x, names(effect)))
  local <- local[lengths(local) >= 10L]
  ans <- as.data.table(suppressWarnings(fgseaMultilevel(
    pathways = local, stats = effect, minSize = 10L, maxSize = 2000L,
    eps = 1e-10, scoreType = "std"
  )))
  ans[, leadingEdge := vapply(leadingEdge, paste, collapse = ";",
                               FUN.VALUE = character(1))]
  setnames(ans, "pathway", "family_id")
  merge(ans, labels, by = "family_id", all.x = TRUE)
}

family_out <- list(); gene_out <- list(); qc <- list()
for (ct in unique(meta$cell_type)) {
  sm <- meta[cell_type == ct & n_nuclei >= 30]
  counts <- readRDS(file.path(out, sprintf("RNA_pseudobulk_counts_%s.rds", ct)))
  counts <- counts[, sm$sample_id, drop = FALSE]
  sm[, diagnosis := factor(diagnosis, levels = c("HC", "ALS", "ALS_FTD"))]
  mod <- model.matrix(~ 0 + diagnosis + sex, sm)
  colnames(mod) <- sub("^diagnosis", "diagnosis", colnames(mod))
  mod0 <- model.matrix(~ 1 + sex, sm)
  y <- DGEList(counts)
  keep <- filterByExpr(y, design = mod)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y)
  logcpm <- cpm(y, log = TRUE, prior.count = 2)
  normalized_counts <- cpm(y, log = FALSE)
  normalized_counts <- normalized_counts[rowMeans(normalized_counts) > 1, , drop = FALSE]

  # This follows the public Ruf code's protected full model and null model,
  # with sex included in both because it is an explicit primary covariate here.
  sv <- svaseq(normalized_counts, mod, mod0, n.sv = 12)$sv
  colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
  design <- cbind(mod, sv)
  if (qr(design)$rank != ncol(design)) stop("Non-full-rank SVA design for ", ct)
  cm <- matrix(0, nrow = ncol(design), ncol = length(contrasts),
               dimnames = list(colnames(design), names(contrasts)))
  for (nm in names(contrasts)) cm[names(contrasts[[nm]]), nm] <- contrasts[[nm]]
  fit <- eBayes(contrasts.fit(lmFit(logcpm, design), cm), trend = TRUE, robust = TRUE)
  qc[[length(qc) + 1L]] <- data.table(
    cell_type = ct, samples = nrow(sm), retained_genes = nrow(logcpm),
    n_sv = ncol(sv), design_rank = qr(design)$rank,
    design_columns = ncol(design), pass = TRUE
  )
  for (j in seq_along(contrasts)) {
    nm <- names(contrasts)[j]
    ranked <- sign(fit$coefficients[, j]) * sqrt(pmax(fit$t[, j]^2, 0))
    fg <- run_fgsea(ranked, rownames(fit))
    fg[, `:=`(cell_type = ct, contrast = nm)]
    family_out[[length(family_out) + 1L]] <- fg
    idx <- match(genes, rownames(fit), nomatch = 0L); idx <- idx[idx > 0L]
    gene_out[[length(gene_out) + 1L]] <- data.table(
      cell_type = ct, contrast = nm, gene = rownames(fit)[idx],
      SVA12_logFC = fit$coefficients[idx, j],
      SVA12_t = fit$t[idx, j], SVA12_PValue = fit$p.value[idx, j],
      SVA12_full_FDR = safe_bh(fit$p.value[, j])[idx]
    )
  }
  rm(counts, y, logcpm, fit, sv); gc()
}

family <- rbindlist(family_out, fill = TRUE)
family[, SVA12_global_fixed_family_FDR := safe_bh(pval)]
setnames(family, "NES", "SVA12_NES")
gene <- rbindlist(gene_out, fill = TRUE)
primary_family <- fread(file.path(out, "RNA_7_family_strict_replication.tsv"))
primary_gene <- fread(file.path(out, "RNA_17_gene_results.tsv"))

family_compare <- merge(
  primary_family,
  family[, .(cell_type, contrast, family_id, SVA12_NES,
              SVA12_p = pval, SVA12_global_fixed_family_FDR)],
  by = c("cell_type", "contrast", "family_id"), all.x = TRUE
)
family_compare[, SVA12_direction_stable :=
                 !is.na(SVA12_NES) & sign(NES) == sign(SVA12_NES)]
family_compare[, strict_with_SVA12 :=
                 strict_external_programme_replication & SVA12_direction_stable]

gene_compare <- merge(primary_gene, gene,
                      by = c("cell_type", "contrast", "gene"), all.x = TRUE)
gene_compare[, SVA12_direction_stable :=
               !is.na(SVA12_logFC) & sign(logFC) == sign(SVA12_logFC)]

fwrite(family_compare, file.path(out, "RNA_7_family_SVA12_sensitivity.tsv"), sep = "\t")
fwrite(gene_compare, file.path(out, "RNA_17_gene_SVA12_sensitivity.tsv"), sep = "\t")
fwrite(rbindlist(qc), file.path(out, "SVA12_model_QC.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "SVA12_sessionInfo.txt"))
