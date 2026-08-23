#!/usr/bin/env Rscript

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
  library(nnls)
})

set.seed(20260810)
sn_path <- file.path(project_root, "data/processed/GSE219280_major_celltype_pseudobulk.rds")
bulk_path <- file.path(project_root, "data/processed/GSE153960_subject_region_counts.rds")
map_path <- file.path(project_root, "data/processed/gencode.v35.gene_id_map.tsv")
primary_path <- file.path(project_root, "results/GSE153960_bulk_DE/all_gene_contrast_results.tsv")
out_dir <- file.path(project_root, "results/GSE153960_composition_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sn_obj <- readRDS(sn_path)
sn_counts <- sn_obj$counts
sn_meta <- as.data.table(sn_obj$samples)
bulk_obj <- readRDS(bulk_path)
bulk_counts <- bulk_obj$counts
bulk_meta <- as.data.table(bulk_obj$samples)
gene_map <- fread(map_path)
primary <- fread(primary_path)

stopifnot(identical(colnames(sn_counts), sn_meta$pseudobulk_id))
stopifnot(identical(colnames(bulk_counts), bulk_meta$subject_region_id))

sn_meta[, broad_cell_type := fcase(
  grepl("^Exc_", annotation_major_cell_type), "Excitatory",
  grepl("^Inh_", annotation_major_cell_type), "Inhibitory",
  annotation_major_cell_type == "Astro", "Astrocyte",
  annotation_major_cell_type == "Oligo", "Oligodendrocyte",
  annotation_major_cell_type == "OPC", "OPC",
  annotation_major_cell_type == "Micro", "Microglia",
  annotation_major_cell_type %chin% c("VLMC", "Endo"), "Vascular",
  default = NA_character_
)]
selected <- which(!is.na(sn_meta$broad_cell_type))
sn_meta_sel <- copy(sn_meta[selected])
sn_counts_sel <- sn_counts[, selected, drop = FALSE]
sn_meta_sel[, aggregate_id := paste(sample, broad_cell_type, sep = "__")]
aggregate_ids <- unique(sn_meta_sel$aggregate_id)

agg_counts <- do.call(cbind, lapply(aggregate_ids, function(id) {
  Matrix::rowSums(sn_counts_sel[, sn_meta_sel$aggregate_id == id, drop = FALSE])
}))
colnames(agg_counts) <- aggregate_ids
rownames(agg_counts) <- rownames(sn_counts_sel)
agg_meta <- sn_meta_sel[match(aggregate_ids, aggregate_id), .(
  aggregate_id, sample, donor_id, diagnosis, brain_region, broad_cell_type
)]

control_idx <- which(agg_meta$diagnosis == "Control")
ref_counts <- agg_counts[, control_idx, drop = FALSE]
ref_meta <- agg_meta[control_idx]
ref_y <- DGEList(counts = ref_counts)
ref_y <- normLibSizes(ref_y, method = "TMM")
ref_logcpm <- cpm(ref_y, log = TRUE, prior.count = 1)
ref_cpm <- cpm(ref_y, log = FALSE, prior.count = 0.25)

cell_types <- c("Excitatory", "Inhibitory", "Astrocyte", "Oligodendrocyte", "OPC", "Microglia", "Vascular")
mean_log <- sapply(cell_types, function(ct) rowMeans(ref_logcpm[, ref_meta$broad_cell_type == ct, drop = FALSE]))
mean_cpm <- sapply(cell_types, function(ct) rowMeans(ref_cpm[, ref_meta$broad_cell_type == ct, drop = FALSE]))
colnames(mean_log) <- colnames(mean_cpm) <- cell_types

marker_list <- lapply(cell_types, function(ct) {
  others <- setdiff(cell_types, ct)
  specificity <- mean_log[, ct] - apply(mean_log[, others, drop = FALSE], 1L, max)
  dt <- data.table(
    gene = rownames(mean_log),
    cell_type = ct,
    mean_logCPM = mean_log[, ct],
    specificity = specificity
  )
  dt[mean_logCPM > 1 & specificity > 0.5][order(-specificity, -mean_logCPM)][1:min(.N, 75L)]
})
markers <- rbindlist(marker_list)
fwrite(markers, file.path(out_dir, "snRNA_control_reference_markers.tsv"), sep = "\t", quote = FALSE)

# Exact GENCODE versioned mapping with unversioned fallback.
map_v <- unique(gene_map[, .(
  ensembl_id_versioned = gene_id_versioned,
  gene_name,
  gene_type
)], by = "ensembl_id_versioned")
bulk_ids <- data.table(
  row_index = seq_len(nrow(bulk_counts)),
  ensembl_id_versioned = rownames(bulk_counts),
  ensembl_id = sub("\\..*$", "", rownames(bulk_counts))
)
bulk_ids <- map_v[bulk_ids, on = "ensembl_id_versioned"]
missing <- is.na(bulk_ids$gene_name)
if (any(missing)) {
  map_u <- unique(gene_map[, .(ensembl_id = gene_id, fallback_gene_name = gene_name, fallback_gene_type = gene_type)], by = "ensembl_id")
  fallback <- map_u[bulk_ids[missing], on = "ensembl_id"]
  bulk_ids[missing, `:=`(gene_name = fallback$fallback_gene_name, gene_type = fallback$fallback_gene_type)]
}

mapped <- bulk_ids[!is.na(gene_name) & nzchar(gene_name) & gene_type == "protein_coding"]
bulk_gene_counts <- rowsum(bulk_counts[mapped$row_index, , drop = FALSE], group = mapped$gene_name, reorder = FALSE)
bulk_y_gene <- DGEList(counts = bulk_gene_counts)
bulk_y_gene <- normLibSizes(bulk_y_gene, method = "TMM")
bulk_cpm_gene <- cpm(bulk_y_gene, log = FALSE, prior.count = 0.25)

marker_union <- unique(markers$gene)
common_markers <- intersect(marker_union, intersect(rownames(mean_cpm), rownames(bulk_cpm_gene)))
if (length(common_markers) < 100L) stop("Too few cross-platform marker genes: ", length(common_markers))

signature <- mean_cpm[common_markers, cell_types, drop = FALSE]
bulk_marker <- bulk_cpm_gene[common_markers, , drop = FALSE]
row_scale <- pmax(rowMeans(signature), 0.5)
A <- signature / row_scale
B <- bulk_marker / row_scale

fit_one <- function(j) {
  fit <- nnls::nnls(A, B[, j])
  coef <- pmax(coef(fit), 0)
  fractions <- if (sum(coef) > 0) coef / sum(coef) else rep(NA_real_, length(coef))
  residual <- sqrt(sum((A %*% coef - B[, j])^2)) / sqrt(sum(B[, j]^2))
  c(fractions, residual = residual)
}
fraction_matrix <- vapply(seq_len(ncol(B)), fit_one, numeric(length(cell_types) + 1L))
fraction_dt <- as.data.table(t(fraction_matrix))
setnames(fraction_dt, c(cell_types, "relative_residual"))
fraction_dt[, subject_region_id := colnames(B)]
setcolorder(fraction_dt, c("subject_region_id", cell_types, "relative_residual"))
fraction_dt <- bulk_meta[fraction_dt, on = "subject_region_id"]
fwrite(fraction_dt, file.path(out_dir, "GSE153960_snRNA_reference_NNLS_composition.tsv"), sep = "\t", quote = FALSE)

fraction_summary <- fraction_dt[, lapply(.SD, median, na.rm = TRUE), by = .(diagnosis_short, region_short), .SDcols = cell_types]
fwrite(fraction_summary, file.path(out_dir, "composition_group_medians.tsv"), sep = "\t", quote = FALSE)

# Use CLR-PCA coordinates as composition covariates. This avoids singularity
# from directly entering seven fractions that sum to one.
frac_mat <- as.matrix(fraction_dt[, ..cell_types])
clr <- log(frac_mat + 1e-4)
clr <- clr - rowMeans(clr)
pca <- prcomp(clr, center = TRUE, scale. = FALSE)
variance <- pca$sdev^2 / sum(pca$sdev^2)
k <- min(4L, max(2L, which(cumsum(variance) >= 0.8)[1L]))
pc_scores <- as.data.table(pca$x[, seq_len(k), drop = FALSE])
setnames(pc_scores, paste0("composition_PC", seq_len(k)))
meta <- cbind(copy(bulk_meta), pc_scores)

meta[, condition_region := factor(
  paste(diagnosis_short, region_short, sep = "_"),
  levels = c("Control_Frontal", "Control_Motor", "ALS_Frontal", "ALS_Motor")
)]
meta[, project := factor(project)]
meta[, library_preparation_method := factor(library_preparation_method)]
pc_terms <- paste(names(pc_scores), collapse = " + ")
design <- model.matrix(
  as.formula(paste("~ 0 + condition_region + project + library_preparation_method +", pc_terms)),
  data = meta
)
colnames(design) <- make.names(sub("^condition_region", "", colnames(design)), unique = TRUE)
stopifnot(qr(design)$rank == ncol(design))

y <- DGEList(counts = bulk_counts)
keep <- filterByExpr(y, group = meta$condition_region, min.count = 10, min.total.count = 20)
y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
v0 <- voom(y, design, plot = FALSE)
corfit <- duplicateCorrelation(v0, design, block = meta$subject_id)
consensus <- corfit$consensus.correlation
if (!is.finite(consensus)) consensus <- 0
v <- voom(y, design, plot = FALSE, block = meta$subject_id, correlation = consensus)
fit <- lmFit(v, design, block = meta$subject_id, correlation = consensus)

contrast_defs <- c(
  ALS_vs_Control_mean = "(ALS_Frontal + ALS_Motor - Control_Frontal - Control_Motor) / 2",
  ALS_region_interaction = "(ALS_Motor - ALS_Frontal) - (Control_Motor - Control_Frontal)",
  ALS_vs_Control_Motor = "ALS_Motor - Control_Motor",
  ALS_vs_Control_Frontal = "ALS_Frontal - Control_Frontal"
)
cm <- makeContrasts(contrasts = unname(contrast_defs), levels = design)
colnames(cm) <- names(contrast_defs)
fitc <- contrasts.fit(fit, cm)
eb <- eBayes(fitc, robust = TRUE)
tr <- treat(fitc, lfc = log2(1.2), robust = TRUE)

adjusted <- rbindlist(lapply(colnames(cm), function(cc) {
  x <- topTable(eb, coef = cc, number = Inf, sort.by = "none")
  z <- topTable(tr, coef = cc, number = Inf, sort.by = "none")
  data.table(
    contrast = cc,
    ensembl_id_versioned = rownames(x),
    ensembl_id = sub("\\..*$", "", rownames(x)),
    log2FC = x$logFC,
    AveExpr = x$AveExpr,
    t = x$t,
    PValue = x$P.Value,
    FDR = x$adj.P.Val,
    treat_t = z$t,
    treat_PValue = z$P.Value,
    treat_FDR = z$adj.P.Val
  )
}))
adjusted[, global_FDR := p.adjust(PValue, method = "BH")]
adjusted[, global_treat_FDR := p.adjust(treat_PValue, method = "BH")]
fwrite(adjusted, file.path(out_dir, "composition_adjusted_all_gene_contrast_results.tsv"), sep = "\t", quote = FALSE)

comparison <- merge(
  primary[, .(
    contrast, ensembl_id_versioned,
    primary_log2FC = log2FC, primary_FDR = FDR, primary_treat_FDR = treat_FDR
  )],
  adjusted[, .(
    contrast, ensembl_id_versioned,
    adjusted_log2FC = log2FC, adjusted_FDR = FDR, adjusted_treat_FDR = treat_FDR
  )],
  by = c("contrast", "ensembl_id_versioned")
)
comparison[, same_direction := sign(primary_log2FC) == sign(adjusted_log2FC)]
comparison[, primary_abs_rank := frank(-abs(primary_log2FC), ties.method = "min"), by = contrast]
comparison_summary <- comparison[, .(
  genes_joined = .N,
  pearson_effect = cor(primary_log2FC, adjusted_log2FC),
  spearman_effect = cor(primary_log2FC, adjusted_log2FC, method = "spearman"),
  direction_all = mean(same_direction),
  direction_primary_top500 = mean(same_direction[primary_abs_rank <= 500]),
  primary_FDR_genes = sum(primary_FDR < 0.05),
  adjusted_FDR_genes = sum(adjusted_FDR < 0.05),
  primary_TREAT_genes = sum(primary_treat_FDR < 0.05),
  adjusted_TREAT_genes = sum(adjusted_treat_FDR < 0.05),
  primary_TREAT_direction_retained = ifelse(sum(primary_treat_FDR < 0.05) > 0, mean(same_direction[primary_treat_FDR < 0.05]), NA_real_),
  primary_TREAT_retain_adjusted_FDR = ifelse(sum(primary_treat_FDR < 0.05) > 0, mean(adjusted_FDR[primary_treat_FDR < 0.05] < 0.05), NA_real_),
  primary_TREAT_retain_adjusted_TREAT = ifelse(sum(primary_treat_FDR < 0.05) > 0, mean(adjusted_treat_FDR[primary_treat_FDR < 0.05] < 0.05), NA_real_)
), by = contrast]

fwrite(comparison, file.path(out_dir, "primary_vs_composition_adjusted_gene_effects.tsv"), sep = "\t", quote = FALSE)
fwrite(comparison_summary, file.path(out_dir, "primary_vs_composition_adjusted_summary.tsv"), sep = "\t", quote = FALSE)
fwrite(data.table(
  marker_genes = length(common_markers),
  composition_PCs = k,
  composition_variance_explained = sum(variance[seq_len(k)]),
  median_NNLS_relative_residual = median(fraction_dt$relative_residual),
  duplicate_correlation = consensus,
  genes_tested = nrow(y)
), file.path(out_dir, "composition_model_QC.tsv"), sep = "\t", quote = FALSE)

print(fraction_summary)
print(comparison_summary)
