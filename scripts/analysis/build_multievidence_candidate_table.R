script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
project_root <- analysis_project_root()
configure_project_library(project_root)
project_lib <- file.path(project_root, "cache/R_library")
.libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages(library(data.table))

bulk <- fread(file.path(
  project_root,
  "results/cross_dataset_validation/GSE153960_bulk_DE_GENCODE_v35_annotated.tsv"
))
sn <- fread(file.path(
  project_root,
  "results/GSE219280_pseudobulk_DE/all_gene_contrast_results.tsv"
))
published <- fread(file.path(
  project_root,
  "data/processed/GSE219280_published_MAST_DE.tsv"
))
panel <- fread(file.path(
  project_root,
  "data/processed/PanelApp_ALS_MND_panel_263_v1.75_genes.tsv"
))
tdp <- fread(file.path(
  project_root,
  "data/processed/Gittings_2023_TDP43_cryptic_exon_targets.tsv"
))

out_dir <- file.path(project_root, "results/multievidence_candidates")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# One well-expressed Ensembl feature per protein-coding symbol.
bulk_mean <- bulk[
  contrast == "ALS_vs_Control_mean" &
    gene_type == "protein_coding" &
    !is.na(gene_name) & nzchar(gene_name)
][order(gene_name, -AveExpr)][, .SD[1L], by = gene_name]
bulk_mean <- bulk_mean[, .(
  gene = gene_name,
  ensembl_id,
  bulk_log2FC = log2FC,
  bulk_t = t,
  bulk_PValue = PValue,
  bulk_FDR = FDR,
  bulk_treat_FDR = treat_FDR,
  bulk_global_FDR = global_FDR,
  bulk_global_treat_FDR = global_treat_FDR,
  bulk_AveExpr = AveExpr
)]
bulk_mean[, bulk_direction := fifelse(bulk_log2FC > 0, "up", "down")]

exc_types <- c("Exc_deep", "Exc_intermediate", "Exc_upper")
sn_exc <- sn[
  analysis == "ALS_focused" &
    contrast == "ALS_vs_Ctrl_mean" &
    cell_type %chin% exc_types,
  .(gene, cell_type, sn_log2FC = log2FC, sn_PValue = PValue, sn_FDR = FDR)
]
sn_effect_wide <- dcast(sn_exc, gene ~ cell_type, value.var = "sn_log2FC")
setnames(
  sn_effect_wide,
  exc_types,
  paste0("sn_", c("exc_deep", "exc_intermediate", "exc_upper"), "_log2FC")
)

published_exc <- published[
  disease == "ALS" & cell_type %chin% exc_types
]
published_exc[, published_significant :=
  FDR < 0.05 &
    abs(model_log2FC) > log2(1.2) &
    model_log2FC_ci_low * model_log2FC_ci_hi > 0 &
    conv_C == TRUE & conv_D == TRUE
]

candidates <- merge(bulk_mean, sn_effect_wide, by = "gene", all.x = TRUE)
sn_columns <- paste0("sn_", c("exc_deep", "exc_intermediate", "exc_upper"), "_log2FC")
candidates[, sn_exc_available := rowSums(!is.na(.SD)), .SDcols = sn_columns]
candidates[, sn_exc_same_direction := rowSums(
  sweep(as.matrix(.SD), 1L, sign(bulk_log2FC), FUN = function(x, s) sign(x) == s),
  na.rm = TRUE
), .SDcols = sn_columns]
candidates[, sn_exc_median_log2FC := apply(as.matrix(.SD), 1L, median, na.rm = TRUE), .SDcols = sn_columns]
candidates[!is.finite(sn_exc_median_log2FC), sn_exc_median_log2FC := NA_real_]

published_summary <- merge(
  published_exc,
  candidates[, .(gene, bulk_sign = sign(bulk_log2FC))],
  by = "gene",
  all = FALSE
)[, .(
  published_exc_tests = .N,
  published_exc_same_direction = sum(sign(model_log2FC) == bulk_sign, na.rm = TRUE),
  published_exc_significant = sum(published_significant, na.rm = TRUE),
  published_exc_significant_same_direction = sum(
    published_significant & sign(model_log2FC) == bulk_sign,
    na.rm = TRUE
  ),
  published_exc_median_log2FC = median(model_log2FC, na.rm = TRUE)
), by = gene]
candidates <- merge(candidates, published_summary, by = "gene", all.x = TRUE)

count_columns <- c(
  "published_exc_tests",
  "published_exc_same_direction",
  "published_exc_significant",
  "published_exc_significant_same_direction"
)
for (column in count_columns) set(candidates, which(is.na(candidates[[column]])), column, 0L)

panel_annotation <- unique(panel[, .(
  gene = gene_symbol,
  panelapp_confidence = confidence_level,
  panelapp_inheritance = mode_of_inheritance,
  panelapp_phenotypes = phenotypes
)], by = "gene")
candidates <- merge(candidates, panel_annotation, by = "gene", all.x = TRUE)
candidates[, panelapp_ALS_gene := !is.na(panelapp_confidence)]
candidates[, panelapp_green := panelapp_confidence == "3"]
candidates[is.na(panelapp_green), panelapp_green := FALSE]

tdp_annotation <- unique(tdp[, .(
  gene,
  tdp43_cryptic_exon_prior_target = TRUE,
  cryptic_exon_detected_in_Gittings_snRNA = ce_detected_in_gittings_snrna
)], by = "gene")
candidates <- merge(candidates, tdp_annotation, by = "gene", all.x = TRUE)
candidates[is.na(tdp43_cryptic_exon_prior_target), tdp43_cryptic_exon_prior_target := FALSE]

candidates[, evidence_class := fcase(
  bulk_treat_FDR < 0.05 & panelapp_green & sn_exc_same_direction >= 2L,
    "genetic_anchor_bulk_and_sn",
  bulk_treat_FDR < 0.05 & sn_exc_same_direction == 3L &
    published_exc_significant_same_direction >= 1L,
    "cross_platform_bulk_sn_published",
  bulk_treat_FDR < 0.05 & sn_exc_same_direction >= 2L,
    "bulk_plus_sn_direction",
  bulk_treat_FDR < 0.05,
    "bulk_TREAT_only",
  bulk_FDR < 0.05 & sn_exc_same_direction >= 2L,
    "bulk_FDR_plus_sn_direction",
  default = "other_tested_gene"
)]

class_order <- c(
  "genetic_anchor_bulk_and_sn",
  "cross_platform_bulk_sn_published",
  "bulk_plus_sn_direction",
  "bulk_TREAT_only",
  "bulk_FDR_plus_sn_direction",
  "other_tested_gene"
)
candidates[, evidence_class := factor(evidence_class, levels = class_order)]
setorder(candidates, evidence_class, bulk_treat_FDR, bulk_FDR)
candidates[, evidence_class := as.character(evidence_class)]

shortlist <- candidates[evidence_class %chin% class_order[1:3]]
summary <- candidates[, .(
  genes = .N,
  median_abs_bulk_log2FC = median(abs(bulk_log2FC)),
  panelapp_green_genes = sum(panelapp_green),
  prior_TDP43_targets = sum(tdp43_cryptic_exon_prior_target)
), by = evidence_class]

fwrite(candidates, file.path(out_dir, "ALS_protein_coding_multievidence_table.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(shortlist, file.path(out_dir, "ALS_multievidence_shortlist.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(summary, file.path(out_dir, "ALS_multievidence_class_summary.tsv"), sep = "\t", quote = FALSE, na = "")

print(summary)
print(shortlist[1:min(.N, 30L), .(
  gene, evidence_class, bulk_log2FC, bulk_treat_FDR,
  sn_exc_same_direction, published_exc_significant_same_direction,
  panelapp_green, tdp43_cryptic_exon_prior_target
)])
