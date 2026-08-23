#!/usr/bin/env Rscript

# Cross-cohort evidence synthesis for the pre-specified seven mechanism
# families and 17 candidate genes. This script deliberately does not pool NES
# with regression coefficients and does not meta-analyse non-identical
# estimands (diagnostic contrasts versus ordinal pTDP burden).

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
project_root <- analysis_project_root()
configure_project_library(project_root)
project_lib <- file.path(project_root, "cache", "R_library")
.libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

out_dir <- file.path(project_root, "results", "submission_cross_cohort_evidence_v1")
fig_dir <- file.path(project_root, "figures", "submission_cross_cohort_evidence_v1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

canonical_family <- function(x) {
  x <- gsub("Synapse/junction", "Synapse / junction assembly", x, fixed = TRUE)
  x
}

canonical_cell <- function(x) {
  fcase(
    x %chin% c("Astro", "ASC"), "Astrocyte",
    x %chin% c("Exc", "EX"), "Excitatory neuron",
    x %chin% c("Inh", "IN"), "Inhibitory neuron",
    x %chin% c("Micro", "MG"), "Microglia",
    x %chin% c("Oligo", "ODC"), "Oligodendrocyte",
    x == "OPC", "OPC",
    x == "Vascular", "Vascular",
    default = x
  )
}

axis_label <- function(x) {
  fcase(
    x %chin% c("ALS_vs_Ctrl_mean", "C9_ALS_vs_Ctrl_Frontal"), "ALS_vs_Control",
    x %chin% c("FTD_vs_Ctrl_mean", "C9_FTD_vs_Ctrl_Frontal"), "FTD_vs_Control",
    x %chin% c("ALS_vs_FTD_mean", "C9_ALS_vs_C9_FTD_Frontal_batch_bridged"), "ALS_vs_FTD",
    default = x
  )
}

read_tsv <- function(...) fread(file.path(project_root, ...), na.strings = c("", "NA"))

main_program <- read_tsv("results", "GSE219280_same_celltype_fixed_programs", "fixed_family_fgsea_results.tsv")
main_program <- main_program[analysis_role == "primary"]
main_loo <- read_tsv("results", "GSE219280_same_celltype_fixed_programs_LOO", "leave_one_donor_out_program_stability.tsv")
main_gene_all <- read_tsv("results", "GSE219280_pseudobulk_DE", "all_gene_contrast_results.tsv")
main_n <- unique(main_gene_all[
  analysis == "all_disease_core_glia" &
    contrast %chin% c("ALS_vs_Ctrl_mean", "FTD_vs_Ctrl_mean", "ALS_vs_FTD_mean"),
  .(cell_type, contrast, n_donors = as.integer(n_donors))
])
main_program <- merge(main_program, main_n, by = c("cell_type", "contrast"), all.x = TRUE)

gittings_program <- read_tsv("results", "Gittings_2023_C9_snRNA", "fixed_7_family_fgsea_results.tsv")
gittings_program <- gittings_program[analysis_role == "primary_frontal"]
gittings_loo <- read_tsv("results", "Gittings_2023_C9_snRNA", "fixed_7_family_LOO_stability_summary.tsv")
gittings_loo <- gittings_loo[contrast %chin% unique(gittings_program$contrast)]

multi_program <- read_tsv("results", "GSE212630_multiome", "GSE212630_RNA_7_family_donor_trend.tsv")

main_program <- merge(
  main_program,
  main_loo[, .(cell_type, contrast, set_id,
               loo_direction_fraction = direction_retained_fraction,
               loo_iterations = leaveouts_tested)],
  by = c("cell_type", "contrast", "set_id"), all.x = TRUE
)
gittings_program <- merge(
  gittings_program,
  gittings_loo[, .(cell_type, contrast, set_id,
                   loo_direction_fraction = direction_concordant_fraction,
                   loo_iterations = valid_iterations)],
  by = c("cell_type", "contrast", "set_id"), all.x = TRUE
)

program_main <- main_program[, .(
  cohort = "GSE219280", cohort_role = "discovery", genetic_background = "C9orf72",
  region = "BA4/BA9 model-averaged", estimand = axis_label(contrast),
  cell_type = canonical_cell(cell_type), family_id = set_id,
  family_label = canonical_family(set_label), effect_scale = "FGSEA_NES",
  effect = NES, p_value = pval, adjusted_p = global_primary_FDR,
  adjusted_p_scope = "global primary family tests", n_donors = as.integer(n_donors),
  loo_direction_fraction, loo_iterations,
  source_file = "results/GSE219280_same_celltype_fixed_programs/fixed_family_fgsea_results.tsv"
)]

program_gittings <- gittings_program[, .(
  cohort = "Gittings_2023", cohort_role = "independent_C9_replication",
  genetic_background = "C9orf72", region = "frontal cortex",
  estimand = axis_label(contrast), cell_type = canonical_cell(cell_type),
  family_id = set_id, family_label = canonical_family(set_label),
  effect_scale = "FGSEA_NES", effect = NES, p_value = pval,
  adjusted_p = global_primary_FDR, adjusted_p_scope = "global primary family tests",
  n_donors = as.integer(n_donors), loo_direction_fraction, loo_iterations,
  source_file = "results/Gittings_2023_C9_snRNA/fixed_7_family_fgsea_results.tsv"
)]

program_multi <- multi_program[, .(
  cohort = "GSE212630_Emory", cohort_role = "paired_multiome_pathology_context",
  genetic_background = "Control/C9orf72 pathology spectrum", region = "BA9",
  estimand = "pTDP_ordinal_trend", cell_type = canonical_cell(cell_type),
  family_id = gsub("[^A-Za-z0-9]+", "_", tolower(family)),
  family_label = canonical_family(family), effect_scale = "RNA_logCPM_beta_per_pTDP_level",
  effect = beta_per_pTDP_level, p_value, adjusted_p = global_FDR,
  adjusted_p_scope = "global 42 cell-type-by-family trends", n_donors = as.integer(n_donors),
  loo_direction_fraction = NA_real_, loo_iterations = NA_integer_,
  source_file = "results/GSE212630_multiome/GSE212630_RNA_7_family_donor_trend.tsv"
)]

program_evidence <- rbindlist(list(program_main, program_gittings, program_multi), fill = TRUE)
program_evidence[, direction := fifelse(effect > 0, "positive", fifelse(effect < 0, "negative", "zero"))]
program_evidence[, significant_global_FDR_0_05 := !is.na(adjusted_p) & adjusted_p < 0.05]
program_evidence[, inference_tier := fcase(
  cohort == "GSE219280", "discovery",
  cohort == "Gittings_2023", "independent replication; no public donor overlap detected",
  cohort == "GSE212630_Emory", "paired-multiome pathology context; different estimand",
  default = "supportive"
)]
setorder(program_evidence, estimand, cell_type, family_label, cohort)
fwrite(program_evidence, file.path(out_dir, "unified_7_family_evidence.tsv"), sep = "\t", quote = FALSE)

# Exact-estimand direction replication is restricted to the two diagnostic
# cohorts, matched by phenotype contrast, canonical cell type and family.
m <- program_evidence[cohort == "GSE219280"]
g <- program_evidence[cohort == "Gittings_2023"]
program_replication <- merge(
  m, g, by = c("estimand", "cell_type", "family_id", "family_label"),
  suffixes = c("_discovery", "_replication")
)
program_replication[, same_direction := sign(effect_discovery) == sign(effect_replication)]
program_replication[, both_global_FDR_0_05 := significant_global_FDR_0_05_discovery & significant_global_FDR_0_05_replication]
program_replication[, replication_class := fcase(
  both_global_FDR_0_05 & same_direction, "bidirectional_global_FDR_concordant",
  significant_global_FDR_0_05_discovery & same_direction, "discovery_significant_direction_supported",
  significant_global_FDR_0_05_replication & same_direction, "replication_significant_direction_supported",
  same_direction, "direction_only",
  default = "discordant_direction"
)]
keep_rep <- c(
  "estimand", "cell_type", "family_id", "family_label",
  "effect_discovery", "adjusted_p_discovery", "loo_direction_fraction_discovery",
  "effect_replication", "adjusted_p_replication", "loo_direction_fraction_replication",
  "same_direction", "both_global_FDR_0_05", "replication_class"
)
program_replication <- program_replication[, ..keep_rep]
setorder(program_replication, estimand, cell_type, family_label)
fwrite(program_replication, file.path(out_dir, "exact_estimand_program_replication.tsv"), sep = "\t", quote = FALSE)

program_summary <- program_replication[, .(
  axes = .N,
  same_direction_n = sum(same_direction),
  same_direction_fraction = mean(same_direction),
  discovery_global_FDR_n = sum(adjusted_p_discovery < 0.05, na.rm = TRUE),
  replication_global_FDR_n = sum(adjusted_p_replication < 0.05, na.rm = TRUE),
  both_global_FDR_same_direction_n = sum(both_global_FDR_0_05 & same_direction)
), by = .(estimand, cell_type)]
fwrite(program_summary, file.path(out_dir, "exact_estimand_program_replication_summary.tsv"), sep = "\t", quote = FALSE)

# Gene-level synthesis. Log2 fold changes can be displayed with model-based
# confidence intervals, but no pooled estimate is produced because the cohort
# models and sampling frames differ.
candidate_genes <- c(
  "BAG3", "HSP90AA1", "HSP90AB1", "CHORDC1", "ARHGAP35", "PARK7",
  "CCT5", "ARHGAP39", "COX5B", "COX7C", "MAP1B", "MAPT", "BCR",
  "CCT8", "COA3", "DOCK7", "STARD13"
)

main_gene <- copy(main_gene_all)
main_gene <- main_gene[
  analysis == "all_disease_core_glia" & cell_type %chin% c("Astro", "OPC", "Oligo") &
    contrast %chin% c("ALS_vs_Ctrl_mean", "FTD_vs_Ctrl_mean", "ALS_vs_FTD_mean") &
    gene %chin% candidate_genes
]
main_gene[, se := abs(log2FC / t)]
main_gene[!is.finite(se), se := NA_real_]
gene_main <- main_gene[, .(
  cohort = "GSE219280", cohort_role = "discovery", genetic_background = "C9orf72",
  region = "BA4/BA9 model-averaged", estimand = axis_label(contrast),
  cell_type = canonical_cell(cell_type), gene, effect_scale = "log2FC",
  effect = log2FC, se, ci_low = log2FC - 1.96 * se, ci_high = log2FC + 1.96 * se,
  p_value = PValue, adjusted_p = global_FDR,
  adjusted_p_scope = "global gene tests within model collection", n_donors = as.integer(n_donors),
  source_file = "results/GSE219280_pseudobulk_DE/all_gene_contrast_results.tsv"
)]

gittings_gene <- read_tsv("results", "Gittings_2023_C9_snRNA", "predefined_17_gene_results.tsv")
gittings_gene <- gittings_gene[analysis_role == "primary_frontal" & tested == TRUE]
gittings_gene[, se := abs(log2FC / t)]
gittings_gene[!is.finite(se), se := NA_real_]
gene_gittings <- gittings_gene[, .(
  cohort = "Gittings_2023", cohort_role = "independent_C9_replication", genetic_background = "C9orf72",
  region = "frontal cortex", estimand = axis_label(contrast),
  cell_type = canonical_cell(cell_type), gene, effect_scale = "log2FC",
  effect = log2FC, se, ci_low = log2FC - 1.96 * se, ci_high = log2FC + 1.96 * se,
  p_value = PValue, adjusted_p = FDR, adjusted_p_scope = "whole-transcriptome FDR per contrast/cell type",
  n_donors = as.integer(n_donors),
  source_file = "results/Gittings_2023_C9_snRNA/predefined_17_gene_results.tsv"
)]

multi_gene <- read_tsv("results", "GSE212630_multiome", "GSE212630_predefined_17_gene_RNA_results.tsv")
multi_gene <- multi_gene[cell_type != "whole_tissue"]
multi_gene[, se := NA_real_]
multi_gene[is.finite(F) & F > 0,
           se := abs(log2FC_per_pTDP_level) / sqrt(F)]
gene_multi <- multi_gene[, .(
  cohort = "GSE212630_Emory", cohort_role = "paired_multiome_pathology_context",
  genetic_background = "Control/C9orf72 pathology spectrum", region = "BA9",
  estimand = "pTDP_ordinal_trend", cell_type = canonical_cell(cell_type), gene = feature,
  effect_scale = "log2FC_per_pTDP_level", effect = log2FC_per_pTDP_level,
  se, ci_low = log2FC_per_pTDP_level - 1.96 * se,
  ci_high = log2FC_per_pTDP_level + 1.96 * se,
  p_value, adjusted_p = target_global_FDR,
  adjusted_p_scope = "global pre-specified target-by-cell-type tests", n_donors = as.integer(n_donors),
  source_file = "results/GSE212630_multiome/GSE212630_predefined_17_gene_RNA_results.tsv"
)]

gene_evidence <- rbindlist(list(gene_main, gene_gittings, gene_multi), fill = TRUE)
gene_evidence[, direction := fifelse(effect > 0, "positive", fifelse(effect < 0, "negative", "zero"))]
gene_evidence[, significant_adjusted_0_05 := !is.na(adjusted_p) & adjusted_p < 0.05]
gene_evidence[, meta_eligible := FALSE]
gene_evidence[, meta_exclusion_reason := fcase(
  cohort == "GSE212630_Emory", "ordinal pathology trend is not a diagnostic contrast",
  default = "displayable common log2FC scale, but cohort/model/region heterogeneity precludes an unqualified pooled estimate"
)]
setorder(gene_evidence, estimand, cell_type, gene, cohort)
fwrite(gene_evidence, file.path(out_dir, "unified_17_gene_evidence_forest_data.tsv"), sep = "\t", quote = FALSE)

gene_diagnostic <- gene_evidence[cohort %chin% c("GSE219280", "Gittings_2023")]
gene_replication <- merge(
  gene_diagnostic[cohort == "GSE219280"],
  gene_diagnostic[cohort == "Gittings_2023"],
  by = c("estimand", "cell_type", "gene"), suffixes = c("_discovery", "_replication")
)
gene_replication[, same_direction := sign(effect_discovery) == sign(effect_replication)]
gene_replication[, replication_class := fcase(
  adjusted_p_discovery < 0.05 & adjusted_p_replication < 0.05 & same_direction,
  "both_adjusted_significant_concordant",
  adjusted_p_discovery < 0.05 & same_direction, "discovery_adjusted_significant_direction_supported",
  adjusted_p_replication < 0.05 & same_direction, "replication_adjusted_significant_direction_supported",
  same_direction, "direction_only", default = "discordant_direction"
)]
gene_replication <- gene_replication[, .(
  estimand, cell_type, gene,
  discovery_log2FC = effect_discovery, discovery_ci_low = ci_low_discovery,
  discovery_ci_high = ci_high_discovery, discovery_adjusted_p = adjusted_p_discovery,
  replication_log2FC = effect_replication, replication_ci_low = ci_low_replication,
  replication_ci_high = ci_high_replication, replication_adjusted_p = adjusted_p_replication,
  same_direction, replication_class
)]
setorder(gene_replication, estimand, cell_type, gene)
fwrite(gene_replication, file.path(out_dir, "exact_estimand_17_gene_replication.tsv"), sep = "\t", quote = FALSE)

# A restrained faceted forest plot: only diagnostic cohort log2FCs are shown;
# no pooled diamond is drawn.
plot_data <- gene_evidence[
  cohort %chin% c("GSE219280", "Gittings_2023") &
    estimand %chin% c("ALS_vs_Control", "FTD_vs_Control", "ALS_vs_FTD") &
    is.finite(ci_low) & is.finite(ci_high)
]
priority <- unique(gene_replication[
  replication_class != "discordant_direction" &
    (discovery_adjusted_p < 0.05 | replication_adjusted_p < 0.05), gene
])
if (!length(priority)) priority <- candidate_genes
plot_data <- plot_data[gene %chin% priority]
plot_data[, gene := factor(gene, levels = rev(candidate_genes[candidate_genes %chin% unique(plot_data$gene)]))]

p <- ggplot(plot_data, aes(x = effect, y = gene, colour = cohort)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "grey60") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.16,
                 position = position_dodge(width = 0.45), linewidth = 0.35) +
  geom_point(position = position_dodge(width = 0.45), size = 1.5) +
  facet_grid(cell_type ~ estimand, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(GSE219280 = "#2B6CB0", Gittings_2023 = "#C05621")) +
  labs(x = "Model-based log2 fold change (95% CI; no pooled estimate)", y = NULL,
       colour = "Cohort") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        strip.text = element_text(size = 7), axis.text.y = element_text(size = 6))

ggsave(file.path(fig_dir, "diagnostic_17_gene_cross_cohort_forest.pdf"), p,
       width = 11, height = 9, units = "in", device = cairo_pdf)
ggsave(file.path(fig_dir, "diagnostic_17_gene_cross_cohort_forest.tiff"), p,
       width = 11, height = 9, units = "in", dpi = 400, compression = "lzw")

field_audit <- data.table(
  source = c("GSE219280 fixed-family FGSEA", "Gittings 2023 fixed-family FGSEA",
             "GSE212630 family trend", "GSE219280 17-gene limma",
             "Gittings 2023 17-gene limma", "GSE212630 17-gene edgeR trend"),
  unit = c("donor-blocked pseudobulk rank", "donor-blocked pseudobulk rank",
           "donor-by-cell-type pseudobulk", "donor-blocked pseudobulk",
           "donor-blocked pseudobulk", "donor-by-cell-type pseudobulk"),
  effect_scale = c("NES", "NES", "beta per ordinal pTDP level", "log2FC", "log2FC",
                   "log2FC per ordinal pTDP level"),
  direct_pooling_allowed = FALSE,
  valid_synthesis = c(
    "same-estimand/cell/family direction replication with Gittings",
    "same-estimand/cell/family direction replication with GSE219280",
    "separate pathology-context panel only",
    "side-by-side forest display and direction replication; no pooled diamond",
    "side-by-side forest display and direction replication; no pooled diamond",
    "separate pathology-context panel only"
  ),
  main_boundary = c(
    "BA4/BA9 model-averaged C9 diagnostic effect",
    "frontal C9 diagnostic effect; ALS-vs-FTD is batch-bridged difference-in-differences",
    "14 Emory donors; pTDP burden is not ALS/FTD diagnosis",
    "three glial cell types only in exact family framework",
    "cell-type-specific frontal estimates",
    "same-nucleus multiome but RNA association does not prove ATAC-to-RNA causality"
  )
)
fwrite(field_audit, file.path(out_dir, "source_scale_and_estimand_audit.tsv"), sep = "\t", quote = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat("Wrote", nrow(program_evidence), "program evidence rows\n")
cat("Wrote", nrow(program_replication), "exact-estimand program replication rows\n")
cat("Wrote", nrow(gene_evidence), "gene evidence rows\n")
cat("Wrote", nrow(gene_replication), "exact-estimand gene replication rows\n")
