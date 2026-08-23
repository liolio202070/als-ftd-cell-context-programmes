#!/usr/bin/env Rscript

# Evidence-controlled leading-edge driver prioritization for the C9orf72
# ALS-FTD spectrum project. This script does not perform a new exploratory
# enrichment. It restricts all inference to the 11 predeclared GO terms grouped
# into seven mechanism families and integrates already-completed analyses.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
project_root <- analysis_project_root()
configure_project_library(project_root)
project_lib <- file.path(project_root, "cache", "R_library")
.libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages(library(data.table))

results_root <- file.path(project_root, "results")
out_dir <- file.path(results_root, "leading_edge_evidence_v2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source_files <- c(
  original_leading_edge = file.path(results_root, "GO_leading_edge_sensitivity", "primary_GO_leading_edge_membership.tsv"),
  original_sensitivity = file.path(results_root, "GO_leading_edge_sensitivity", "GO_BP_priority_pathway_sensitivity.tsv"),
  same_cell_family = file.path(results_root, "GSE219280_same_celltype_fixed_programs", "fixed_family_fgsea_results.tsv"),
  same_cell_sensitivity = file.path(results_root, "GSE219280_same_celltype_fixed_programs", "threshold_sensitivity_by_program.tsv"),
  same_cell_loo = file.path(results_root, "GSE219280_same_celltype_fixed_programs_LOO", "leave_one_donor_out_program_stability.tsv"),
  same_cell_de = file.path(results_root, "GSE219280_pseudobulk_DE", "all_gene_contrast_results.tsv"),
  atac_gene_activity = file.path(results_root, "snATAC_GSE219279", "donor_level_gene_activity_priority20", "candidate_gene_activity_results.tsv"),
  gse330_pathways = file.path(results_root, "GSE330_orthogonal_validation", "pathway_validation", "predefined_story_module_GSE330_evidence.tsv"),
  gse330_overlap = file.path(results_root, "GSE330_orthogonal_validation", "donor_overlap_audit", "GSE330_vs_GSE153960_donor_overlap_summary.tsv"),
  proteomics = file.path(results_root, "PXD037393_synaptic_proteomics", "PXD037393_group_control_protein_ratios.tsv")
)
missing_files <- source_files[!file.exists(source_files)]
if (length(missing_files)) stop("Missing required source file(s): ", paste(missing_files, collapse = "; "))

term_definitions <- data.table(
  go_id = c(
    "GO:0006119", "GO:0042775", "GO:0006457", "GO:0009408",
    "GO:0050821", "GO:0051056", "GO:0007264", "GO:0007416",
    "GO:0034329", "GO:0007409", "GO:0061564"
  ),
  term_label = c(
    "Oxidative phosphorylation", "Mitochondrial ATP synthesis",
    "Protein folding", "Heat response", "Protein stabilization",
    "Small GTPase regulation", "Small GTPase signaling",
    "Synapse assembly", "Cell junction assembly", "Axonogenesis",
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
    "Small GTPase", "Synapse / junction assembly",
    "Synapse / junction assembly", "Axon development", "Axon development"
  )
)

original_comparisons <- c(
  "ALS_bulk", "ALS_Exc_deep", "ALS_Exc_intermediate", "ALS_Exc_upper",
  "FTD_Astro", "FTD_OPC", "FTD_Oligo"
)
external_ftd_comparisons <- c("Independent_C9_FTD", "Independent_all_FTD")

expand_leading_edge <- function(dt, column, separator, group_columns) {
  if (!nrow(dt)) return(data.table())
  x <- copy(dt)
  x[, gene_list__ := strsplit(get(column), separator, fixed = TRUE)]
  x <- x[, .(gene = unlist(gene_list__)), by = group_columns]
  x[!is.na(gene) & nzchar(gene)]
}

comparison_id_from_spec <- function(dataset, cell_type, contrast) {
  fcase(
    dataset == "GSE153960", "ALS_bulk",
    dataset == "GSE219280" & contrast == "ALS_vs_Ctrl_mean", paste0("ALS_", cell_type),
    dataset == "GSE219280" & contrast == "FTD_vs_Ctrl_mean", paste0("FTD_", cell_type),
    dataset == "GSE195872" & contrast == "C9orf72_vs_Control", "Independent_C9_FTD",
    dataset == "GSE195872" & contrast == "FTD_vs_Control_all", "Independent_all_FTD",
    default = NA_character_
  )
}

# -------------------------------------------------------------------------
# 1. Original seven comparisons: only fixed terms with comparison-level FDR.
# -------------------------------------------------------------------------
leading <- fread(source_files[["original_leading_edge"]])
leading <- merge(leading, term_definitions, by.x = "ID", by.y = "go_id", all.x = TRUE)
original_membership <- unique(leading[
  !is.na(family_id) & comparison_id %chin% original_comparisons & p.adjust < 0.05,
  .(gene, go_id = ID, term_label, family_id, family_label, comparison_id, NES, pathway_FDR = p.adjust)
])
external_ftd_membership <- unique(leading[
  !is.na(family_id) & comparison_id %chin% external_ftd_comparisons & p.adjust < 0.05,
  .(gene, go_id = ID, term_label, family_id, family_label, comparison_id, NES, pathway_FDR = p.adjust)
])

original_summary <- original_membership[, .(
  original7_comparisons = uniqueN(comparison_id),
  original7_terms = uniqueN(go_id),
  original7_families = uniqueN(family_id),
  original7_ALS_comparisons = uniqueN(comparison_id[grepl("^ALS", comparison_id)]),
  original7_FTD_comparisons = uniqueN(comparison_id[grepl("^FTD", comparison_id)]),
  ALS_bulk_leading_edge = any(comparison_id == "ALS_bulk"),
  original7_comparison_labels = paste(sort(unique(comparison_id)), collapse = ";"),
  original7_family_labels = paste(sort(unique(family_label)), collapse = ";")
), by = gene]

# Gene membership robustness under four alternative GSEA settings.
original_sensitivity <- fread(source_files[["original_sensitivity"]])
original_sensitivity[, comparison_id := comparison_id_from_spec(dataset, cell_type, contrast)]
original_sensitivity <- original_sensitivity[
  pathway %chin% term_definitions$go_id & comparison_id %chin% original_comparisons
]
original_sensitivity_genes <- expand_leading_edge(
  original_sensitivity,
  "leadingEdge", "/",
  c("comparison_id", "pathway", "setting", "same_direction", "sensitivity_significant")
)
original_pairs <- unique(original_membership[, .(gene, comparison_id, pathway = go_id)])
original_setting_hits <- original_sensitivity_genes[
  original_pairs,
  on = .(gene, comparison_id, pathway),
  .(gene = i.gene, comparison_id = i.comparison_id, pathway = i.pathway,
    settings_with_gene = uniqueN(setting)),
  by = .EACHI
]
original_setting_hits[, membership_fraction := settings_with_gene / 4]
original_gene_sensitivity <- original_setting_hits[, .(
  original_LE_sensitivity_mean = mean(membership_fraction),
  original_LE_sensitivity_median = median(membership_fraction),
  original_LE_sensitivity_minimum = min(membership_fraction)
), by = gene]

# -------------------------------------------------------------------------
# 2. Same-cell-type comparisons: global-FDR families, threshold and LOO QC.
# -------------------------------------------------------------------------
same_cell_family <- fread(source_files[["same_cell_family"]])
same_cell_primary <- same_cell_family[
  threshold == "min30_primary" & global_primary_FDR < 0.05
]
same_cell_membership <- expand_leading_edge(
  same_cell_primary,
  "leadingEdge", ";",
  c("cell_type", "contrast", "set_id", "set_label", "NES", "global_primary_FDR")
)
setnames(same_cell_membership, c("set_id", "set_label"), c("family_id", "family_label"))
same_cell_membership <- unique(same_cell_membership)

threshold_results <- fread(source_files[["same_cell_sensitivity"]])
threshold_direction <- threshold_results[, .(
  threshold_directions_retained = all(direction_retained),
  threshold_max_abs_delta_NES = max(absolute_NES_difference)
), by = .(cell_type, contrast, family_id = set_id)]
threshold_gene_membership <- expand_leading_edge(
  threshold_results,
  "leadingEdge", ";",
  c("cell_type", "contrast", "set_id", "threshold")
)
setnames(threshold_gene_membership, "set_id", "family_id")
threshold_gene_hits <- threshold_gene_membership[
  same_cell_membership[, .(gene, cell_type, contrast, family_id)],
  on = .(gene, cell_type, contrast, family_id),
  .(gene = i.gene, cell_type = i.cell_type, contrast = i.contrast,
    family_id = i.family_id, thresholds_with_gene = uniqueN(threshold)),
  by = .EACHI
]
threshold_gene_hits[, membership_fraction := thresholds_with_gene / 2]

loo <- fread(source_files[["same_cell_loo"]])
loo <- loo[, .(
  cell_type, contrast, family_id = set_id,
  loo_direction_retained_fraction = direction_retained_fraction,
  loo_nominal_fraction = nominal_P_0_05_fraction
)]
same_cell_membership <- merge(
  same_cell_membership, threshold_direction,
  by = c("cell_type", "contrast", "family_id"), all.x = TRUE
)
same_cell_membership <- merge(
  same_cell_membership, loo,
  by = c("cell_type", "contrast", "family_id"), all.x = TRUE
)

same_cell_summary <- same_cell_membership[, .(
  same_cell_significant_comparisons = uniqueN(paste(cell_type, contrast)),
  same_cell_significant_families = uniqueN(family_id),
  same_cell_ALS_celltypes = uniqueN(cell_type[contrast == "ALS_vs_Ctrl_mean"]),
  same_cell_FTD_celltypes = uniqueN(cell_type[contrast == "FTD_vs_Ctrl_mean"]),
  direct_ALS_vs_FTD_celltypes = uniqueN(cell_type[contrast == "ALS_vs_FTD_mean"]),
  same_cell_program_stable = all(
    threshold_directions_retained & loo_direction_retained_fraction >= 0.94
  ),
  minimum_LOO_direction_fraction = min(loo_direction_retained_fraction),
  same_cell_family_labels = paste(sort(unique(family_label)), collapse = ";")
), by = gene]
same_cell_gene_sensitivity <- threshold_gene_hits[, .(
  same_cell_LE_threshold_mean = mean(membership_fraction),
  same_cell_LE_threshold_median = median(membership_fraction),
  same_cell_LE_threshold_minimum = min(membership_fraction)
), by = gene]

# -------------------------------------------------------------------------
# 3. External/orthogonal and cross-modal evidence.
# -------------------------------------------------------------------------
external_ftd_summary <- external_ftd_membership[, .(
  external_FTD_comparisons = uniqueN(comparison_id),
  external_FTD_terms = uniqueN(go_id),
  external_FTD_family_labels = paste(sort(unique(family_label)), collapse = ";")
), by = gene]

gse330 <- fread(source_files[["gse330_pathways"]])[FDR < 0.05]
gse330 <- merge(gse330, term_definitions[, .(go_id, family_id, family_label)], by.x = "ID", by.y = "go_id")
gse330_membership <- expand_leading_edge(
  gse330,
  "leadingEdge", "/",
  c("dataset", "unit", "ID", "module", "family_id", "family_label", "NES", "FDR")
)
gse330_summary <- gse330_membership[, .(
  GSE330_significant_units = uniqueN(paste(dataset, unit)),
  GSE330_datasets = uniqueN(dataset),
  GSE330_families = uniqueN(family_id),
  GSE330_family_labels = paste(sort(unique(family_label)), collapse = ";")
), by = gene]

sn <- fread(source_files[["same_cell_de"]])
atac <- fread(source_files[["atac_gene_activity"]])[
  minimum_cells == 50 & model == "diagnosis_only" &
    role %chin% c("primary_axis", "supportive_context") & global_primary_FDR < 0.05
]
atac_map <- data.table(
  comparison = c(
    "ALS_motor_Exc_deep", "ALS_motor_Exc_intermediate", "ALS_motor_Exc_upper",
    "FTD_frontal_Astro", "FTD_frontal_Oligo", "ALS_motor_Astro_context"
  ),
  analysis = c(
    rep("ALS_focused", 3), rep("all_disease_core_glia", 3)
  ),
  cell_type = c("Exc_deep", "Exc_intermediate", "Exc_upper", "Astro", "Oligo", "Astro"),
  contrast = c(
    rep("ALS_vs_Ctrl_mean", 3), rep("FTD_vs_Ctrl_mean", 2), "ALS_vs_Ctrl_mean"
  )
)
atac <- merge(atac, atac_map, by = "comparison", all.x = TRUE)
rna_effects <- sn[, .(analysis, cell_type, contrast, gene, RNA_log2FC = log2FC, RNA_FDR = FDR)]
atac <- merge(atac, rna_effects, by = c("analysis", "cell_type", "contrast", "gene"), all.x = TRUE)
atac[, RNA_ATAC_direction_concordant := sign(effect) == sign(RNA_log2FC)]
atac_summary <- atac[, .(
  same_donor_ATAC_significant_axes = uniqueN(comparison),
  same_donor_ATAC_concordant_axes = uniqueN(comparison[RNA_ATAC_direction_concordant == TRUE]),
  same_donor_ATAC_ALS_axes = uniqueN(comparison[grepl("^ALS", comparison)]),
  same_donor_ATAC_FTD_axes = uniqueN(comparison[grepl("^FTD", comparison)]),
  minimum_ATAC_global_FDR = min(global_primary_FDR),
  ATAC_axis_effects = paste(
    paste0(comparison, ":", sprintf("%.3f", effect),
           "/RNA:", sprintf("%.3f", RNA_log2FC)),
    collapse = ";"
  )
), by = gene]

protein <- fread(source_files[["proteomics"]])
protein[, protein_max_abs_log2ratio := pmax(
  abs(ALS_BA9_log2ratio), abs(ALS_BA4_log2ratio), abs(C9pos_BA9_log2ratio),
  na.rm = TRUE
)]
protein[, proteomics_20pct_descriptive := protein_max_abs_log2ratio >= log2(1.2)]
protein <- protein[, .(
  gene, proteomics_detected = TRUE, proteomics_20pct_descriptive,
  protein_max_abs_log2ratio, ALS_BA9_log2ratio, ALS_BA4_log2ratio,
  C9pos_BA9_log2ratio
)]

# -------------------------------------------------------------------------
# 4. Assemble evidence, apply explicit gates, and score.
# -------------------------------------------------------------------------
evidence <- Reduce(
  function(x, y) merge(x, y, by = "gene", all = TRUE),
  list(
    original_summary, original_gene_sensitivity, same_cell_summary,
    same_cell_gene_sensitivity, external_ftd_summary, gse330_summary,
    atac_summary, protein
  )
)

integer_zero <- c(
  "original7_comparisons", "original7_terms", "original7_families",
  "original7_ALS_comparisons", "original7_FTD_comparisons",
  "same_cell_significant_comparisons", "same_cell_significant_families",
  "same_cell_ALS_celltypes", "same_cell_FTD_celltypes",
  "direct_ALS_vs_FTD_celltypes", "external_FTD_comparisons",
  "external_FTD_terms", "GSE330_significant_units", "GSE330_datasets",
  "GSE330_families", "same_donor_ATAC_significant_axes",
  "same_donor_ATAC_concordant_axes", "same_donor_ATAC_ALS_axes",
  "same_donor_ATAC_FTD_axes"
)
for (column in intersect(integer_zero, names(evidence))) {
  set(evidence, which(is.na(evidence[[column]])), column, 0L)
}
logical_false <- c(
  "ALS_bulk_leading_edge", "same_cell_program_stable",
  "proteomics_detected", "proteomics_20pct_descriptive"
)
for (column in intersect(logical_false, names(evidence))) {
  set(evidence, which(is.na(evidence[[column]])), column, FALSE)
}
for (column in c(
  "original_LE_sensitivity_mean", "original_LE_sensitivity_median",
  "original_LE_sensitivity_minimum", "same_cell_LE_threshold_mean",
  "same_cell_LE_threshold_median", "same_cell_LE_threshold_minimum"
)) {
  set(evidence, which(is.na(evidence[[column]])), column, 0)
}

evidence[, external_or_crossmodal_layers :=
  (external_FTD_comparisons > 0L) + (GSE330_significant_units > 0L) +
  (same_donor_ATAC_concordant_axes > 0L)]

evidence[, shared_gate :=
  original7_comparisons >= 4L & original7_ALS_comparisons >= 2L &
  original7_FTD_comparisons >= 2L & original_LE_sensitivity_median >= 0.5 &
  same_cell_FTD_celltypes >= 2L & same_cell_significant_comparisons >= 2L &
  same_cell_program_stable & same_cell_LE_threshold_median >= 0.5 &
  (external_or_crossmodal_layers >= 1L |
     (direct_ALS_vs_FTD_celltypes >= 2L & same_cell_FTD_celltypes >= 2L))]

evidence[, ALS_specific_gate :=
  ALS_bulk_leading_edge & original7_ALS_comparisons >= 3L &
  original7_FTD_comparisons <= 1L & original_LE_sensitivity_median >= 0.5 &
  (GSE330_significant_units > 0L | same_donor_ATAC_ALS_axes > 0L)]

evidence[, FTD_specific_gate :=
  original7_FTD_comparisons >= 2L & original7_ALS_comparisons <= 1L &
  same_cell_FTD_celltypes >= 2L & direct_ALS_vs_FTD_celltypes >= 2L &
  same_cell_program_stable & same_cell_LE_threshold_median >= 0.5 &
  (external_FTD_comparisons > 0L | same_donor_ATAC_FTD_axes > 0L)]

evidence[, phenotype_differentiation_gate :=
  original7_FTD_comparisons >= 2L & original7_comparisons >= 3L &
  same_cell_FTD_celltypes >= 2L & direct_ALS_vs_FTD_celltypes >= 2L &
  same_cell_program_stable & same_cell_LE_threshold_median >= 0.5 &
  external_or_crossmodal_layers >= 1L]

evidence[, strict_candidate :=
  shared_gate | ALS_specific_gate | FTD_specific_gate | phenotype_differentiation_gate]

evidence[, spectrum_class := fcase(
  ALS_specific_gate & !shared_gate, "ALS-associated",
  FTD_specific_gate & !shared_gate, "FTD-glial-biased",
  shared_gate & direct_ALS_vs_FTD_celltypes >= 2L,
    "Shared program, FTD-amplified",
  shared_gate, "Shared recurrent program",
  phenotype_differentiation_gate, "Phenotype differentiation (FTD higher)",
  original7_ALS_comparisons >= 2L & original7_FTD_comparisons >= 2L,
    "Shared, insufficient stringent support",
  original7_ALS_comparisons > original7_FTD_comparisons,
    "ALS-leaning exploratory",
  original7_FTD_comparisons > original7_ALS_comparisons,
    "FTD-leaning exploratory",
  default = "Unclassified exploratory"
)]

evidence[, evidence_score :=
  pmin(original7_comparisons, 5L) + pmin(original7_families, 3L) +
  as.integer(ALS_bulk_leading_edge) +
  as.integer(original7_ALS_comparisons >= 2L) +
  as.integer(original7_FTD_comparisons >= 2L) +
  pmin(same_cell_significant_comparisons, 4L) +
  2L * as.integer(direct_ALS_vs_FTD_celltypes >= 2L) +
  as.integer(original_LE_sensitivity_median >= 0.5) +
  as.integer(same_cell_LE_threshold_median >= 0.5) +
  as.integer(same_cell_program_stable) +
  2L * as.integer(external_FTD_comparisons > 0L) +
  2L * as.integer(GSE330_significant_units > 0L) +
  2L * as.integer(same_donor_ATAC_concordant_axes > 0L) +
  as.integer(proteomics_20pct_descriptive)
]

evidence[, evidence_tier := fcase(
  strict_candidate & evidence_score >= 20 & external_or_crossmodal_layers >= 2L,
    "Tier 1: convergent multi-layer driver",
  strict_candidate,
    "Tier 2: robust transcriptomic driver",
  default = "Tier 3: exploratory/not shortlisted"
)]

evidence[, gate_failure_reason := fifelse(
  strict_candidate, "",
  paste0(
    fifelse(original7_comparisons < 3L, "original7_recurrence;", ""),
    fifelse(original_LE_sensitivity_median < 0.5, "original_LE_sensitivity;", ""),
    fifelse(same_cell_significant_comparisons < 2L, "same_cell_recurrence;", ""),
    fifelse(!same_cell_program_stable, "LOO_or_threshold_direction;", ""),
    fifelse(external_or_crossmodal_layers < 1L, "orthogonal_or_crossmodal_support;", "")
  )
)]

# -------------------------------------------------------------------------
# 5. Mechanism-balanced shortlist: top three strict genes per family plus the
#    top two representatives of each spectrum class. No gate is relaxed to
#    reach a target count; the final list may therefore contain fewer than 30.
# -------------------------------------------------------------------------
family_membership <- unique(rbindlist(list(
  original_membership[, .(gene, family_id, family_label)],
  same_cell_membership[, .(gene, family_id, family_label)],
  external_ftd_membership[, .(gene, family_id, family_label)],
  gse330_membership[, .(gene, family_id, family_label)]
)))
family_counts <- rbindlist(list(
  original_membership[, .(original_family_comparisons = uniqueN(comparison_id)),
                      by = .(gene, family_id, family_label)],
  same_cell_membership[, .(same_cell_family_comparisons = uniqueN(paste(cell_type, contrast))),
                         by = .(gene, family_id, family_label)],
  external_ftd_membership[, .(external_FTD_family_comparisons = uniqueN(comparison_id)),
                          by = .(gene, family_id, family_label)],
  gse330_membership[, .(GSE330_family_units = uniqueN(paste(dataset, unit))),
                    by = .(gene, family_id, family_label)]
), use.names = TRUE, fill = TRUE)
family_counts <- family_counts[, lapply(.SD, function(z) sum(z, na.rm = TRUE)),
                               by = .(gene, family_id, family_label)]
family_rank <- merge(
  family_membership, family_counts,
  by = c("gene", "family_id", "family_label"), all.x = TRUE
)
family_rank <- merge(
  family_rank,
  evidence[strict_candidate == TRUE, .(gene, evidence_score, evidence_tier, spectrum_class)],
  by = "gene"
)
for (column in c(
  "original_family_comparisons", "same_cell_family_comparisons",
  "external_FTD_family_comparisons", "GSE330_family_units"
)) set(family_rank, which(is.na(family_rank[[column]])), column, 0)
family_rank[, family_specific_score :=
  original_family_comparisons + same_cell_family_comparisons +
  2L * (external_FTD_family_comparisons > 0L) +
  2L * (GSE330_family_units > 0L)]
setorder(family_rank, family_id, -evidence_score, -family_specific_score, gene)
family_top <- family_rank[, head(.SD, 3L), by = family_id]
family_selection <- family_top[, .(
  selected_by = paste0("family_top3:", paste(sort(unique(family_label)), collapse = ","))
), by = gene]

candidate_classes <- c(
  "ALS-associated", "FTD-glial-biased", "Shared program, FTD-amplified",
  "Shared recurrent program", "Phenotype differentiation (FTD higher)"
)
class_top <- evidence[
  strict_candidate == TRUE & spectrum_class %chin% candidate_classes
][order(spectrum_class, -evidence_score, -original7_comparisons, gene), head(.SD, 2L),
  by = spectrum_class]
class_selection <- class_top[, .(
  gene,
  selected_by = paste0("class_top2:", spectrum_class)
)]
selection <- rbindlist(list(family_selection, class_selection), use.names = TRUE)
selection <- selection[, .(selected_by = paste(sort(unique(selected_by)), collapse = ";")), by = gene]

shortlist <- merge(evidence[strict_candidate == TRUE], selection, by = "gene")
setorder(shortlist, -evidence_score, evidence_tier, gene)
if (nrow(shortlist) > 30L) shortlist <- shortlist[1:30]
shortlist[, shortlist_rank := seq_len(.N)]
setcolorder(shortlist, c(
  "shortlist_rank", "gene", "spectrum_class", "evidence_tier",
  "evidence_score", "selected_by", setdiff(names(shortlist), c(
    "shortlist_rank", "gene", "spectrum_class", "evidence_tier",
    "evidence_score", "selected_by"
  ))
))

# Long-form traceability for shortlisted genes.
short_genes <- shortlist$gene
evidence_long <- rbindlist(list(
  original_membership[gene %chin% short_genes, .(
    gene, evidence_layer = "original_7_GSEA", dataset_or_modality = "GSEA",
    comparison = comparison_id, mechanism_family = family_label,
    direction = fifelse(NES > 0, "positive", "negative"),
    statistic = NES, FDR = pathway_FDR,
    interpretation = "Fixed-term significant leading edge"
  )],
  same_cell_membership[gene %chin% short_genes, .(
    gene, evidence_layer = "same_cell_direct_GSEA", dataset_or_modality = "GSE219280_snRNA",
    comparison = paste(cell_type, contrast, sep = ":"), mechanism_family = family_label,
    direction = fifelse(NES > 0, "positive", "negative"),
    statistic = NES, FDR = global_primary_FDR,
    interpretation = "Global-FDR fixed-family leading edge; donor-level"
  )],
  external_ftd_membership[gene %chin% short_genes, .(
    gene, evidence_layer = "external_FTD_GSEA", dataset_or_modality = "GSE195872_bulk_PFC",
    comparison = comparison_id, mechanism_family = family_label,
    direction = fifelse(NES > 0, "positive", "negative"),
    statistic = NES, FDR = pathway_FDR,
    interpretation = "External FTD fixed-term leading edge"
  )],
  gse330_membership[gene %chin% short_genes, .(
    gene, evidence_layer = "orthogonal_ALS_GSE330", dataset_or_modality = dataset,
    comparison = unit, mechanism_family = family_label,
    direction = fifelse(NES > 0, "positive", "negative"),
    statistic = NES, FDR,
    interpretation = "Orthogonal public-data leading edge; modest donor overlap audited"
  )],
  atac[gene %chin% short_genes, .(
    gene, evidence_layer = "same_donor_snATAC", dataset_or_modality = "GSE219279_snATAC_GeneScore",
    comparison, mechanism_family = NA_character_,
    direction = fifelse(effect > 0, "positive", "negative"),
    statistic = effect, FDR = global_primary_FDR,
    interpretation = fifelse(
      RNA_ATAC_direction_concordant,
      "Same-donor cross-modal RNA-ATAC direction concordant",
      "Same-donor chromatin signal; RNA direction discordant"
    )
  )]
), use.names = TRUE, fill = TRUE)

scoring_rules <- data.table(
  rule_id = c(
    "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9",
    "R10", "R11", "R12", "G_shared", "G_ALS", "G_FTD", "G_diff",
    "Tier1", "Tier2", "Shortlist"
  ),
  rule = c(
    "Original seven significant leading-edge comparisons: +1 each, capped at 5",
    "Original seven mechanism-family breadth: +1 each, capped at 3",
    "ALS bulk leading-edge support: +1",
    "At least two ALS original comparisons: +1",
    "At least two FTD original comparisons: +1",
    "Same-cell significant comparison recurrence: +1 each, capped at 4",
    "Direct ALS-vs-FTD leading edge in >=2 glial cell types: +2",
    "Original leading-edge membership retained in >=50% sensitivity settings: +1",
    "Same-cell leading-edge membership retained at median >=50% across min20/min50: +1",
    "All supporting same-cell programs retain direction at min20/min50 and LOO >=0.94: +1",
    "External FTD / GSE330 / concordant same-donor ATAC: +2 per layer",
    "Group-mean synaptic proteomics absolute change >=20%: +1 descriptive only",
    "Shared gate: recurrent ALS+FTD original evidence plus stable same-cell FTD evidence",
    "ALS gate: ALS bulk + >=3 ALS comparisons, <=1 FTD comparison, sensitivity and GSE330/ATAC support",
    "FTD gate: >=2 FTD comparisons, stable FTD-vs-control and ALS-vs-FTD in >=2 glial types, external FTD/ATAC support",
    "Differentiation gate: stable FTD-vs-control and direct ALS-vs-FTD in >=2 glial types plus an external/cross-modal layer",
    "Strict gate, score >=20, and >=2 external/cross-modal layers",
    "Any strict gate not meeting Tier 1",
    "Top 3 strict genes per mechanism family plus top 2 per spectrum class; union capped at 30; gates are never relaxed"
  )
)

source_manifest <- data.table(
  source_id = names(source_files),
  path = unname(source_files),
  md5 = unname(tools::md5sum(source_files)),
  bytes = file.info(source_files)$size
)
gse330_overlap <- fread(source_files[["gse330_overlap"]])
qc <- data.table(
  metric = c(
    "fixed_GO_terms", "mechanism_families", "original_comparisons",
    "original_significant_leading_edge_rows", "same_cell_global_FDR_leading_edge_rows",
    "genes_in_full_evidence_table", "strict_candidates", "tier1_candidates",
    "shortlist_genes", "shortlist_min_score", "shortlist_max_score",
    "shortlist_duplicate_genes", "GSE330_max_overlap_fraction"
  ),
  value = as.character(c(
    nrow(term_definitions), uniqueN(term_definitions$family_id), length(original_comparisons),
    nrow(original_membership), nrow(same_cell_membership), nrow(evidence),
    sum(evidence$strict_candidate), sum(evidence$evidence_tier == "Tier 1: convergent multi-layer driver"),
    nrow(shortlist), min(shortlist$evidence_score), max(shortlist$evidence_score),
    nrow(shortlist) - uniqueN(shortlist$gene), max(gse330_overlap$overlap_fraction)
  ))
)

setorder(evidence, -strict_candidate, -evidence_score, gene)
fwrite(evidence, file.path(out_dir, "all_fixed_program_gene_evidence.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(shortlist, file.path(out_dir, "leading_edge_driver_shortlist.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(evidence_long, file.path(out_dir, "shortlist_evidence_long.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(family_rank, file.path(out_dir, "family_specific_candidate_ranking.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(scoring_rules, file.path(out_dir, "evidence_scoring_and_gating_rules.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(qc, file.path(out_dir, "analysis_QC.tsv"), sep = "\t", quote = FALSE, na = "")
fwrite(source_manifest, file.path(out_dir, "source_file_manifest.tsv"), sep = "\t", quote = FALSE, na = "")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

class_counts <- shortlist[, .N, by = spectrum_class][order(-N)]
tier_counts <- shortlist[, .N, by = evidence_tier][order(evidence_tier)]
summary_lines <- c(
  "# Leading-edge关键驱动基因证据分级（v2）",
  "",
  sprintf("最终短名单包含%d个基因；没有为达到数量目标而放宽任何门槛。", nrow(shortlist)),
  sprintf("其中Tier 1为%d个，Tier 2为%d个。",
          sum(shortlist$evidence_tier == "Tier 1: convergent multi-layer driver"),
          sum(shortlist$evidence_tier == "Tier 2: robust transcriptomic driver")),
  "",
  "## 方法学边界",
  "",
  "- 仅使用事先固定的11个GO条目及其7个机制家族，不重新搜索或扩增通路。",
  "- 原7组GSEA只在比较内FDR<0.05时计入leading-edge证据。",
  "- 同细胞类型结果只在63个预定义主检验的全局FDR<0.05时计入，并要求20/30/50阈值方向及供体留一方向稳定。",
  "- GSE219279 snATAC与snRNA来自同一研究且供体重叠，严格表述为同供体跨模态支持，不作为独立外部验证。",
  "- GSE330存在1–6名供体与GSE153960重叠（最高14.3%），因此作为正交公共数据支持而非完全独立复现。",
  "- PXD037393为组均值蛋白比值，只作为描述性加分，不能替代个体级统计验证。",
  "",
  "## 短名单分层",
  "",
  paste0("- ", class_counts$spectrum_class, "：", class_counts$N, "个"),
  "",
  "## 结论",
  "",
  "严格证据收敛支持以蛋白稳态/伴侣蛋白为共享主轴，并保留线粒体、小GTPase及轴突/突触相关驱动基因。多数共享驱动在FTD胶质中呈更强的通路层效应；ALS特异和FTD特异候选采用各自独立门槛，不因故事需要强行补足。该表适合用于主图leading-edge网络和后续有限候选验证，不应被解释为单基因因果证据。",
  "",
  "## 文件",
  "",
  "- `leading_edge_driver_shortlist.tsv`：10–30个主候选及证据分级。",
  "- `all_fixed_program_gene_evidence.tsv`：所有候选与未通过门槛原因。",
  "- `shortlist_evidence_long.tsv`：主候选逐数据层溯源。",
  "- `evidence_scoring_and_gating_rules.tsv`：评分与硬门槛。",
  "- `analysis_QC.tsv`和`source_file_manifest.tsv`：QC及输入文件指纹。"
)
writeLines(summary_lines, file.path(out_dir, "README_zh.md"))

cat("Completed leading-edge evidence v2\n")
print(qc)
print(shortlist[, .(
  shortlist_rank, gene, spectrum_class, evidence_tier, evidence_score,
  original7_comparisons, same_cell_significant_comparisons,
  external_FTD_comparisons, GSE330_significant_units,
  same_donor_ATAC_concordant_axes
)])
