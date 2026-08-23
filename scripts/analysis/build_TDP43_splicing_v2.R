#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(xml2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
project_dir <- analysis_project_root()
configure_project_library(project_dir)
output_dir <- file.path(project_dir, "results/TDP43_splicing_v2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gittings_dir <- file.path(project_dir, "data/raw/Gittings_2023")
source_dir <- file.path(project_dir, "data/raw/TDP43_splicing_v2_sources")
gittings_xml <- file.path(source_dir, "PMC10412668.xml")
gittings_supp3 <- file.path(gittings_dir, "401_2023_2599_MOESM3_ESM.xlsx")
gittings_supp4 <- file.path(gittings_dir, "401_2023_2599_MOESM4_ESM.xlsx")
ma_supp1 <- file.path(
  project_dir,
  "data/raw/TDP43_reference_sets/Ma_2022_Nature_Supplementary_Table1.xlsx"
)
gse219280_gene_matrix <- file.path(
  project_dir,
  "data/raw/snRNA_geneByCell_dgCMatrix_RNA_raw_count_clean_for_manuscript.rds"
)
gse219280_published <- file.path(project_dir, "data/raw/published_results")

required_files <- c(gittings_xml, gittings_supp3, gittings_supp4, ma_supp1, gse219280_gene_matrix)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required inputs: ", paste(missing_files, collapse = "; "))
}

clean_number <- function(x) {
  x <- trimws(as.character(x))
  x[x %chin% c("", "-", "–", "—", "NA")] <- NA_character_
  suppressWarnings(as.numeric(gsub(",", "", x, fixed = TRUE)))
}

yes_flag <- function(x) {
  !is.na(x) & grepl("^yes$", trimws(as.character(x)), ignore.case = TRUE)
}

# -----------------------------------------------------------------------------
# 1. Audit whether the primary project data can support a de novo splice test.
# -----------------------------------------------------------------------------
audit_roots <- c(
  file.path(project_dir, "data/raw"),
  file.path(project_dir, "data/processed"),
  file.path(project_dir, "results")
)
all_files <- unlist(lapply(audit_roots, function(root) {
  list.files(root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
}), use.names = FALSE)

file_type <- function(path) {
  p <- tolower(path)
  fifelse(grepl("\\.(bam|cram|sam)$", p), "alignment",
    fifelse(grepl("sj\\.out\\.tab|junction", p), "junction_table",
      fifelse(grepl("\\.(fastq|fq)(\\.gz)?$", p), "raw_reads",
        fifelse(grepl("transcript|isoform|psi|leafcutter|majiq", p), "transcript_or_splicing_summary",
          fifelse(grepl("gene.*(count|matrix)|genebycell|\\.h5$|\\.h5ad$", p), "gene_level_matrix", "other")
        )
      )
    )
  )
}

file_audit <- data.table(path = all_files)
file_audit[, type := file_type(path)]
file_audit <- file_audit[type != "other"]
file_audit[, bytes := file.info(path)$size]
file_audit[, belongs_to_primary_GSE219280 := grepl(
  "GSE219280|snRNA_geneByCell",
  path,
  ignore.case = TRUE
)]
fwrite(file_audit, file.path(output_dir, "local_splicing_input_file_audit.tsv"), sep = "\t")

primary_splice_inputs <- file_audit[
  belongs_to_primary_GSE219280 == TRUE & type %chin% c(
    "alignment", "junction_table", "raw_reads", "transcript_or_splicing_summary"
  )
]

data_availability <- data.table(
  dataset = c("GSE219280 primary snRNA", "Gittings 2023 independent C9 spectrum", "Ma 2022 TDP-43-negative human neuronal nuclei"),
  local_input = c(
    "gene-by-cell count matrix plus gene-level published differential-expression tables",
    "published event coordinates, snRNA detection calls, cell-type/group junction counts and author-reported statistics",
    "published MAJIQ/LeafCutter event table"
  ),
  alignment_or_junction_input_available = c(nrow(primary_splice_inputs) > 0, TRUE, TRUE),
  de_novo_splicing_inference_possible = c(FALSE, FALSE, FALSE),
  reproducible_secondary_extraction_possible = c(FALSE, TRUE, TRUE),
  interpretation = c(
    "Not testable from the current local GSE219280 inputs: a gene-count matrix cannot recover exon-exon junctions, transcript isoforms or PSI.",
    "Formal event-level published evidence can be extracted and summarized, but this is independent published evidence rather than a reanalysis of raw reads.",
    "Formal TDP-43-loss event metrics can be extracted from the published supplementary table, but this is reference evidence rather than a new cohort analysis."
  )
)
fwrite(data_availability, file.path(output_dir, "data_availability_and_claim_audit.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 2. Extract the Gittings et al. 2023 event catalogue (Supplementary Table 3).
# -----------------------------------------------------------------------------
gittings_events <- as.data.table(read_excel(
  gittings_supp3,
  sheet = "Supplementary Table 3",
  .name_repair = "unique"
))
setnames(gittings_events, c(
  "gene",
  "ce_detected_in_liu_2019",
  "ce_junction_coordinates_liu_2019",
  "ce_detected_in_gittings_snrna",
  "distance_to_3prime_bp",
  "within_1805_bp_of_3prime"
))
gittings_events[, `:=`(
  gene = toupper(trimws(as.character(gene))),
  liu_2019_event_detected = yes_flag(ce_detected_in_liu_2019),
  gittings_snrna_event_detected = yes_flag(ce_detected_in_gittings_snrna),
  distance_to_3prime_bp = clean_number(distance_to_3prime_bp),
  source = "Gittings et al., Acta Neuropathologica, 2023",
  doi = "10.1007/s00401-023-02599-5",
  pmid = "37466726",
  pmcid = "PMC10412668",
  evidence_unit = "published splice-junction detection"
)]
fwrite(gittings_events, file.path(output_dir, "Gittings_2023_cryptic_exon_event_catalog.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 3. Parse article Table 2: exact group/cell-type CE counts for STMN2 and KALRN.
# -----------------------------------------------------------------------------
doc <- read_xml(gittings_xml)
table2_rows <- xml_find_all(doc, "//table-wrap[@id='Tab2']//tr[position() > 2]")
row_values <- lapply(table2_rows, function(node) {
  xml_text(xml_find_all(node, "./th|./td"), trim = TRUE)
})
if (!all(lengths(row_values) == 16L)) {
  stop("Unexpected Gittings Table 2 schema; expected 16 values per data row.")
}

wide <- rbindlist(lapply(row_values, function(v) as.list(v)))
setnames(wide, c(
  "cell_type",
  paste0("total_cells_", c("Control", "C9_ALS", "C9_ALS_FTD", "C9_FTD")),
  paste0("STMN2_CE_count_", c("Control", "C9_ALS", "C9_ALS_FTD", "C9_FTD")),
  paste0("KALRN_CE_count_", c("Control", "C9_ALS", "C9_ALS_FTD", "C9_FTD")),
  paste0("author_reported_CE_per_cell_", c("C9_ALS", "C9_ALS_FTD", "C9_FTD"))
))

numeric_columns <- setdiff(names(wide), "cell_type")
wide[, (numeric_columns) := lapply(.SD, clean_number), .SDcols = numeric_columns]

disease_levels <- c("Control", "C9_ALS", "C9_ALS_FTD", "C9_FTD")
counts_long <- rbindlist(lapply(disease_levels, function(group) {
  rate_col <- paste0("author_reported_CE_per_cell_", group)
  data.table(
    cell_type = wide$cell_type,
    disease_group = group,
    n_donors_article_table2_header = c(Control = 13, C9_ALS = 10, C9_ALS_FTD = 6, C9_FTD = 9)[[group]],
    n_donors_results_text_frontal = c(Control = 12, C9_ALS = 10, C9_ALS_FTD = 6, C9_FTD = 9)[[group]],
    total_cells = wide[[paste0("total_cells_", group)]],
    STMN2_CE_junction_count = wide[[paste0("STMN2_CE_count_", group)]],
    KALRN_CE_junction_count = wide[[paste0("KALRN_CE_count_", group)]],
    author_reported_CE_per_cell = if (rate_col %chin% names(wide)) wide[[rate_col]] else NA_real_
  )
}), fill = TRUE)
counts_long[, descriptive_combined_junctions_per_1000_nuclei :=
  1000 * (STMN2_CE_junction_count + KALRN_CE_junction_count) / total_cells]
counts_long[, `:=`(
  inference_status = "descriptive_only_no_new_p_value",
  caution = "Junction counts are published aggregates; nuclei and junctions are not independent donor-level observations, and a nucleus may contribute more than one event.",
  source_count_note = fifelse(
    disease_group == "Control",
    "Article Table 2 header states Control n=13, whereas Figure 1/Results and local Supplementary Table 2 list 12 frontal-cortex controls; both values are retained rather than silently reconciled.",
    "Article Table 2 header and Figure 1/Results agree on group size."
  ),
  source_table = "Gittings 2023 article Table 2",
  doi = "10.1007/s00401-023-02599-5"
)]
fwrite(counts_long, file.path(output_dir, "Gittings_2023_CE_counts_by_celltype_and_disease.tsv"), sep = "\t")

group_totals <- counts_long[cell_type == "Total", .(
  disease_group,
  n_donors_article_table2_header,
  n_donors_results_text_frontal,
  total_cells,
  STMN2_CE_junction_count,
  KALRN_CE_junction_count,
  descriptive_combined_junctions_per_1000_nuclei
)]
group_totals[, interpretation :=
  "Descriptive group aggregate only; use author-reported donor-level Wilcoxon tests for inference."]
fwrite(group_totals, file.path(output_dir, "Gittings_2023_CE_group_totals.tsv"), sep = "\t")

# Deep-sequencing sensitivity example from the same paper, Table 3.
table3_rows <- xml_find_all(doc, "//table-wrap[@id='Tab3']//tr[position() > 2]")
row_values3 <- lapply(table3_rows, function(node) {
  xml_text(xml_find_all(node, "./th|./td"), trim = TRUE)
})
if (!all(lengths(row_values3) == 8L)) {
  stop("Unexpected Gittings Table 3 schema; expected 8 values per data row.")
}
deep <- rbindlist(lapply(row_values3, function(v) as.list(v)))
setnames(deep, c(
  "cell_type", "number_of_cells", "baseline_STMN2_CE_junctions",
  "baseline_KALRN_CE_junctions", "baseline_CE_per_cell",
  "deep_STMN2_CE_junctions", "deep_KALRN_CE_junctions", "deep_CE_per_cell"
))
deep[, setdiff(names(deep), "cell_type") := lapply(.SD, clean_number), .SDcols = setdiff(names(deep), "cell_type")]
deep[, `:=`(
  subject = "C9-FTD 4",
  interpretation = "Technical-depth sensitivity example in one donor; not an independent biological replicate.",
  source_table = "Gittings 2023 article Table 3"
)]
fwrite(deep, file.path(output_dir, "Gittings_2023_deep_sequencing_sensitivity.tsv"), sep = "\t")

# Author-reported donor-level significance from the article Figure 1 caption.
reported_tests <- data.table(
  gene = c("STMN2", "STMN2", "STMN2", "KALRN"),
  brain_region = c("frontal cortex", "frontal cortex", "occipital cortex", "frontal cortex"),
  comparison = c("C9-ALS-FTD vs Control", "C9-FTD vs Control", "C9-FTD vs Control", "C9-FTD vs Control"),
  p_value = c(0.00048, 0.00079, 0.044, 0.00042),
  test = "Wilcoxon test",
  outcome = "average CE junctions detected per subject",
  statistic_source = "author-reported in Gittings 2023 Figure 1 caption",
  independently_recomputed = FALSE,
  doi = "10.1007/s00401-023-02599-5"
)
fwrite(reported_tests, file.path(output_dir, "Gittings_2023_author_reported_donor_tests.tsv"), sep = "\t")

# Numbers of frontal-cortex subjects with a detected CE, as reported in the
# article Results (subject-level detection, not a newly fitted statistical test).
subject_detection <- data.table(
  gene = rep(c("STMN2", "KALRN"), each = 4),
  disease_group = rep(c("Control", "C9_ALS", "C9_ALS_FTD", "C9_FTD"), 2),
  subjects_with_detected_CE = c(1, 4, 6, 7, 1, 1, 1, 7),
  subjects_assessed = c(12, 10, 6, 9, 12, 10, 6, 9),
  brain_region = "frontal cortex",
  statistic_source = "author-reported in Gittings 2023 Results/Figure 1",
  independently_recomputed = FALSE,
  doi = "10.1007/s00401-023-02599-5"
)
subject_detection[, proportion_subjects_with_detected_CE :=
  subjects_with_detected_CE / subjects_assessed]
fwrite(subject_detection, file.path(output_dir, "Gittings_2023_subject_level_CE_detection.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 4. Extract Ma et al. 2022 human neuronal TDP-43-loss splice-event metrics.
# -----------------------------------------------------------------------------
ma_events <- as.data.table(read_excel(ma_supp1, sheet = "Sheet1", .name_repair = "unique"))
setnames(ma_events, c(
  "gene", "chromosome", "coordinates_hg38", "MAJIQ_abs_max_delta_PSI",
  "LeafCutter_FDR", "LeafCutter_abs_max_delta_PSI"
))
ma_events[, gene := toupper(trimws(as.character(gene)))]
ma_events[!nzchar(gene), gene := NA_character_]
observed_gene_rows <- which(!is.na(ma_events$gene))
ma_events[, gene := gene[observed_gene_rows[
  findInterval(seq_len(.N), observed_gene_rows)
]]]
ma_events[, `:=`(
  source = "Ma et al., Nature, 2022, Supplementary Table 1",
  comparison = "TDP-43-negative vs TDP-43-positive human neuronal nuclei",
  doi = "10.1038/s41586-022-04424-7",
  pmid = "35197626",
  evidence_unit = "published differential splice event"
)]
fwrite(ma_events, file.path(output_dir, "Ma_2022_human_neuronal_TDP43_loss_splice_events.tsv"), sep = "\t")

ma_gene <- ma_events[, .(
  n_events = .N,
  max_MAJIQ_abs_delta_PSI = max(MAJIQ_abs_max_delta_PSI, na.rm = TRUE),
  min_LeafCutter_FDR = min(LeafCutter_FDR, na.rm = TRUE),
  max_LeafCutter_abs_delta_PSI = max(LeafCutter_abs_max_delta_PSI, na.rm = TRUE),
  coordinates_hg38 = paste(unique(coordinates_hg38), collapse = ";")
), by = gene]

# -----------------------------------------------------------------------------
# 5. Predefined key-gene evidence, with explicit claim boundaries.
# -----------------------------------------------------------------------------
key_genes <- data.table(gene = c("STMN2", "UNC13A", "KALRN"))
key_evidence <- merge(key_genes, ma_gene, by = "gene", all.x = TRUE)
key_evidence <- merge(
  key_evidence,
  gittings_events[, .(
    gene,
    liu_2019_event_detected,
    ce_junction_coordinates_liu_2019,
    gittings_snrna_event_detected,
    distance_to_3prime_bp
  )],
  by = "gene",
  all.x = TRUE
)

key_evidence[, primary_GSE219280_direct_splicing_test := "not_testable"]
key_evidence[, primary_GSE219280_reason :=
  "Only gene-level count/expression inputs are local; no read alignments, junction counts, transcript counts or PSI are available."]
key_evidence[, evidence_tier := fifelse(
  gene %chin% c("STMN2", "UNC13A"),
  "Tier A: established TDP-43-dependent cryptic-splicing target in human neurons/brain",
  "Tier B: strong human-neuronal TDP-43-loss splice event plus independent C9 snRNA junction detection"
)]
key_evidence[, project_specific_conclusion := fcase(
  gene == "STMN2",
  "Independent published C9ORF72-spectrum snRNA data provide event-level STMN2 CE support, including donor-level group tests; the present GSE219280 reanalysis itself did not measure splicing.",
  gene == "UNC13A",
  "UNC13A is a pre-established mechanistic TDP-43 target. It was detected in prior bulk junction data but not in Gittings 3-prime-biased snRNA, and it was not directly testable in the present GSE219280 inputs.",
  gene == "KALRN",
  "Independent published C9ORF72-spectrum snRNA data detected KALRN CE junctions, especially in C9-FTD. This supports a candidate TDP-43 dysfunction target but is not equivalent to an individual causal perturbation experiment in this project."
)]
key_evidence[, permitted_manuscript_label := fcase(
  gene == "STMN2", "independent published event-level cryptic-exon support",
  gene == "UNC13A", "established literature-based TDP-43 cryptic-splicing target; not directly validated here",
  gene == "KALRN", "independent published event-level candidate cryptic-exon support"
)]
key_evidence[, prohibited_manuscript_claim :=
  "Do not state that our GSE219280 gene-expression or GSE219279 snATAC analysis demonstrated cryptic-exon inclusion or TDP-43 loss of function."]
fwrite(key_evidence, file.path(output_dir, "predefined_key_gene_splicing_evidence.tsv"), sep = "\t")

claim_boundaries <- data.table(
  item = c(
    "Primary GSE219280 splicing",
    "Gittings 2023 STMN2",
    "Gittings 2023 KALRN",
    "UNC13A",
    "Gene-expression module enrichment",
    "snATAC GeneScore/accessibility"
  ),
  status = c(
    "not testable with current inputs",
    "formal independent published event-level support",
    "formal independent published event-level support",
    "established mechanistic literature evidence; not detected in Gittings snRNA and not tested in primary input",
    "supportive context only",
    "supportive chromatin context only"
  ),
  allowed_wording = c(
    "The available GSE219280 gene-count matrix did not permit junction- or isoform-level testing.",
    "An independent C9ORF72-spectrum snRNA study detected STMN2 CE junctions and reported donor-level group differences.",
    "An independent C9ORF72-spectrum snRNA study detected KALRN CE junctions, most prominently in C9-FTD.",
    "UNC13A was treated as a predefined, experimentally established TDP-43 cryptic-splicing target.",
    "The transcriptional data were consistent with perturbation of a predefined TDP-43 target program.",
    "Chromatin accessibility at candidate loci provided orthogonal regulatory context."
  ),
  disallowed_wording = c(
    "We detected cryptic exons in GSE219280.",
    "Our cohort proves TDP-43 loss of function through STMN2 splicing.",
    "KALRN CE is an ALS-versus-FTD causal switch.",
    "UNC13A cryptic splicing was replicated in GSE219280.",
    "Expression change proves mis-splicing.",
    "ATAC GeneScore proves cryptic-exon inclusion or RNA splicing."
  )
)
fwrite(claim_boundaries, file.path(output_dir, "manuscript_claim_boundaries.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 6. Source manifest and concise decision report.
# -----------------------------------------------------------------------------
source_manifest <- data.table(
  source = c(
    "Gittings 2023 full-text JATS XML",
    "Gittings 2023 Supplementary Table 3",
    "Gittings 2023 Supplementary Table 4",
    "Ma 2022 Supplementary Table 1",
    "Primary GSE219280 gene-by-cell matrix"
  ),
  local_path = c(gittings_xml, gittings_supp3, gittings_supp4, ma_supp1, gse219280_gene_matrix),
  source_url = c(
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10412668/fullTextXML",
    "https://doi.org/10.1007/s00401-023-02599-5",
    "https://doi.org/10.1007/s00401-023-02599-5",
    "https://doi.org/10.1038/s41586-022-04424-7",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219280"
  ),
  doi = c(
    "10.1007/s00401-023-02599-5",
    "10.1007/s00401-023-02599-5",
    "10.1007/s00401-023-02599-5",
    "10.1038/s41586-022-04424-7",
    NA_character_
  ),
  retrieved_or_available_date = as.character(Sys.Date())
)
source_manifest[, bytes := file.info(local_path)$size]
source_manifest[, md5 := unname(tools::md5sum(local_path))]
fwrite(source_manifest, file.path(output_dir, "source_manifest.tsv"), sep = "\t")

report <- c(
  "# TDP-43 RNA-processing / cryptic-exon evidence audit",
  "",
  "## Decision",
  "",
  "The current local GSE219280 inputs cannot support a de novo cryptic-exon, splice-junction, transcript-isoform or PSI analysis. The main input is a gene-by-cell count matrix; the required read-alignment or junction-level information is absent.",
  "",
  "Formal TDP-43 splicing evidence can nevertheless be included as a clearly labelled, independent published validation layer:",
  "",
  "- STMN2: Tier A. Gittings 2023 provides direct C9ORF72-spectrum snRNA junction evidence and author-reported donor-level group tests; Ma 2022 provides a strong human neuronal TDP-43-loss splice event.",
  "- UNC13A: Tier A mechanistic literature evidence, but not a direct validation in the present project. Gittings detected the event in earlier bulk junction data but not in its 3-prime-biased snRNA data; absence of detection is not evidence of absence.",
  "- KALRN: Tier B. Ma 2022 provides strong human-neuronal TDP-43-loss differential splicing, and Gittings 2023 independently detected CE junctions in C9ORF72 snRNA, especially C9-FTD. It should not be called a proven ALS-versus-FTD causal switch.",
  "",
  "Source audit note: Gittings article Table 2 labels the control group as n=13, whereas its Figure 1/Results and local Supplementary Table 2 indicate 12 frontal-cortex controls. Both source values are retained in the machine-readable tables; donor-level narrative follows the article Results (1/12 controls).",
  "",
  "## Manuscript placement",
  "",
  "Use this as independent published event-level support for the TDP-43 axis, not as a new primary-cohort splicing result. The article title and abstract should not claim that this project directly measured cryptic exons unless read-level GSE219280 data are later obtained and analysed.",
  "",
  "## Key sources",
  "",
  "- Gittings et al. Acta Neuropathologica (2023), DOI: 10.1007/s00401-023-02599-5, PMID: 37466726.",
  "- Ma et al. Nature (2022), DOI: 10.1038/s41586-022-04424-7, PMID: 35197626.",
  "- Brown et al. Nature (2022), DOI: 10.1038/s41586-022-04436-3, PMID: 35197628. Brown et al. independently established UNC13A TDP-43-dependent mis-splicing; its metadata are cited for evidence triangulation but no Brown source-data table was reanalysed here."
)
writeLines(report, file.path(output_dir, "TDP43_splicing_evidence_decision.md"), useBytes = TRUE)

manuscript_text <- c(
  "# Manuscript-ready bounded text",
  "",
  "## Methods",
  "",
  "Because the primary GSE219280 input available for reanalysis was a gene-by-cell count matrix rather than read alignments or splice-junction counts, cryptic-exon inclusion was not inferred from gene expression. We instead predefined STMN2, UNC13A and KALRN as TDP-43-associated splicing targets and extracted event-level evidence from published human neuronal TDP-43-loss data (Ma et al., Nature 2022; DOI: 10.1038/s41586-022-04424-7) and an independent C9ORF72-spectrum single-nucleus study (Gittings et al., Acta Neuropathologica 2023; DOI: 10.1007/s00401-023-02599-5). Published aggregate junction counts were summarized descriptively; inferential p values were retained only when the original authors reported donor-level tests.",
  "",
  "## Results",
  "",
  "Human neuronal nuclei depleted of nuclear TDP-43 showed large differential-splicing effects for STMN2 (MAJIQ |delta PSI|=0.682; LeafCutter |delta PSI|=0.706, FDR=4.82e-6), UNC13A (0.779; 0.836, FDR=1.17e-8) and KALRN (0.709; 0.831, FDR=2.26e-12). In an independent C9ORF72-spectrum snRNA-seq study, STMN2 cryptic-exon junctions were detected in frontal cortex from 4/10 C9-ALS, 6/6 C9-ALS-FTD and 7/9 C9-FTD subjects, compared with 1/12 controls; the authors reported increases for C9-ALS-FTD versus control (P=4.8e-4) and C9-FTD versus control (P=7.9e-4). KALRN cryptic-exon junctions were detected in 7/9 C9-FTD subjects versus 1/12 controls (author-reported P=4.2e-4), but only 1/10 C9-ALS and 1/6 C9-ALS-FTD subjects. UNC13A cryptic-exon inclusion is experimentally established in human neurons and brain, but was not detected in the independent 10x 3-prime snRNA dataset and could not be tested directly in our gene-level input. Accordingly, these data provide independent published event-level support for a TDP-43 RNA-processing axis, rather than direct cryptic-exon validation in the primary cohort.",
  "",
  "## Required limitation",
  "",
  "The available primary single-nucleus expression matrix does not retain exon-junction information. Therefore, transcriptional module enrichment and chromatin accessibility were interpreted as orthogonal context and not as evidence of cryptic-exon inclusion or TDP-43 nuclear loss."
)
writeLines(manuscript_text, file.path(output_dir, "manuscript_ready_bounded_text.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

cat("TDP43 splicing evidence v2 complete\n")
cat("Primary direct splicing inputs found:", nrow(primary_splice_inputs), "\n")
cat("Gittings event catalogue rows:", nrow(gittings_events), "\n")
cat("Ma event rows:", nrow(ma_events), "\n")
cat("Output:", output_dir, "\n")
