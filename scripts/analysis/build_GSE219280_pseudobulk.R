script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
project_root <- analysis_project_root()
configure_project_library(project_root)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

matrix_path <- file.path(project_root, "data", "raw", "snRNA_geneByCell_dgCMatrix_RNA_raw_count_clean_for_manuscript.rds")
metadata_path <- file.path(project_root, "data", "metadata", "GSE219280_nuclei_metadata.tsv")
output_rds <- file.path(project_root, "data", "processed", "GSE219280_major_celltype_pseudobulk.rds")
output_meta <- file.path(project_root, "data", "processed", "GSE219280_major_celltype_pseudobulk_samples.tsv")
output_qc <- file.path(project_root, "results", "GSE219280_pseudobulk_QC_summary.tsv")
output_audit <- file.path(project_root, "logs", "GSE219280_pseudobulk_build_audit.txt")

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(matrix_path)
meta <- fread(metadata_path)

stopifnot(inherits(counts, "dgCMatrix"))
stopifnot(ncol(counts) == nrow(meta))
matrix_barcode <- sub("^.*_([ACGT]+-[0-9]+)$", "\\1", colnames(counts))
stopifnot(identical(matrix_barcode, meta$cell_barcode))
stopifnot(!anyNA(meta$annotation_major_cell_type))

meta[, pseudobulk_id := paste(sample, annotation_major_cell_type, sep = "__")]
pseudobulk_levels <- unique(meta$pseudobulk_id)
group_index <- match(meta$pseudobulk_id, pseudobulk_levels)

membership <- sparseMatrix(
  i = seq_len(nrow(meta)),
  j = group_index,
  x = 1,
  dims = c(nrow(meta), length(pseudobulk_levels))
)

pseudobulk_counts <- counts %*% membership
colnames(pseudobulk_counts) <- pseudobulk_levels
rownames(pseudobulk_counts) <- rownames(counts)

pseudobulk_meta <- meta[, .(
  sample = first(sample),
  donor_id = first(donor_id),
  diagnosis = first(diagnosis),
  brain_region = first(brain_region),
  sex = first(sex),
  age = first(age),
  sequencing_batch = first(sequencing_batch),
  annotation_cell_class = first(annotation_cell_class),
  annotation_major_cell_type = first(annotation_major_cell_type),
  n_nuclei = .N,
  median_UMI_raw = median(num_UMI_raw),
  median_genes_raw = median(num_genes_raw),
  median_percent_mito = median(percent_mitochondrial_reads)
), by = pseudobulk_id]

setkey(pseudobulk_meta, pseudobulk_id)
pseudobulk_meta <- pseudobulk_meta[pseudobulk_levels]
pseudobulk_meta[, library_size_cellbender := as.numeric(colSums(pseudobulk_counts))]
stopifnot(identical(pseudobulk_meta$pseudobulk_id, colnames(pseudobulk_counts)))

qc <- pseudobulk_meta[, .(
  pseudobulk_columns = .N,
  donors_any_nuclei = uniqueN(donor_id),
  donors_ge_20_nuclei = uniqueN(donor_id[n_nuclei >= 20]),
  donors_ge_30_nuclei = uniqueN(donor_id[n_nuclei >= 30]),
  donors_ge_50_nuclei = uniqueN(donor_id[n_nuclei >= 50]),
  minimum_nuclei = as.numeric(min(n_nuclei)),
  median_nuclei = as.numeric(median(n_nuclei)),
  maximum_nuclei = as.numeric(max(n_nuclei))
), by = .(annotation_major_cell_type, diagnosis, brain_region)]
setorder(qc, annotation_major_cell_type, diagnosis, brain_region)

integer_like <- all(abs(counts@x - round(counts@x)) < 1e-8)
nonnegative <- all(counts@x >= 0)
total_original <- sum(counts)
total_pseudobulk <- sum(pseudobulk_counts)

saveRDS(
  list(
    counts = pseudobulk_counts,
    samples = as.data.frame(pseudobulk_meta),
    source_matrix_md5 = "6571e4dc07e42e5d5007c49f59ab82b8",
    aggregation_level = "sample x annotation_major_cell_type"
  ),
  output_rds,
  compress = "gzip"
)
fwrite(pseudobulk_meta, output_meta, sep = "\t", quote = FALSE, na = "")
fwrite(qc, output_qc, sep = "\t", quote = FALSE, na = "")

audit <- c(
  sprintf("input_dimensions\t%d genes x %d nuclei", nrow(counts), ncol(counts)),
  sprintf("output_dimensions\t%d genes x %d pseudobulk columns", nrow(pseudobulk_counts), ncol(pseudobulk_counts)),
  sprintf("integer_like_nonzero_counts\t%s", integer_like),
  sprintf("nonnegative_counts\t%s", nonnegative),
  sprintf("total_input_counts\t%.0f", total_original),
  sprintf("total_pseudobulk_counts\t%.0f", total_pseudobulk),
  sprintf("count_conservation\t%s", isTRUE(all.equal(total_original, total_pseudobulk))),
  sprintf("minimum_nuclei_per_pseudobulk\t%d", min(pseudobulk_meta$n_nuclei)),
  sprintf("median_nuclei_per_pseudobulk\t%.1f", median(pseudobulk_meta$n_nuclei)),
  sprintf("maximum_nuclei_per_pseudobulk\t%d", max(pseudobulk_meta$n_nuclei))
)
writeLines(audit, output_audit)
cat(paste(audit, collapse = "\n"), "\n")
