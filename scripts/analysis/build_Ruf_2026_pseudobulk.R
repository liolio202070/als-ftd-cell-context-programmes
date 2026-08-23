suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
root <- analysis_project_root()
configure_project_library(root)
raw <- file.path(root, "data/raw/Ruf_2026_multiome")
out <- file.path(root, "results/Ruf_2026_frozen_validation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

metadata <- fread(file.path(raw, "Multiome_Dataset_Metadata.txt"))
rna_barcodes <- fread(file.path(raw, "Multiome_Dataset_RNA_barcodes.tsv"),
                      header = FALSE)[[1]]
atac_barcodes <- fread(file.path(raw, "Multiome_Dataset_ATAC_barcodes.tsv"),
                       header = FALSE)[[1]]
stopifnot(identical(rna_barcodes, atac_barcodes),
          uniqueN(metadata$CellId) == nrow(metadata))
metadata <- metadata[match(rna_barcodes, CellId)]
stopifnot(!anyNA(metadata$CellId), identical(metadata$CellId, rna_barcodes))

metadata[, group := paste(ID, WNN_L2, sep = "|")]
groups <- sort(unique(metadata$group))
group_code <- match(metadata$group, groups)
fwrite(data.table(group_index = group_code),
       file.path(out, "cell_to_pseudobulk_group.txt"), col.names = FALSE)

covariates <- fread(file.path(out, "donor_covariates.tsv"))
sample_meta <- unique(metadata[, .(
  sample_id = group, donor_id = ID, cell_type = WNN_L2,
  diagnosis = Case, case_type = Case_Type, C9ORF72, sex = Sex
)])
technical_meta <- metadata[, .(
  median_RNA_UMI = median(nCount_RNA, na.rm = TRUE),
  median_RNA_features = median(nFeature_RNA, na.rm = TRUE),
  median_percent_mito = median(percent_mito, na.rm = TRUE),
  median_ATAC_counts = median(nCount_ATAC, na.rm = TRUE),
  median_ATAC_features = median(nFeature_ATAC, na.rm = TRUE),
  median_ATAC_fragments = median(atac_fragments, na.rm = TRUE),
  median_ATAC_peak_region_fragments = median(atac_peak_region_fragments, na.rm = TRUE)
), by = .(sample_id = group)]
sample_meta <- merge(sample_meta, covariates[, .(donor_id, age, RIN, PMI_hours)],
                     by = "donor_id", all.x = TRUE, sort = FALSE)
sample_meta <- merge(sample_meta, technical_meta, by = "sample_id",
                     all.x = TRUE, sort = FALSE)
sample_meta <- sample_meta[match(groups, sample_id)]
stopifnot(nrow(sample_meta) == length(groups),
          identical(sample_meta$sample_id, groups))
sample_meta[, binary_group_index := seq_len(.N)]
sample_meta[, n_nuclei := tabulate(group_code, nbins = length(groups))]
fwrite(sample_meta, file.path(out, "pseudobulk_sample_metadata.tsv"), sep = "\t")

cpp <- file.path(root, "scripts/stream_matrixmarket_pseudobulk.cpp")
exe <- file.path(root, "cache/Ruf_2026_tmp/stream_matrixmarket_pseudobulk")
dir.create(dirname(exe), recursive = TRUE, showWarnings = FALSE)
compile <- system2("clang++", c("-O3", "-std=c++17", shQuote(cpp), "-o", shQuote(exe)))
if (compile != 0L || !file.exists(exe)) stop("C++ streaming aggregator compilation failed")

feature_files <- c(
  RNA = file.path(raw, "Multiome_Dataset_RNA_features.tsv"),
  ATAC = file.path(raw, "Multiome_Dataset_ATAC_features.tsv")
)
expected_features <- c(RNA = 35367L, ATAC = 200791L)
binary_files <- character()

for (modality in names(feature_files)) {
  features <- fread(feature_files[[modality]], header = FALSE)[[1]]
  stopifnot(length(features) == expected_features[[modality]], !anyDuplicated(features))
  matrix_file <- file.path(raw, sprintf("Multiome_Dataset_%s_counts_raw.mtx", modality))
  binary_file <- file.path(out, sprintf("%s_pseudobulk_counts.uint32.bin", modality))
  args <- c(shQuote(matrix_file),
            shQuote(file.path(out, "cell_to_pseudobulk_group.txt")),
            shQuote(binary_file))
  stream_log <- file.path(out, sprintf("%s_stream_aggregation.log", modality))
  command <- paste(shQuote(exe), paste(args, collapse = " "),
                   ">", shQuote(stream_log), "2>&1")
  status <- system(command)
  if (status != 0L) stop(modality, " streaming aggregation failed")
  expected_bytes <- as.double(length(features)) * length(groups) * 4
  if (file.info(binary_file)$size != expected_bytes)
    stop(modality, " binary pseudobulk size mismatch")
  binary_files[[modality]] <- binary_file
}

read_group_columns <- function(path, n_features, group_indices) {
  con <- file(path, "rb")
  on.exit(close(con))
  result <- matrix(0L, nrow = n_features, ncol = length(group_indices))
  for (j in seq_along(group_indices)) {
    seek(con, where = as.double(group_indices[j] - 1L) * n_features * 4,
         origin = "start")
    value <- readBin(con, integer(), n = n_features, size = 4,
                     signed = FALSE, endian = "little")
    if (length(value) != n_features) stop("Short binary read")
    if (anyNA(value)) stop("Pseudobulk count exceeded R signed-integer range")
    result[, j] <- value
  }
  result
}

for (modality in names(feature_files)) {
  features <- fread(feature_files[[modality]], header = FALSE)[[1]]
  for (ct in unique(sample_meta$cell_type)) {
    sm <- sample_meta[cell_type == ct]
    mat <- read_group_columns(binary_files[[modality]], length(features),
                              sm$binary_group_index)
    dimnames(mat) <- list(features, sm$sample_id)
    saveRDS(mat, file.path(out, sprintf("%s_pseudobulk_counts_%s.rds", modality, ct)),
            compress = FALSE)
    rm(mat); gc()
  }
}

qc <- data.table(
  metric = c("n_nuclei", "n_donors", "n_donor_celltype_samples",
             "rna_features", "atac_features", "cell_types",
             "streaming_uint32_no_local_tmp"),
  value = c(nrow(metadata), uniqueN(metadata$ID), length(groups), 35367L,
            200791L, uniqueN(metadata$WNN_L2), 1L),
  pass = c(nrow(metadata) == 180016, uniqueN(metadata$ID) == 79,
           length(groups) == 79L * 6L, TRUE, TRUE,
           uniqueN(metadata$WNN_L2) == 6L, TRUE)
)
fwrite(qc, file.path(out, "RNA_pseudobulk_QC.tsv"), sep = "\t")
if (!all(qc$pass)) stop("Ruf 2026 streaming pseudobulk QC failed")
writeLines(capture.output(sessionInfo()), file.path(out, "pseudobulk_sessionInfo.txt"))
