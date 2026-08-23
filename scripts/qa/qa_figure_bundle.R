#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
project_root <- analysis_project_root()
configure_project_library(project_root)

suppressPackageStartupMessages(library(data.table))

figure_dir <- file.path(project_root, "figures")
source_dir <- file.path(project_root, "data", "source_data")
qa_dir <- file.path(project_root, "outputs", "qa")
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

main_stems <- paste0("Figure", 1:6)
extended_stems <- c("ExtendedDataFigure13", "ExtendedDataFigure14")
required <- rbindlist(list(
  CJ(stem = main_stems, extension = c("pdf", "svg", "tiff")),
  CJ(stem = extended_stems, extension = c("pdf", "svg", "tiff"))
))
# Figure 1 is the author-final raster overview and has no locked vector master.
required <- required[!(stem == "Figure1" & extension %chin% c("pdf", "svg"))]
required[, path := file.path("figures", sprintf("%s.%s", stem, extension))]
required[, full_path := file.path(project_root, path)]
required[, `:=`(
  exists = file.exists(full_path),
  bytes = fifelse(file.exists(full_path), file.info(full_path)$size, NA_real_)
)]
required[, full_path := NULL]

source_files <- list.files(source_dir, pattern = "source_data\\.tsv$|dataset_roles\\.tsv$|evidence_roles\\.tsv$",
                           recursive = TRUE, full.names = TRUE)
source_qc <- rbindlist(lapply(source_files, function(path) {
  x <- fread(path, nrows = 1)
  relative_path <- substring(normalizePath(path), nchar(normalizePath(project_root)) + 2L)
  data.table(path = relative_path, columns = ncol(x), bytes = file.info(path)$size)
}), fill = TRUE)

fwrite(required, file.path(qa_dir, "figure_export_file_QC.tsv"), sep = "\t")
fwrite(source_qc, file.path(qa_dir, "figure_source_data_QC.tsv"), sep = "\t")
summary <- data.table(
  metric = c("final_figure_files_expected", "final_figure_files_present", "source_tables", "empty_source_tables"),
  value = c(nrow(required), sum(required$exists), length(source_files), sum(source_qc$bytes <= 0))
)
fwrite(summary, file.path(qa_dir, "figure_bundle_QC_summary.tsv"), sep = "\t")

if (!all(required$exists) || any(required$bytes <= 0, na.rm = TRUE)) stop("Missing or empty final figure export.")
if (length(source_files) < 33L || any(source_qc$bytes <= 0)) stop("Figure Source Data bundle is incomplete.")
message("Release figure bundle QA passed")
