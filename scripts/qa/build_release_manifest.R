#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
project_root <- analysis_project_root()
configure_project_library(project_root)
suppressPackageStartupMessages(library(data.table))

metadata_dir <- file.path(project_root, "metadata")
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(metadata_dir, "FILE_MANIFEST.tsv")
checksum_path <- file.path(metadata_dir, "SHA256SUMS")

files <- list.files(project_root, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
files <- files[!dir.exists(files)]
files <- files[!grepl("(^|/)\\.git(/|$)", files)]
files <- files[!grepl("(^|/)__pycache__(/|$)|\\.pyc$", files)]
files <- files[!grepl("/outputs/figures/", files)]
files <- setdiff(normalizePath(files, mustWork = TRUE), normalizePath(c(manifest_path, checksum_path), mustWork = FALSE))
relative <- substring(files, nchar(normalizePath(project_root)) + 2L)

sha_command <- if (nzchar(Sys.which("shasum"))) {
  c(Sys.which("shasum"), "-a", "256")
} else if (nzchar(Sys.which("sha256sum"))) {
  c(Sys.which("sha256sum"))
} else {
  stop("Neither shasum nor sha256sum is available.")
}

sha256 <- vapply(files, function(path) {
  output <- system2(sha_command[1], c(sha_command[-1], shQuote(path)), stdout = TRUE)
  sub("[[:space:]].*$", "", output[1])
}, character(1))

category <- sub("/.*$", "", relative)
category[!grepl("/", relative)] <- "root"
manifest <- data.table(
  path = relative,
  bytes = as.numeric(file.info(files)$size),
  sha256 = sha256,
  category = category
)[order(path)]

fwrite(manifest, manifest_path, sep = "\t", quote = FALSE)
writeLines(sprintf("%s  %s", manifest$sha256, manifest$path), checksum_path)
message("Release manifest written for ", nrow(manifest), " files.")
