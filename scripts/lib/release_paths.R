# Shared path helpers for the public release.

analysis_project_root <- function() {
  configured <- Sys.getenv("ALS_FTD_PROJECT_ROOT", unset = "")
  if (nzchar(configured)) return(normalizePath(configured, mustWork = FALSE))

  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script_arg)) {
    stop("Run with Rscript or set ALS_FTD_PROJECT_ROOT to the reconstructed project directory.")
  }
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
  normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
}

configure_project_library <- function(project_root) {
  configured <- Sys.getenv("ALS_FTD_R_LIBRARY", unset = "")
  candidate <- if (nzchar(configured)) configured else file.path(project_root, "cache", "R_library")
  if (dir.exists(candidate)) .libPaths(c(candidate, .libPaths()))
  invisible(candidate)
}
