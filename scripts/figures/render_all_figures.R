#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
script_dir <- dirname(script_path)

scripts <- c(
  "redesign_Figure2_2026.R",
  "redesign_Figure3_2026.R",
  "redesign_Figure4_2026.R",
  "redesign_Figure5_2026.R",
  "redesign_Figure6_2026.R",
  "render_ExtendedDataFigure13.R",
  "redesign_ExtendedDataFigure14_2026.R"
)

for (script in scripts) {
  path <- file.path(script_dir, script)
  message("Rendering with ", script)
  code <- system2(file.path(R.home("bin"), "Rscript"), path)
  if (!identical(code, 0L)) stop("Figure rendering failed: ", script)
}

message("All scripted figures rendered successfully.")

