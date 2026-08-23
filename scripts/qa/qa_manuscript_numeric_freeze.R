#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
base <- analysis_project_root()
configure_project_library(base)
suppressPackageStartupMessages(library(data.table))

sd <- file.path(base, "data", "source_data", "main_figures")
manu <- file.path(base, "manuscript")
checks <- list()
add_check <- function(name, observed, expected, tolerance = 0) {
  ok <- if (is.numeric(observed) && is.numeric(expected)) {
    length(observed) == length(expected) && all(abs(observed - expected) <= tolerance)
  } else identical(as.character(observed), as.character(expected))
  checks[[length(checks) + 1L]] <<- data.table(
    check = name,
    observed = paste(observed, collapse = ";"),
    expected = paste(expected, collapse = ";"),
    pass = ok
  )
}

f1 <- fread(file.path(sd, "Figure1_dataset_roles.tsv"))
add_check("Figure1 dataset-role rows", nrow(f1), 8)

f2a <- fread(file.path(sd, "Figure2_panel_a_source_data.tsv"))
add_check("Figure2 discovery axes", nrow(f2a), 63)
add_check("Figure2 global-FDR axes", sum(f2a$global_primary_FDR < 0.05), 35)
f2c <- fread(file.path(sd, "Figure2_panel_c_source_data.tsv"))
add_check("Figure2 threshold rows direction retained", sum(f2c$direction_retained), 126)

f3a <- fread(file.path(sd, "Figure3_panel_a_source_data.tsv"))
add_check("Figure3 exact axes", sum(f3a$axes), 63)
add_check("Figure3 same-direction axes", sum(f3a$same_direction_n), 39)
add_check("Figure3 strict replicated axes", sum(f3a$both_global_FDR_same_direction_n), 13)
f3c <- fread(file.path(sd, "Figure3_panel_c_source_data.tsv"))
add_check("Figure3 Gittings strict axes", nrow(f3c), 68)
add_check("Figure3 Gittings >=90% LOO", sum(f3c$direction_concordant_fraction >= 0.9), 68)
add_check("Figure3 Gittings 100% LOO", sum(f3c$direction_concordant_fraction == 1), 66)
add_check("Figure3 Gittings donor omissions", unique(f3c$total_donors), 37)
f3d <- fread(file.path(sd, "Figure3_panel_d_source_data.tsv"))
add_check("Figure3 GRN axes", nrow(f3d), 3)
add_check("Figure3 GRN donors", unique(f3d$n_donors), 9)
add_check("Figure3 GRN LOO retention", sum(f3d$direction_retained_fraction == 1), 3)

f4a <- fread(file.path(sd, "Figure4_panel_a_source_data.tsv"))
add_check("Figure4 candidate genes", uniqueN(f4a$gene), 17)
f4c <- fread(file.path(sd, "Figure4_panel_c_source_data.tsv"))
event <- setNames(f4c$max_LeafCutter_abs_delta_PSI, f4c$gene)
add_check("Figure4 STMN2 delta PSI", event[["STMN2"]], 0.706, 1e-12)
add_check("Figure4 KALRN delta PSI", event[["KALRN"]], 0.831, 1e-12)
add_check("Figure4 UNC13A delta PSI", event[["UNC13A"]], 0.836, 1e-12)

f5a <- fread(file.path(sd, "Figure5_panel_a_source_data.tsv"))
prot <- f5a[cell_type == "ASC" & family == "Protein stabilization"]
add_check("Figure5 astrocyte beta", prot$beta_per_pTDP_level, -0.840209208987552, 1e-12)
add_check("Figure5 astrocyte CI low", prot$ci_low, -1.18204396526168, 1e-12)
add_check("Figure5 astrocyte CI high", prot$ci_high, -0.498374452713428, 1e-12)
add_check("Figure5 astrocyte global FDR", prot$global_FDR, 0.039908208056746, 1e-12)
f5b <- fread(file.path(sd, "Figure5_panel_b_source_data.tsv"))
hits <- f5b[target_global_FDR < 0.05]
add_check("Figure5 candidate RNA hits", nrow(hits), 8)
add_check("Figure5 candidate RNA 100% LOO", sum(hits$loo_same_direction_fraction == 1), 8)
bag3 <- hits[feature == "BAG3" & cell_type == "MG"]
add_check("Figure5 microglial BAG3 transcriptome FDR", bag3$FDR, 0.04426075, 1e-8)
f5c <- fread(file.path(sd, "Figure5_panel_c_source_data.tsv"))
dar <- f5c[, .N, by = cell_type][order(cell_type)]
add_check("Figure5 total DAR", nrow(f5c), 409)
add_check("Figure5 MG DAR", dar[cell_type == "MG", N], 3)
add_check("Figure5 ODC DAR", dar[cell_type == "ODC", N], 406)
f5d <- fread(file.path(sd, "Figure5_panel_d_source_data.tsv"))
add_check("Figure5 peak-gene links", nrow(f5d), 7)
add_check("Figure5 rho minimum", min(f5d$rho), -0.854945054945055, 1e-12)
add_check("Figure5 rho maximum", max(f5d$rho), 0.89010989010989, 1e-12)
f5e <- fread(file.path(sd, "Figure5_panel_e_source_data.tsv"))
add_check("Figure5 complete chains", sum(f5e$criterion == "Full three-part chain" & f5e$passed), 0)

md <- list.files(manu, pattern = "\\.md$", full.names = TRUE)
content <- paste(unlist(lapply(md, readLines, warn = FALSE)), collapse = "\n")
add_check("No FINAL_VALUE_PENDING in manuscript markdown", grepl("FINAL_VALUE_PENDING", content, fixed = TRUE), FALSE)
add_check("RNA-ATAC causal boundary present", grepl("not enhancer-to-gene regulation or an ATAC-to-RNA causal cascade", content, fixed = TRUE), TRUE)
add_check("LOO boundary present", grepl("internal robustness test and not an additional validation cohort", content, fixed = TRUE), TRUE)

out <- rbindlist(checks)
dir.create(file.path(base, "outputs", "qa"), recursive = TRUE, showWarnings = FALSE)
fwrite(out, file.path(base, "outputs", "qa", "manuscript_numeric_freeze_QC.tsv"), sep = "\t")
cat(sprintf("Manuscript numeric freeze QA: %d/%d checks passed\n", sum(out$pass), nrow(out)))
if (!all(out$pass)) {
  print(out[pass == FALSE])
  quit(status = 1)
}
