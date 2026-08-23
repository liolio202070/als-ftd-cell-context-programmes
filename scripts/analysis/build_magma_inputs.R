#!/usr/bin/env Rscript
# 构建 MAGMA 输入：基因集 .gsa 文件（7 机制家族 + 17 冻结候选 + 13 条严格复制轴 leading-edge）
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
root <- analysis_project_root()
configure_project_library(root)
.libPaths(c(file.path(root, "cache/R_library"), .libPaths()))
suppressPackageStartupMessages(library(data.table))
out <- file.path(root, "external/magma")

mem <- fread(file.path(root, "results/GSE219280_same_celltype_fixed_programs/predefined_gene_set_membership.tsv"))
fam_sets <- mem[, .(genes = list(unique(gene))), by = family_id]

cand17 <- c("BAG3","HSP90AA1","HSP90AB1","CHORDC1","ARHGAP35","PARK7","CCT5","ARHGAP39",
  "COX5B","COX7C","MAP1B","MAPT","BCR","CCT8","COA3","DOCK7","STARD13")

rep <- fread(file.path(root, "results/submission_cross_cohort_evidence_v1/exact_estimand_program_replication.tsv"))
strict <- rep[replication_class == "bidirectional_global_FDR_concordant",
              .(estimand, cell_type, family_id)]
# 与 build_cross_cohort_evidence.R 相同的过滤与映射
fg1 <- fread(file.path(root, "results/GSE219280_same_celltype_fixed_programs/fixed_family_fgsea_results.tsv"))
fg1 <- fg1[analysis_role == "primary"]
fg1[, `:=`(estimand = fcase(contrast == "ALS_vs_Ctrl_mean", "ALS_vs_Control",
                            contrast == "FTD_vs_Ctrl_mean", "FTD_vs_Control",
                            contrast == "ALS_vs_FTD_mean", "ALS_vs_FTD"),
          cell_type_h = fcase(cell_type == "Astro", "Astrocyte",
                              cell_type == "Oligo", "Oligodendrocyte", default = cell_type))]
fg2 <- fread(file.path(root, "results/Gittings_2023_C9_snRNA/fixed_7_family_fgsea_results.tsv"))
fg2 <- fg2[analysis_role == "primary_frontal"]
fg2[, `:=`(estimand = fcase(contrast == "C9_ALS_vs_Ctrl_Frontal", "ALS_vs_Control",
                            contrast == "C9_FTD_vs_Ctrl_Frontal", "FTD_vs_Control",
                            contrast == "C9_ALS_vs_C9_FTD_Frontal_batch_bridged", "ALS_vs_FTD"),
          cell_type_h = fcase(cell_type == "Astro", "Astrocyte",
                              cell_type == "Oligo", "Oligodendrocyte", default = cell_type))]
le_union <- function(fg) {
  m <- merge(strict, fg[, .(estimand, cell_type = cell_type_h, family_id = set_id, leadingEdge)],
             by = c("estimand", "cell_type", "family_id"))
  unique(unlist(strsplit(m$leadingEdge[!is.na(m$leadingEdge)], "/")))
}
le1 <- le_union(fg1); le2 <- le_union(fg2)
le_genes <- union(le1, le2)
cat("严格复制轴 leading-edge：发现", length(le1), "个 + Gittings", length(le2), "个 = 并集", length(le_genes), "\n")

gsa_lines <- c()
add_set <- function(name, genes) {
  g <- unique(na.omit(genes)); g <- g[nzchar(g)]
  if (length(g) >= 5) gsa_lines <<- c(gsa_lines, paste(c(name, g), collapse = " "))
  else cat("跳过（<5 基因）：", name, "n =", length(g), "\n")
}
for (i in seq_len(nrow(fam_sets))) add_set(paste0("FAM_", fam_sets$family_id[i]), fam_sets$genes[[i]])
add_set("CAND17_frozen", cand17)
add_set("LE_strict13_pooled", le_genes)
# 预设合并层（与二项检验一致）
hi <- c("axon_development","small_GTPase","synapse_junction")
add_set("FAM_pool3_transferable", unlist(fam_sets[family_id %in% hi, genes]))

writeLines(gsa_lines, file.path(out, "program_gene_sets.gsa"))
cat("写入", length(gsa_lines), "个基因集 → program_gene_sets.gsa\n")
for (l in gsa_lines) cat(" ", strsplit(l, " ")[[1]][1], "n =", length(strsplit(l, " ")[[1]]) - 1, "\n")
