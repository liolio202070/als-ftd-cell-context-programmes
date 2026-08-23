#!/usr/bin/env Rscript
# 统计修复③：三队列迁移门槛可达性/功效分析 + GSE212630 t 分布 CI
# 问题1（审稿）：三队列零迁移是否只是功效不足？——在"Ruf 效应为真"的乐观情景下计算门槛可达概率
# 问题2（审稿）：14–17 供体模型的正态近似 95%CI 偏乐观——用 t 分布重建（df 由 F/p 数值反推）
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
root <- analysis_project_root()
configure_project_library(root)
.libPaths(c(file.path(root, "cache/R_library"), .libPaths()))
suppressPackageStartupMessages(library(data.table))

rf <- file.path(root, "results/Ruf_2026_frozen_validation")
three <- fread(file.path(rf, "three_cohort_ALS_programme_replication.tsv"))
ruf <- fread(file.path(rf, "RNA_7_family_results.tsv"))
out <- file.path(root, "results/statistical_robustness_v1")

## ---------- Part A：三队列门槛可达性 ----------
zfrom <- function(nes, fdr) {
  fdr <- pmin(pmax(fdr, 1e-300), 1)
  sign(nes) * qnorm(1 - fdr / 2)   # 用全局FDR作p的保守代理（门槛本身即FDR门槛）
}
# 可估计子集 = 3 种胶质 × 7 家族（GSE219280 有 NES 的轴）
glia <- three[!is.na(GSE219280_NES) & !is.na(Gittings_NES) & !is.na(Ruf_NES)]
glia[, `:=`(z_disc = zfrom(GSE219280_NES, GSE219280_global_FDR),
            z_git  = zfrom(Gittings_NES, Gittings_global_FDR),
            z_ruf  = zfrom(Ruf_NES, Ruf_global_FDR))]
# 队列有效供体数（正文口径）：发现 17、Gittings 37、Ruf 按对比/细胞类型的 n_donors
ruf_n <- ruf[contrast == "ALS_vs_Control", .(cell_type, n = unique(n_donors))]
n_ruf_map <- setNames(ruf_n$n, ruf_n$cell_type)
glia[, n_ruf := n_ruf_map[cell_type]]
glia[, `:=`(n_disc = 17, n_git = 37)]

# 情景A：真实 z 在三队列均等于 Ruf 观测 z（最乐观；对负向效应取 |z|）
glia[, powerA := pnorm(abs(z_ruf) - 1.96)^3]
# 情景B：标准化效应恒定，z 按 √(n/n_Ruf) 缩放（更现实）
glia[, powerB := pnorm(abs(z_ruf) * sqrt(n_disc / n_ruf) - 1.96) *
                  pnorm(abs(z_ruf) * sqrt(n_git / n_ruf) - 1.96) *
                  pnorm(abs(z_ruf) - 1.96)]

# Ruf 最终通过 5 门槛的 ALS 轴（可估计子集中的 4 条）
s <- fread(file.path(rf, "RNA_7_family_SVA12_sensitivity.tsv"))
tf <- function(v) !is.na(v) & v
s[, pass_all := (!is.na(global_fixed_family_FDR) & global_fixed_family_FDR < 0.05) & tf(LOO_stable_90pct) &
  tf(technical_direction_stable) & tf(age_direction_stable) & tf(SVA12_direction_stable) &
  (!is.na(SVA12_global_fixed_family_FDR) & SVA12_global_fixed_family_FDR < 0.05)]
final_als <- s[contrast == "ALS_vs_Control" & pass_all == TRUE, .(cell_type, family_id)]
glia[, in_final10 := paste(cell_type, family_id) %in% paste(final_als$cell_type, final_als$family_id)]

summA <- glia[, .(median_power = round(median(powerA), 3), max_power = round(max(powerA), 3),
                  axes_with_power50pct = sum(powerA > 0.5), axes_with_power80pct = sum(powerA > 0.8))]
summB <- glia[, .(median_power = round(median(powerB), 3), max_power = round(max(powerB), 3),
                  axes_with_power50pct = sum(powerB > 0.5), axes_with_power80pct = sum(powerB > 0.8))]
setorder(glia, -powerB)
fwrite(glia[, .(cell_type_label, family_id, z_disc, z_git, z_ruf, n_disc, n_git, n_ruf,
                powerA = round(powerA, 3), powerB = round(powerB, 3), in_final10)],
       file.path(out, "three_cohort_gate_attainability.tsv"), sep = "\t")

## ---------- Part B：GSE212630 t 分布 CI ----------
g6 <- file.path(root, "results/GSE212630_multiome")
trend <- fread(file.path(g6, "GSE212630_RNA_celltype_donor_trend.tsv"))
fam <- fread(file.path(g6, "GSE212630_RNA_7_family_donor_trend.tsv"))

# 家族表自带 t 与 beta；趋势表为 F 与 log2FC。统一反推 df2 后重建 t CI
tci <- function(dt) {
  if ("beta_per_pTDP_level" %in% names(dt)) {
    setnames(dt, c("beta_per_pTDP_level", "t"), c("beta", "t_stat_raw"))
    dt[, Fstat := t_stat_raw^2]
  } else {
    setnames(dt, c("log2FC_per_pTDP_level", "F"), c("beta", "Fstat"))
  }
  dt <- dt[!is.na(Fstat) & Fstat > 0 & !is.na(beta) & beta != 0]
  dt[, df2 := mapply(function(f, p) {
    tryCatch(uniroot(function(d) pf(f, 1, d, lower.tail = FALSE) - p, lower = 1, upper = 500, tol = 1e-8)$root,
             error = function(e) NA_real_)
  }, Fstat, p_value)]
  dt[, `:=`(t_stat = sqrt(Fstat), se = abs(beta) / sqrt(Fstat))]
  dt[!is.na(df2), `:=`(tcrit = qt(0.975, df2),
                       ci_low_t = beta - qt(0.975, df2) * se,
                       ci_high_t = beta + qt(0.975, df2) * se)]
  dt[, `:=`(ci_low_norm = beta - 1.96 * se, ci_high_norm = beta + 1.96 * se)]
  dt[, log2FC_per_pTDP_level := beta][]
}
tc_trend <- tci(trend)
tc_fam <- tci(fam)
fwrite(tc_trend[, .(feature, cell_type, log2FC = log2FC_per_pTDP_level, F, p_value, FDR, n_donors, df2,
                    se, ci_low_norm, ci_high_norm, ci_low_t, ci_high_t)],
       file.path(out, "GSE212630_trend_t_based_CI.tsv"), sep = "\t")
fwrite(tc_fam[, .(family, cell_type, log2FC = log2FC_per_pTDP_level, F, p_value, global_FDR, n_donors, df2,
                  se, ci_low_norm, ci_high_norm, ci_low_t, ci_high_t)],
       file.path(out, "GSE212630_family_t_based_CI.tsv"), sep = "\t")

cat("=== A. 三队列门槛可达性（21 条胶质 ALS 轴） ===\n")
cat("情景A(Ruf效应为真, 最乐观): 中位功效", summA$median_power, " 最大", summA$max_power,
    " | >50%功效的轴:", summA$axes_with_power50pct, "/21 | >80%:", summA$axes_with_power80pct, "/21\n")
cat("情景B(效应恒定,√n缩放):      中位功效", summB$median_power, " 最大", summB$max_power,
    " | >50%:", summB$axes_with_power50pct, "/21 | >80%:", summB$axes_with_power80pct, "/21\n")
cat("\n可估计且在 Ruf 最终 10 条中的轴（4条）:\n")
print(glia[in_final10 == TRUE, .(cell_type_label, family_id, z_ruf = round(z_ruf, 2),
                                 powerA = round(powerA, 3), powerB = round(powerB, 3),
                                 GSE_disc = GSE219280_NES, Git = Gittings_NES, Ruf = Ruf_NES)])
cat("\n=== B. 关键结果 t 分布 CI（正态近似 → t） ===\n")
key_fam <- tc_fam[grepl("Protein stabilization", family) & cell_type == "ASC"]
cat("星形胶质细胞蛋白稳定家族：\n")
print(key_fam[, .(family, cell_type, log2FC = round(log2FC_per_pTDP_level, 4), df2 = round(df2, 1),
                  ci_norm = paste0("[", round(ci_low_norm, 4), ", ", round(ci_high_norm, 4), "]"),
                  ci_t = paste0("[", round(ci_low_t, 4), ", ", round(ci_high_t, 4), "]"))])
key_b <- tc_trend[feature == "BAG3" & cell_type == "MG"]
cat("小胶质 BAG3：\n")
print(key_b[, .(feature, cell_type, log2FC = round(log2FC_per_pTDP_level, 4), df2 = round(df2, 1),
                ci_norm = paste0("[", round(ci_low_norm, 4), ", ", round(ci_high_norm, 4), "]"),
                ci_t = paste0("[", round(ci_low_t, 4), ", ", round(ci_high_t, 4), "]"))])
cat("\nCI_t 覆盖检查（FDR<0.05 的行中，t-CI 是否仍不含 0）：",
    tc_trend[FDR < 0.05 & !is.na(ci_low_t), sum(sign(ci_low_t) == sign(ci_high_t))], "/",
    tc_trend[FDR < 0.05 & !is.na(ci_low_t), .N], "行不含0\n")
