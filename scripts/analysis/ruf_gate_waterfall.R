#!/usr/bin/env Rscript
# 统计修复②：Ruf 五门槛瀑布分解 + 三队列零迁移归因
# 目的：回答审稿问题——126→63→40 的每一步淘汰发生在哪个门槛；三队列零结果是方向翻转还是单门槛功效损失。
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
root <- analysis_project_root()
configure_project_library(root)
.libPaths(c(file.path(root, "cache/R_library"), .libPaths()))
suppressPackageStartupMessages(library(data.table))

rf <- file.path(root, "results/Ruf_2026_frozen_validation")
s <- fread(file.path(rf, "RNA_7_family_SVA12_sensitivity.tsv"))
three <- fread(file.path(rf, "three_cohort_ALS_programme_replication.tsv"))
out <- file.path(root, "results/statistical_robustness_v1")
stopifnot(nrow(s) == 126)

# 五门槛（与正文一致，SVA12 表已并入前四门槛结果；NA 一律视为未通过）
tf <- function(v) !is.na(v) & v
s[, g1_global_FDR := !is.na(global_fixed_family_FDR) & global_fixed_family_FDR < 0.05]
s[, g2_LOO := tf(LOO_stable_90pct)]
s[, g3_tech_age := tf(technical_direction_stable) & tf(age_direction_stable)]
s[, g4_SVA12_dir := tf(SVA12_direction_stable)]
s[, g5_SVA12_FDR := !is.na(SVA12_global_fixed_family_FDR) & SVA12_global_fixed_family_FDR < 0.05]
s[, pass_all := g1_global_FDR & g2_LOO & g3_tech_age & g4_SVA12_dir & g5_SVA12_FDR]
s[, strict_base := tf(strict_external_programme_replication)]

# 每条轴的首个失败门槛（顺序：G1全局FDR → G2 LOO → G3技术+年龄 → G4 SVA12方向 → G5 SVA12 FDR）
s[, first_failed_gate := fifelse(!g1_global_FDR, "G1_全局家族FDR",
  fifelse(!g2_LOO, "G2_供体LOO稳定",
  fifelse(!g3_tech_age, "G3_技术+年龄方向",
  fifelse(!g4_SVA12_dir, "G4_SVA12方向",
  fifelse(!g5_SVA12_FDR, "G5_SVA12全局FDR", "全部通过")))))]

# 瀑布表：逐门槛累计
wf <- data.table(
  gate = c("全部轴", "G1_全局家族FDR<0.05", "G2_供体LOO≥90%(累计)", "G3_技术+年龄同向(累计)",
           "G4_SVA12同向(累计)", "G5_SVA12全局FDR<0.05(累计)", "strict_base字段对照"),
  pass_n = c(nrow(s),
    s[, sum(g1_global_FDR)],
    s[g1_global_FDR == TRUE, sum(g2_LOO)],
    s[g1_global_FDR == TRUE & g2_LOO == TRUE, sum(g3_tech_age)],
    s[g1_global_FDR == TRUE & g2_LOO == TRUE & g3_tech_age == TRUE, sum(g4_SVA12_dir)],
    s[, sum(pass_all)],
    s[, sum(strict_base)]),
  drop_n = c(NA, NA, NA, NA, NA, NA, NA))
wf[, cumulative_fraction := round(pass_n / nrow(s), 3)]
wf[1, drop_n := NA]
fwrite(wf, file.path(out, "ruf_gate_waterfall.tsv"), sep = "\t")

# 首失败门槛归因（126 条全轴 + ALS 对对照子集）
attr_all <- s[, .N, by = first_failed_gate][order(-N)]
setnames(attr_all, c("first_failed_gate", "n_axes"))
attr_all[, subset := "全部126轴"]
attr_als <- s[contrast == "ALS_vs_Control", .N, by = first_failed_gate][order(-N)]
setnames(attr_als, c("first_failed_gate", "n_axes"))
attr_als[, subset := "ALS_vs_Control 21轴"]
fwrite(rbind(attr_all, attr_als), file.path(out, "ruf_gate_first_failure_attribution.tsv"), sep = "\t")

# 单门槛边际通过率（不分先后，诊断哪个门槛最致命）
marg <- rbindlist(list(
  data.table(gate = "G1_全局家族FDR", pass = s[, sum(g1_global_FDR)]),
  data.table(gate = "G2_供体LOO≥90%", pass = s[, sum(g2_LOO)]),
  data.table(gate = "G3_技术+年龄同向", pass = s[, sum(g3_tech_age)]),
  data.table(gate = "G4_SVA12方向", pass = s[, sum(g4_SVA12_dir)]),
  data.table(gate = "G5_SVA12全局FDR", pass = s[, sum(g5_SVA12_FDR)])))
marg[, pass_fraction := round(pass / 126, 3)]
fwrite(marg, file.path(out, "ruf_gate_marginal_pass.tsv"), sep = "\t")

# 三队列零迁移归因：ALS 对对照且 Ruf 通过全部 5 门槛的 10 条轴，在另两队列死在哪
# 注意：GSE219280 机制矩阵仅覆盖 3 种胶质；神经元/小胶质轴在三队列门槛上不可估计（NA）
als_final <- s[contrast == "ALS_vs_Control" & pass_all == TRUE,
  .(cell_type, family_id, Ruf_NES = NES, Ruf_global_FDR = global_fixed_family_FDR)]
tc <- merge(three, als_final, by = c("cell_type", "family_id"), suffixes = c("", ".final"), all.x = FALSE)
tc[, GSE219280_estimable := !is.na(GSE219280_NES)]
tc[, GSE219280_FDR_ok := GSE219280_estimable & GSE219280_global_FDR < 0.05]
tc[, Gittings_FDR_ok := !is.na(Gittings_global_FDR) & Gittings_global_FDR < 0.05]
tc[, dir_GSE219280_ok := !GSE219280_estimable | (Ruf_vs_GSE219280_direction_concordant == TRUE)]
tc[, dir_Gittings_ok := Ruf_vs_Gittings_direction_concordant == TRUE]
tc[, failure_cause := fifelse(!GSE219280_estimable, "三队列门槛不可估计(发现队列无该细胞类型矩阵)",
  fifelse(dir_GSE219280_ok & dir_Gittings_ok & GSE219280_FDR_ok & Gittings_FDR_ok, "三队列同向双显著(应为阳性)",
  fifelse(!dir_GSE219280_ok & !dir_Gittings_ok, "两队列方向均翻转",
  fifelse(!dir_GSE219280_ok, "GSE219280方向翻转",
  fifelse(!dir_Gittings_ok, "Gittings方向翻转",
  fifelse(!GSE219280_FDR_ok & !Gittings_FDR_ok, "两队列同向但FDR均不显著",
  fifelse(!GSE219280_FDR_ok, "GSE219280同向但FDR不显著", "Gittings同向但FDR不显著")))))))]
keep <- c("cell_type_label", "family_id", "Ruf_NES", "Ruf_global_FDR", "GSE219280_NES", "GSE219280_global_FDR",
          "Gittings_NES", "Gittings_global_FDR", "dir_GSE219280_ok", "GSE219280_FDR_ok",
          "dir_Gittings_ok", "Gittings_FDR_ok", "failure_cause")
fwrite(tc[, ..keep], file.path(out, "three_cohort_ALS_zero_attribution.tsv"), sep = "\t")

# 附加：全部 21 条可估计胶质 ALS 轴的三队列状态汇总（不只是 Ruf 通过的 10 条）
tc_all <- three[!is.na(GSE219280_NES),
  .(cell_type_label, family_id, Ruf_strict_internal_gate, Ruf_global_FDR,
    GSE219280_global_FDR, Gittings_global_FDR,
    dir_all3 = (strict_three_cohort_same_direction == TRUE))]
cat("\n=== 21条胶质ALS轴三队列同向标记 ===\n"); print(tc_all[, .N, by = dir_all3])

cat("=== 瀑布 ===\n"); print(wf[, .(gate, pass_n, cumulative_fraction)])
cat("\n=== 首失败门槛归因（全部126） ===\n"); print(attr_all)
cat("\n=== ALS 子集首失败 ===\n"); print(attr_als)
cat("\n=== 单门槛边际通过 ===\n"); print(marg)
cat("\n=== 三队列零迁移归因（Ruf 最终 10 条 ALS 轴） ===\n")
print(tc[, .N, by = failure_cause])
cat("\nRuf 最终通过(pass_all) =", s[, sum(pass_all)], "；ALS 对对照 =", s[contrast == "ALS_vs_Control", sum(pass_all)], "\n")
