#!/usr/bin/env Rscript
# 统计修复①：跨队列程序轴方向一致性的精确二项检验
# 输入：results/submission_cross_cohort_evidence_v1/exact_estimand_program_replication.tsv
# 输出：results/statistical_robustness_v1/
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
root <- analysis_project_root()
configure_project_library(root)
.libPaths(c(file.path(root, "cache/R_library"), .libPaths()))
suppressPackageStartupMessages(library(data.table))

inp <- file.path(root, "results/submission_cross_cohort_evidence_v1/exact_estimand_program_replication.tsv")
out <- file.path(root, "results/statistical_robustness_v1")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- fread(inp)
stopifnot(nrow(d) == 63)
d[, same_direction := same_direction == "TRUE" | same_direction == TRUE]

# 效应量为 NES 尺度的程序得分；两队列 NES 直接 Spearman（支持性分析）
sp <- cor.test(d$effect_discovery, d$effect_replication, method = "spearman", exact = FALSE)

btest <- function(k, n, label, stratum) {
  if (n == 0) return(data.table(stratum = stratum, label = label, n = 0L, k = 0L,
    fraction = NA_real_, p_two_sided = NA_real_, ci95_low = NA_real_, ci95_high = NA_real_))
  bt <- binom.test(k, n, p = 0.5, alternative = "two.sided", conf.level = 0.95)
  data.table(stratum = stratum, label = label, n = n, k = k,
    fraction = k / n, p_two_sided = unname(bt$p.value),
    ci95_low = bt$conf.int[1], ci95_high = bt$conf.int[2])
}

res <- list()
res$overall <- btest(sum(d$same_direction), nrow(d), "全部63条匹配轴", "overall")

# 分层：机制家族（7，各9条）、estimand（3，各21）、细胞类型（3，各21）
res$family <- d[, btest(sum(same_direction), .N, paste(unique(family_label)), "family"), by = family_id]
res$family[, family_id := NULL]
res$estimand <- d[, btest(sum(same_direction), .N, unique(estimand), "estimand"), by = estimand]$V1
res$cell_type <- d[, btest(sum(same_direction), .N, unique(cell_type), "cell_type"), by = cell_type]$V1

# 预设合并层：三个高一致家族（轴突发育/小GTP酶/突触连接）与线粒体家族对照
hi <- d[family_id %in% c("axon_development", "small_GTPase", "synapse_junction")]
res$pool3 <- btest(sum(hi$same_direction), nrow(hi), "轴突+小GTP酶+突触合并(预设高一致层)", "pooled")
mito <- d[family_id == "mitochondrial_energy"]
res$mito <- btest(sum(mito$same_direction), nrow(mito), "线粒体能量家族(方向翻转层)", "pooled")

# 敏感性1：仅发现队列全局显著的轴（方向更可信的子集）
d1 <- d[adjusted_p_discovery < 0.05]
res$discovery_sig <- btest(sum(d1$same_direction), nrow(d1), "发现队列FDR<0.05子集", "sensitivity")

# 敏感性2：任一队列全局显著的轴
d2 <- d[adjusted_p_discovery < 0.05 | adjusted_p_replication < 0.05]
res$either_sig <- btest(sum(d2$same_direction), nrow(d2), "任一队列FDR<0.05子集", "sensitivity")

tab <- rbindlist(res, use.names = TRUE)
tab[, fraction := round(fraction, 4)]
tab[, p_two_sided := signif(p_two_sided, 4)]
tab[, `:=`(ci95_low = round(ci95_low, 4), ci95_high = round(ci95_high, 4))]
setorder(tab, stratum, label)
fwrite(tab, file.path(out, "binomial_concordance_results.tsv"), sep = "\t")

# 严格复制轴的一致性下限（13条双显著轴中同向数，供正文引用）
strict <- d[replication_class == "bidirectional_global_FDR_concordant"]
summ <- data.table(
  metric = c("matched_axes", "same_direction_n", "exact_binomial_p",
             "cp95_low", "cp95_high", "spearman_NES_rho", "spearman_p",
             "strict_both_FDR_axes", "strict_same_direction_n"),
  value = c(nrow(d), sum(d$same_direction), tab[stratum == "overall"]$p_two_sided,
            tab[stratum == "overall"]$ci95_low, tab[stratum == "overall"]$ci95_high,
            round(unname(sp$estimate), 4), signif(sp$p.value, 4),
            nrow(strict), sum(strict$same_direction)))
fwrite(summ, file.path(out, "binomial_concordance_summary.tsv"), sep = "\t")

cat("=== 总体 ===\n"); print(tab[stratum == "overall"])
cat("\n=== 分机制家族 ===\n"); print(res$family)
cat("\n=== Spearman(NES discovery vs replication, n=63) rho =",
    round(unname(sp$estimate), 4), " p =", signif(sp$p.value, 4), "===\n")
cat("\n输出：", out, "/binomial_concordance_results.tsv\n", sep = "")
