#!/usr/bin/env Rscript
# TDP-43 扰动干验证：人 iPSC 皮层神经元 TDP-43 敲低（GSE296710 Zanovello 4v4；GSE296714 Humphrey 6v6）
# 检验：TDP-43 缺失本身是否使 7 机制家族 / 17 冻结候选 / BAG3–HSP90 蛋白稳态基因
# 沿与 GSE212630 pTDP 趋势一致的方向变化（程序层同向性 + 基因层相关）
# 注意：扰动在神经元、疾病趋势在胶质——这是"机制一致性检验"，非细胞类型复制。
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
root <- analysis_project_root()
configure_project_library(root)
.libPaths(c(file.path(root, "cache/R_library"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(edgeR)})

ext <- file.path(root, "external")
out <- file.path(root, "results/TDP43_perturbation_validation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

## 1) ENSG→symbol（GENCODE v19，GRCh37）
gtf <- fread(cmd = paste("gunzip -c", shQuote(file.path(ext, "magma/gencode.v19.annotation.gtf.gz")),
  "| awk '$3==\"gene\"'"), sep = "\t", header = FALSE, select = 9, quote = "")
idm <- gtf[[1]]
ensg <- sub('^(ENSG[0-9]+)\\..*$', '\\1', regmatches(idm, regexpr('gene_id "[^"]+"', idm)))
symb <- regmatches(idm, regexpr('gene_name "[^"]+"', idm))
map_dt <- data.table(ensg = gsub('gene_id "|"', "", ensg),
                     symbol = gsub('gene_name "|"', "", symb))
map_dt[, ensg := sub("\\..*$", "", ensg)]
map_dt <- map_dt[!duplicated(ensg)]
cat("map_dt:", nrow(map_dt), "genes; TARDBP present:", "TARDBP" %in% map_dt$symbol, "\n")

## 2) 两个扰动数据集的 KD vs 对照 DE
run_de <- function(path, ctrl_pat, kd_pat, label) {
  m <- as.data.frame(fread(cmd = paste("gunzip -c", shQuote(file.path(ext, path)))))
  rownames(m) <- sub("\\..*$", "", m$gene_id); m$gene_id <- NULL
  m <- m[!duplicated(rownames(m)), ]
  ctrl <- grep(ctrl_pat, colnames(m), value = TRUE); kd <- grep(kd_pat, colnames(m), value = TRUE)
  stopifnot(length(ctrl) >= 2, length(kd) >= 2)
  y <- DGEList(m[, c(ctrl, kd)])
  y <- y[rowSums(cpm(y) > 1) >= 2, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  grp <- factor(c(rep("ctrl", length(ctrl)), rep("kd", length(kd))), levels = c("ctrl", "kd"))
  fit <- glmQLFit(estimateDisp(y, model.matrix(~grp), robust = TRUE))
  tt <- glmQLFTest(fit, coef = 2)
  res <- as.data.table(topTags(tt, Inf)$table)
  res[, ensg := rownames(topTags(tt, Inf)$table)]
  if ("logFC" %in% names(res)) setnames(res, "logFC", "log2FC")
  res[, dataset := label]
  out_dt <- merge(res, map_dt, by = "ensg", all.x = TRUE)[, .(dataset, ensg, symbol, log2FC, logCPM,
    p = PValue, FDR)][order(p)]
  cat("  [", label, "] DE rows:", nrow(out_dt), "| symbol NA:", sum(is.na(out_dt$symbol)),
      "| TARDBP row:", paste(out_dt[symbol == "TARDBP", round(log2FC, 3)], collapse=","), "\n")
  out_dt
}
z <- run_de("TDP43_perturbation_GSE296710_counts.csv.gz", "^i3Neuron_CTRL_CTRL", "^i3Neuron_TDP43_CTRL", "GSE296710_zanovello")
h <- run_de("TDP43_perturbation_GSE296714_counts.csv.gz", "^CTRL", "^TDP", "GSE296714_humphrey")
z <- z[!is.na(symbol)]; h <- h[!is.na(symbol)]
de <- rbind(z, h)
fwrite(de, file.path(out, "TDP43_KD_vs_ctrl_DE.tsv"), sep = "\t")

# sanity：TARDBP 自身应显著下调
cat("=== sanity: TARDBP ===\n"); print(de[symbol == "TARDBP", .(dataset, symbol, log2FC = round(log2FC, 3), FDR = signif(FDR, 3))])

## 3) 两数据集 KD 签名一致性（内部信度）
both <- merge(z[, .(symbol, z = log2FC)], h[, .(symbol, h = log2FC)], by = "symbol")
sp <- cor.test(both$z, both$h, method = "spearman", exact = FALSE)
sgn <- both[abs(z) > 0.5 & abs(h) > 0.5]
bt <- binom.test(sum(sign(sgn$z) == sign(sgn$h)), nrow(sgn), 0.5)
cat("\n=== 两 KD 数据集一致性 ===\nSpearman rho =", round(unname(sp$estimate), 3), " p =", signif(sp$p.value, 3),
    "| |log2FC|>0.5 同向", sum(sign(sgn$z) == sign(sgn$h)), "/", nrow(sgn), " 二项 p =", signif(bt$p.value, 4), "\n")

## 4) 与疾病程序的对齐（GSE212630 pTDP 趋势）
tr <- fread(file.path(root, "results/GSE212630_multiome/GSE212630_RNA_celltype_donor_trend.tsv"))
fam <- fread(file.path(root, "results/GSE212630_multiome/GSE212630_RNA_7_family_donor_trend.tsv"))
mem <- fread(file.path(root, "results/GSE219280_same_celltype_fixed_programs/predefined_gene_set_membership.tsv"))
cand17 <- c("BAG3","HSP90AA1","HSP90AB1","CHORDC1","ARHGAP35","PARK7","CCT5","ARHGAP39",
  "COX5B","COX7C","MAP1B","MAPT","BCR","CCT8","COA3","DOCK7","STARD13")

align <- function(kd_dt, lab) {
  # 家族层：KD 中家族成员的均值 log2FC vs 疾病侧家族趋势方向
  fam_kd <- kd_dt[symbol %in% unique(mem$gene), .(mean_log2FC = mean(log2FC), n = .N, 
    n_up = sum(log2FC > 0), n_down = sum(log2FC < 0)), by = .(family_id = mem[match(symbol, gene), family_id])]
  fam_kd[, binom_p := mapply(function(k, n) binom.test(k, n, .5)$p.value, pmax(n_up, n_down), n)]
  # 基因层：17 候选与 pTDP 趋势的相关（每细胞类型）
  cor_cells <- rbindlist(lapply(unique(tr$cell_type), function(ct) {
    x <- merge(kd_dt[, .(symbol, kd = log2FC)], tr[cell_type == ct, .(symbol = feature, disease = log2FC_per_pTDP_level)], by = "symbol")
    x <- x[symbol %in% cand17]
    if (nrow(x) < 8) return(NULL)
    cc <- suppressWarnings(cor.test(x$kd, x$disease, method = "spearman", exact = FALSE))
    data.table(dataset = lab, cell_type = ct, n = nrow(x), rho = round(unname(cc$estimate), 3), p = signif(cc$p.value, 3))
  }))
  list(fam = fam_kd, cor = cor_cells)
}
a_z <- align(z, "GSE296710_zanovello"); a_h <- align(h, "GSE296714_humphrey")
fam_tab <- rbind(a_z$fam[, dataset := "GSE296710_zanovello"], a_h$fam[, dataset := "GSE296714_humphrey"])
cor_tab <- rbind(a_z$cor, a_h$cor)
fwrite(fam_tab, file.path(out, "KD_family_mean_log2FC.tsv"), sep = "\t")
fwrite(cor_tab, file.path(out, "KD_vs_pTDP_trend_spearman_17candidates.tsv"), sep = "\t")

cat("\n=== KD 家族均值 log2FC（Zanovello） ===\n"); print(a_z$fam[order(-mean_log2FC)])
cat("\n=== KD vs pTDP 趋势相关（17 候选） ===\n"); print(cor_tab)

## 5) BAG3–HSP90 蛋白稳态核心的直接检验
core <- c("BAG3","HSP90AA1","HSP90AB1","HSP90B1","CCT5","CCT8","CHORDC1","DNAJB1","HSPH1")
direct <- rbind(z[symbol %in% core, .(dataset = "GSE296710_zanovello", symbol, KD_log2FC = round(log2FC, 3), FDR = signif(FDR, 3))],
                h[symbol %in% core, .(dataset = "GSE296714_humphrey", symbol, KD_log2FC = round(log2FC, 3), FDR = signif(FDR, 3))])
fwrite(direct, file.path(out, "proteostasis_core_genes_KD_effect.tsv"), sep = "\t")
cat("\n=== 蛋白稳态核心基因 KD 效应 ===\n"); print(direct)
