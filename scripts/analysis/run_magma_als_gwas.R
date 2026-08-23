#!/usr/bin/env Rscript
# GWAS 锚定：MAGMA 基因层关联 + 程序基因集竞争性富集
# 输入：ALS GWAS GCST90027164（van Rheenen 2021，29,612/122,656）+ g1000_eur 面板 + GENCODE v19
# 输出：results/GWAS_genetic_anchoring/
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
root <- analysis_project_root()
configure_project_library(root)
.libPaths(c(file.path(root, "cache/R_library"), .libPaths()))
suppressPackageStartupMessages(library(data.table))
mag <- file.path(root, "external/magma")
out <- file.path(root, "results/GWAS_genetic_anchoring")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
magma_bin <- file.path(mag, "magma")

## 1) gene.loc（GENCODE v19 → MAGMA 格式：symbol chr start end strand）
gtf <- fread(cmd = paste("gunzip -c", shQuote(file.path(mag, "gencode.v19.annotation.gtf.gz")),
  "| awk '$3==\"gene\"'"), sep = "\t", header = FALSE, select = c(1, 4, 5, 7, 9), quote = "")
setnames(gtf, c("chr", "start", "end", "strand", "attr"))
idm <- gtf$attr
gtf[, symbol := gsub('gene_name "|"', "", regmatches(idm, regexpr('gene_name "[^"]+"', idm)))]
gtf[, biotype := gsub('gene_type "|"', "", regmatches(idm, regexpr('gene_type "[^"]+"', idm)))]
gtf[, chr := sub("^chr", "", chr)]
gl <- unique(gtf[grepl("protein_coding", biotype) & chr %in% c(1:22, "X") & !is.na(symbol) & nzchar(symbol),
  .(symbol, chr, start, end, strand)])
gl <- gl[!duplicated(symbol)]
fwrite(gl, file.path(mag, "gencode_v19_gene.loc"), sep = "\t", quote = FALSE)
cat("gene.loc:", nrow(gl), "protein-coding genes\n")

## 2) 注释（±10 kb 窗口；已存在则跳过）
if (!file.exists(file.path(mag, "als_annot.genes.annot"))) {
  system2(magma_bin, c("--annotate", "window=10",
    "--snp-loc", shQuote(file.path(mag, "g1000_eur.bim")),
    "--gene-loc", shQuote(file.path(mag, "gencode_v19_gene.loc")),
    "--out", shQuote(file.path(mag, "als_annot"))))
}

## 3) 基因层关联（SNP-wise 主模型，per-SNP N；每个修饰符单独传）
if (!file.exists(file.path(mag, "als_gene.genes.out"))) {
  system2(magma_bin, c("--bfile", shQuote(file.path(mag, "g1000_eur")),
    "--pval", shQuote(file.path(mag, "ALS_GCST90027164.pval.tab")),
    "ncol=N",
    "--gene-annot", shQuote(file.path(mag, "als_annot.genes.annot")),
    "--out", shQuote(file.path(mag, "als_gene"))))
}

genes <- fread(file.path(mag, "als_gene.genes.out"))
setorder(genes, P)
fwrite(genes, file.path(out, "ALS_GWAS_gene_level_P.tsv"), sep = "\t")
cat("\n=== ALS GWAS 基因层 Top20 ===\n"); print(head(genes, 20))

## 4) 程序基因集竞争性富集（10 个预冻结基因集）
system2(magma_bin, c("--gene-results", shQuote(file.path(mag, "als_gene.genes.raw")),
  "--set-annot", shQuote(file.path(mag, "program_gene_sets.gsa")),
  "--out", shQuote(file.path(mag, "als_sets"))))
sets <- fread(file.path(mag, "als_sets.gsa.out"), skip = "VARIABLE")
setorder(sets, P)
fwrite(sets, file.path(out, "ALS_GWAS_program_gene_set_enrichment.tsv"), sep = "\t")
cat("\n=== 基因集竞争性富集 ===\n"); print(sets[])

## 5) 17 候选基因的基因层 P 值直接提取
cand17 <- c("BAG3","HSP90AA1","HSP90AB1","CHORDC1","ARHGAP35","PARK7","CCT5","ARHGAP39",
  "COX5B","COX7C","MAP1B","MAPT","BCR","CCT8","COA3","DOCK7","STARD13")
cand_tab <- genes[GENE %in% cand17][order(P)]
fwrite(cand_tab, file.path(out, "ALS_GWAS_17_candidates_gene_P.tsv"), sep = "\t")
cat("\n=== 17 冻结候选基因层关联 ===\n"); print(cand_tab)
