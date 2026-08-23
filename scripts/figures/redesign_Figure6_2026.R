script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
project_lib <- Sys.getenv("ALS_FTD_R_LIBRARY", unset = "")
if (nzchar(project_lib) && dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(svglite)
  library(ragg)
})

source_root <- file.path(project_root, "data", "source_data", "main_figures")
figure6_inputs <- file.path(project_root, "data", "figure_inputs", "Figure6")
draft_root <- file.path(project_root, "outputs", "figures", "Figure6")
draft_source <- file.path(draft_root, "source_data")
dir.create(draft_source, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  discovery = file.path(source_root, "Figure2_panel_a_source_data.tsv"),
  replication = file.path(figure6_inputs, "Figure6_panel_c_source_data.tsv"),
  boundary = file.path(figure6_inputs, "Figure6_panel_b_source_data.tsv"),
  rna = file.path(source_root, "Figure5_panel_b_source_data.tsv"),
  dar = file.path(source_root, "Figure5_panel_c_source_data.tsv"),
  peak_gene = file.path(source_root, "Figure5_panel_d_source_data.tsv"),
  joint = file.path(source_root, "Figure5_panel_e_source_data.tsv"),
  perturbation = file.path(
    figure6_inputs, "KD_vs_pTDP_trend_spearman_17candidates.tsv"
  ),
  gwas_sets = file.path(
    figure6_inputs, "ALS_GWAS_program_gene_set_enrichment.tsv"
  ),
  gwas_genes = file.path(
    figure6_inputs, "ALS_GWAS_17_candidates_gene_P.tsv"
  )
)
if (!all(file.exists(unlist(paths)))) {
  stop("One or more frozen Figure 6 evidence inputs are missing.")
}

discovery <- fread(paths$discovery)
replication <- fread(paths$replication)
boundary <- fread(paths$boundary)
rna <- fread(paths$rna)
dar <- fread(paths$dar)
peak_gene <- fread(paths$peak_gene)
joint <- fread(paths$joint)
perturbation <- fread(paths$perturbation)
gwas_sets <- fread(paths$gwas_sets)
gwas_genes <- fread(paths$gwas_genes)

# Figure contract:
# Core conclusion: the recurrent ALS–FTD programmes are robust under donor-aware
# and external tests, but their direction remains context dependent and neither a
# universal ALS signature nor a complete RNA–ATAC causal chain is established.
# Archetype: schematic-led composite with a full-width conceptual hero panel,
# quantitative evidence cards, orthogonal anchors and an interpretation boundary.
# Backend: R only. Final output: editable SVG/PDF, 600-dpi TIFF, PNG preview and TSV data.

pal <- c(
  ink = "#172033",
  muted = "#667085",
  line = "#CCD5E0",
  soft = "#F5F7FA",
  blue = "#2867A8",
  blue_soft = "#DFEAF5",
  teal = "#278F88",
  teal_soft = "#DCEFEB",
  mauve = "#A85674",
  mauve_soft = "#F1DEE6",
  amber = "#CF872C",
  amber_soft = "#F7E9D5",
  violet = "#6955A5",
  violet_soft = "#E9E4F3",
  neutral = "#98A4B3",
  neutral_soft = "#EEF2F6"
)

font_body <- "Times New Roman"
font_display <- "Times New Roman"

write_clean_tsv <- function(x, path) {
  y <- copy(as.data.table(x))
  text_columns <- names(y)[vapply(y, function(z) is.character(z) || is.factor(z), logical(1))]
  for (column in text_columns) {
    value <- gsub("[\r\n\t]", " ", as.character(y[[column]]), perl = TRUE)
    set(y, j = column, value = value)
  }
  fwrite(y, path, sep = "\t", quote = FALSE, na = "")
}

round_box <- function(xmin, xmax, ymin, ymax, fill, colour, linewidth = 0.45, radius_mm = 2.4) {
  annotation_custom(
    grid::roundrectGrob(
      r = grid::unit(radius_mm, "mm"),
      gp = grid::gpar(fill = fill, col = colour, lwd = linewidth * 2.85)
    ),
    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax
  )
}

add_text <- function(p, x, y, label, size = 2.2, colour = pal[["ink"]],
                     family = font_body, face = "plain", hjust = 0.5,
                     vjust = 0.5, lineheight = 1.0) {
  p + annotate(
    "text", x = x, y = y, label = label, size = size, colour = colour,
    family = family, fontface = face, hjust = hjust, vjust = vjust,
    lineheight = lineheight
  )
}

add_panel_heading <- function(p, tag, title, subtitle, x, y) {
  p <- add_text(p, x - 2.5, y + 0.3, tag, size = 3.5, family = font_display,
                face = "bold", hjust = 1)
  p <- add_text(p, x, y, title, size = 3.25, family = font_display,
                face = "bold", hjust = 0)
  add_text(p, x, y - 3.0, subtitle, size = 2.05, colour = pal[["muted"]], hjust = 0)
}

# Frozen quantitative anchors.
n_discovery_total <- nrow(discovery)
n_discovery_pass <- discovery[global_primary_FDR < 0.05, .N]
n_strict_replication <- 13L
n_ruf_total <- 126L
n_ruf_robust <- 40L
n_three_cohort <- 0L
n_ptdp_donors <- max(rna$n_donors, na.rm = TRUE)
n_dars <- nrow(dar)
n_odc_dars <- dar[cell_type == "ODC", .N]
n_peak_gene <- nrow(peak_gene)
n_complete_chain <- joint[criterion == "Full three-part chain" & passed == TRUE, .N]

glial_perturbation <- perturbation[cell_type %chin% c("ASC", "ODC", "OPC")]
n_glial_opposite <- glial_perturbation[rho < 0, .N]
n_glial_tests <- nrow(glial_perturbation)
glial_binom_p <- binom.test(n_glial_opposite, n_glial_tests, p = 0.5)$p.value

n_gwas_sets <- nrow(gwas_sets)
n_gwas_set_hits <- gwas_sets[P < 0.05, .N]
bonf_threshold <- 0.05 / nrow(gwas_genes)
gwas_gene_hits <- gwas_genes[P < bonf_threshold]

stopifnot(
  n_discovery_total == 63L,
  n_discovery_pass == 35L,
  n_strict_replication == 13L,
  n_ruf_robust == 40L,
  n_three_cohort == 0L,
  n_ptdp_donors == 14L,
  n_dars == 409L,
  n_odc_dars == 406L,
  n_peak_gene == 7L,
  n_complete_chain == 0L,
  n_glial_tests == 6L,
  n_glial_opposite == 6L,
  abs(glial_binom_p - 0.03125) < 1e-10,
  n_gwas_sets == 9L,
  n_gwas_set_hits == 0L,
  nrow(gwas_gene_hits) == 1L,
  gwas_gene_hits$GENE == "MAPT"
)

# Panel source tables represent exactly what is displayed in the redesigned figure.
panel_a_source <- data.table(
  element = c(
    "Pathological context", "Conditioned programme model",
    "Mitochondrial/proteostasis family", "Small-GTPase/synaptic/axonal family",
    "Diagnostic observation", "Pathology-burden observation"
  ),
  quantitative_anchor = c(
    "C9orf72/TDP-43", "cell, region and cohort dependent",
    "recurrent across C9 cohorts", "contextual across external testing",
    "13 strict two-cohort axes; 40/126 robust external-test axes; 0 complete three-cohort ALS axes",
    "14 donors; 409 DARs (406 ODC); 7 peak-gene associations; 0 complete chains"
  ),
  interpretation = c(
    "pathological context", "core supported model", "recurrent programme family",
    "recurrent programme family", "not a universal signature",
    "bounded chromatin context"
  )
)

panel_b_source <- data.table(
  challenge = c(
    "Donor-aware discovery", "Two-cohort C9 replication",
    "Frozen diagnostic portability", "Paired multiome joint gate"
  ),
  numerator = c(n_discovery_pass, n_strict_replication, n_ruf_robust, n_complete_chain),
  denominator = c(n_discovery_total, NA_integer_, n_ruf_total, NA_integer_),
  displayed_metric = c("35/63", "13", "40/126", "0"),
  detail = c(
    "global-FDR programme axes",
    "dual-FDR + donor-LOO axes",
    "robust axes in 79 donors",
    "complete RNA–ATAC chains"
  ),
  secondary_detail = c(
    "biological donors: 16–17",
    "same estimand across two C9 cohorts",
    "zero complete three-cohort ALS axes",
    "14 donors; 409 DARs; 7 peak–gene associations"
  )
)

panel_c_source <- data.table(
  anchor = c("TDP-43 perturbation", "ALS GWAS gene sets", "ALS GWAS candidate genes"),
  displayed_metric = c("6/6", "0/9", "1/17"),
  detail = c(
    sprintf("glial directions opposite; exact binomial P = %.3f", glial_binom_p),
    "prespecified programme sets enriched at P < 0.05",
    sprintf("MAPT passes Bonferroni; P = %.1e", gwas_gene_hits$P)
  ),
  interpretation = c(
    "cross-context perturbation support",
    "programmes are not enriched for inherited ALS risk",
    "limited gene-level genetic anchoring"
  )
)

panel_d_source <- data.table(
  category = rep(c("Supported", "Not established", "Next validation"), each = 4),
  item = c(
    "Recurrent programme families",
    "13 strict C9 replication axes",
    "40 robust external-test axes",
    "Bounded chromatin context",
    "Universal ALS direction (0 axes)",
    "Complete RNA–ATAC chain (0)",
    "Enhancer-to-gene causality",
    "Motif / TF-activity mechanism",
    "Genetically matched C9 cohort",
    "Functional BAG3–HSP90 perturbation",
    "Cell-type enhancer validation",
    "Junction-resolved RNA analysis"
  )
)

write_clean_tsv(panel_a_source, file.path(draft_source, "Figure6_panel_a_conceptual_model.tsv"))
write_clean_tsv(panel_b_source, file.path(draft_source, "Figure6_panel_b_evidence_challenges.tsv"))
write_clean_tsv(panel_c_source, file.path(draft_source, "Figure6_panel_c_orthogonal_anchors.tsv"))
write_clean_tsv(panel_d_source, file.path(draft_source, "Figure6_panel_d_boundaries.tsv"))

# Single vector canvas avoids panel-alignment drift and keeps all objects editable.
p <- ggplot() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 112), expand = FALSE, clip = "off") +
  theme_void(base_family = font_body) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(8, 10, 8, 10, unit = "pt")
  )

# Panel a — hero conceptual model with quantitative anchors.
p <- add_panel_heading(
  p, "a", "Evidence-constrained model of ALS–FTD programmes",
  "Recurrent biological systems survive multiple tests, but their direction remains conditional on cell, region and cohort",
  4, 109
)

p <- p +
  round_box(3.5, 19.5, 88.0, 99.0, pal[["violet_soft"]], pal[["violet"]]) +
  round_box(26.0, 48.5, 85.5, 101.5, pal[["neutral_soft"]], pal[["neutral"]], linewidth = 0.55) +
  round_box(55.0, 72.5, 93.5, 102.5, pal[["mauve_soft"]], pal[["mauve"]]) +
  round_box(55.0, 72.5, 82.0, 91.0, pal[["blue_soft"]], pal[["blue"]]) +
  round_box(79.0, 97.0, 93.5, 102.5, pal[["teal_soft"]], pal[["teal"]]) +
  round_box(79.0, 97.0, 80.0, 92.0, pal[["teal_soft"]], pal[["teal"]]) +
  annotate("segment", x = 19.5, xend = 26.0, y = 93.5, yend = 93.5,
           linetype = 2, linewidth = 0.6, colour = pal[["neutral"]]) +
  annotate("segment", x = 48.5, xend = 55.0, y = 96.5, yend = 98.0,
           linetype = 2, linewidth = 0.6, colour = pal[["neutral"]]) +
  annotate("segment", x = 48.5, xend = 55.0, y = 90.5, yend = 86.5,
           linetype = 2, linewidth = 0.6, colour = pal[["neutral"]]) +
  annotate("segment", x = 72.5, xend = 79.0, y = 98.0, yend = 98.0,
           linetype = 2, linewidth = 0.6, colour = pal[["neutral"]]) +
  annotate("segment", x = 72.5, xend = 79.0, y = 86.5, yend = 86.0,
           linetype = 2, linewidth = 0.6, colour = pal[["neutral"]])

p <- add_text(p, 11.5, 96.0, "PATHOLOGICAL\nCONTEXT", 1.75, pal[["violet"]], font_display, "bold")
p <- add_text(p, 11.5, 91.5, "C9orf72 / TDP-43", 2.35, pal[["ink"]], font_display, "bold")

p <- add_text(p, 37.25, 98.2, "CORE MODEL", 1.75, pal[["muted"]], font_display, "bold")
p <- add_text(p, 37.25, 93.9, "CELL- & REGION-CONDITIONED\nPROGRAMMES", 2.55, pal[["ink"]], font_display, "bold", lineheight = 0.95)
p <- add_text(p, 37.25, 88.7, "Direction varies across cell type,\nregion and cohort", 1.85, pal[["muted"]], lineheight = 1.0)

p <- add_text(p, 63.75, 99.9, "MITOCHONDRIAL /\nPROTEOSTASIS", 2.0, pal[["mauve"]], font_display, "bold", lineheight = 0.95)
p <- add_text(p, 63.75, 95.3, "Recurrent across C9 cohorts", 1.70, pal[["ink"]])
p <- add_text(p, 63.75, 88.4, "SMALL-GTPASE / SYNAPTIC /\nAXONAL", 1.95, pal[["blue"]], font_display, "bold", lineheight = 0.95)
p <- add_text(p, 63.75, 83.9, "Contextual across external testing", 1.62, pal[["ink"]])

p <- add_text(p, 88.0, 100.0, "DIAGNOSTIC CONTRASTS", 1.85, pal[["teal"]], font_display, "bold")
p <- add_text(p, 88.0, 97.3, "13 strict axes · 40/126 robust", 1.65, pal[["ink"]])
p <- add_text(p, 88.0, 94.7, "0 universal three-cohort ALS axes", 1.60, pal[["amber"]], face = "bold")
p <- add_text(p, 88.0, 89.4, "ORDINAL pTDP BURDEN", 1.85, pal[["teal"]], font_display, "bold")
p <- add_text(p, 88.0, 86.4, "14 donors · 409 DARs (406 ODC)", 1.65, pal[["ink"]])
p <- add_text(p, 88.0, 83.6, "7 peak–gene links · 0 complete chains", 1.65, pal[["amber"]], face = "bold")
p <- add_text(p, 50, 77.9, "Dashed links denote evidence relationships — no causal direction is implied", 1.75, pal[["muted"]], hjust = 0.5)

# Section separators.
p <- p +
  annotate("segment", x = 2, xend = 98, y = 75.3, yend = 75.3,
           linewidth = 0.35, colour = pal[["line"]]) +
  annotate("segment", x = 63.5, xend = 63.5, y = 40.5, yend = 73.3,
           linewidth = 0.35, colour = pal[["line"]]) +
  annotate("segment", x = 2, xend = 98, y = 38.2, yend = 38.2,
           linewidth = 0.35, colour = pal[["line"]])

# Panel b — quantitative challenge cards.
p <- add_panel_heading(
  p, "b", "What remains after each analytical challenge",
  "Counts come from distinct testing universes and are not a sequential funnel",
  4, 72.1
)

b_cards <- data.table(
  xmin = c(4.0, 18.7, 33.4, 48.1),
  xmax = c(17.7, 32.4, 47.1, 61.8),
  fill = c(pal[["blue_soft"]], pal[["teal_soft"]], pal[["teal_soft"]], pal[["amber_soft"]]),
  border = c(pal[["blue"]], pal[["teal"]], pal[["teal"]], pal[["amber"]]),
  metric = panel_b_source$displayed_metric,
  label = c("DISCOVERY", "C9 REPLICATION", "FROZEN TEST", "JOINT GATE"),
  detail = c(
    "global-FDR\nprogramme axes",
    "dual-FDR + donor-LOO\naxes",
    "robust axes\nin 79 donors",
    "complete RNA–ATAC\nchains"
  ),
  foot = c(
    "16–17 donors", "two C9 cohorts",
    "0 three-cohort ALS axes", "14 donors · 409 DARs · 7 links"
  )
)

for (i in seq_len(nrow(b_cards))) {
  row <- b_cards[i]
  p <- p + round_box(row$xmin, row$xmax, 44.0, 65.6, row$fill, row$border)
  p <- add_text(p, (row$xmin + row$xmax) / 2, 62.5, row$label, 1.62, row$border, font_display, "bold")
  p <- add_text(p, (row$xmin + row$xmax) / 2, 57.3, row$metric, 4.35, row$border, font_display, "bold")
  p <- add_text(p, (row$xmin + row$xmax) / 2, 51.8, row$detail, 1.72, pal[["ink"]], lineheight = 0.95)
  p <- add_text(p, (row$xmin + row$xmax) / 2, 46.3, row$foot, 1.48, pal[["muted"]], lineheight = 0.95)
}

# Panel c — perturbation and genetic anchors.
p <- add_panel_heading(
  p, "c", "Orthogonal anchors",
  "Independent evidence constrains mechanism without implying causality",
  67.0, 72.1
)

p <- p +
  round_box(66.5, 97.0, 55.2, 65.6, pal[["mauve_soft"]], pal[["mauve"]]) +
  round_box(66.5, 97.0, 43.0, 53.4, pal[["violet_soft"]], pal[["violet"]])

p <- add_text(p, 70.7, 62.6, "6/6", 3.7, pal[["mauve"]], font_display, "bold")
p <- add_text(p, 76.0, 63.1, "TDP-43 PERTURBATION", 1.75, pal[["mauve"]], font_display, "bold", hjust = 0)
p <- add_text(p, 76.0, 59.8, "Glial disease-trend correlations point opposite\nin both independent iPSC-neuron lines", 1.62, pal[["ink"]], hjust = 0, lineheight = 0.98)
p <- add_text(p, 76.0, 56.7, "Exact binomial P = 0.031", 1.5, pal[["muted"]], hjust = 0)

p <- add_text(p, 70.7, 50.4, "0/9", 3.45, pal[["violet"]], font_display, "bold")
p <- add_text(p, 76.0, 51.1, "ALS GWAS PROGRAMME SETS", 1.7, pal[["violet"]], font_display, "bold", hjust = 0)
p <- add_text(p, 76.0, 47.9, "No prespecified set enriched at P < 0.05", 1.62, pal[["ink"]], hjust = 0)
p <- add_text(p, 76.0, 44.8, "MAPT alone passes 0.05/17 (P = 4.8 × 10⁻⁴)", 1.5, pal[["muted"]], hjust = 0)

# Panel d — interpretation boundaries and next tests.
p <- add_panel_heading(
  p, "d", "Interpretation boundary",
  "The final model is deliberately narrower than a universal signature, causal cascade or biomarker claim",
  4, 35.2
)

d_cards <- data.table(
  xmin = c(4.0, 35.5, 67.0),
  xmax = c(32.5, 64.0, 97.0),
  fill = c(pal[["teal_soft"]], pal[["amber_soft"]], pal[["violet_soft"]]),
  border = c(pal[["teal"]], pal[["amber"]], pal[["violet"]]),
  heading = c("SUPPORTED", "NOT ESTABLISHED", "NEXT VALIDATION"),
  body = c(
    "• Recurrent programme families\n• 13 strict C9 replication axes\n• 40 robust external-test axes\n• Bounded chromatin context",
    "• Universal ALS direction (0 axes)\n• Complete RNA–ATAC chain (0)\n• Enhancer-to-gene causality\n• Motif / TF-activity mechanism",
    "• Genetically matched C9 cohort\n• Functional BAG3–HSP90 perturbation\n• Cell-type enhancer validation\n• Junction-resolved RNA analysis"
  )
)

for (i in seq_len(nrow(d_cards))) {
  row <- d_cards[i]
  p <- p + round_box(row$xmin, row$xmax, 5.0, 28.8, row$fill, row$border)
  p <- add_text(p, row$xmin + 2.0, 25.4, row$heading, 2.05, row$border, font_display, "bold", hjust = 0)
  p <- p + annotate("segment", x = row$xmin + 2.0, xend = row$xmax - 2.0,
                    y = 22.9, yend = 22.9, linewidth = 0.4, colour = alpha(row$border, 0.55))
  p <- add_text(p, row$xmin + 2.0, 20.7, row$body, 1.82, pal[["ink"]],
                hjust = 0, vjust = 1, lineheight = 1.38)
}

width_mm <- 190
height_mm <- 205
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
prefix <- file.path(draft_root, "Figure6_redesign_draft")

svglite::svglite(paste0(prefix, ".svg"), width = width_in, height = height_in, bg = "white")
print(p)
dev.off()

grDevices::cairo_pdf(
  paste0(prefix, ".pdf"), width = width_in, height = height_in,
  family = font_body, bg = "white"
)
print(p)
dev.off()

ragg::agg_tiff(
  paste0(prefix, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, background = "white", compression = "lzw"
)
print(p)
dev.off()

ragg::agg_png(
  paste0(prefix, "_preview.png"), width = width_in, height = height_in,
  units = "in", res = 300, background = "white"
)
print(p)
dev.off()

qa_lines <- c(
  "Figure 6 redesign draft — evidence-constrained synthesis",
  "Core conclusion: recurrent ALS–FTD programmes are robust under donor-aware and external testing, but direction remains context dependent and no universal ALS signature or complete RNA–ATAC causal chain is established.",
  sprintf("Panel a: conceptual relationships are undirected; anchors are 13 strict two-cohort axes, 40/126 robust external-test axes, zero complete three-cohort ALS axes, and the 14-donor pTDP multiome results (409 DARs, 406 oligodendroglial; seven peak-gene associations; zero complete chains)."),
  sprintf("Panel b: discovery %d/%d; strict two-cohort replication %d axes; frozen diagnostic test %d/%d; complete RNA–ATAC chains %d. Counts come from distinct universes and are not a sequential funnel.", n_discovery_pass, n_discovery_total, n_strict_replication, n_ruf_robust, n_ruf_total, n_complete_chain),
  sprintf("Panel c: all %d/%d glial correlation directions are opposite across two independent TDP-43 knockdown lines (two-sided exact binomial P %.5f). No ALS GWAS programme set passes P < 0.05 (%d/%d); MAPT alone passes Bonferroni among 17 candidates (P %.6g; threshold %.6g).", n_glial_opposite, n_glial_tests, glial_binom_p, n_gwas_set_hits, n_gwas_sets, gwas_gene_hits$P, bonf_threshold),
  "Panel d: supported statements, unestablished claims and future validation priorities are separated explicitly. Motif/TF-activity interpretation remains a manuscript-level methodological boundary, not a data panel.",
  "No panel treats nuclei as biological replicates, pools incompatible effect sizes, infers causal arrows, or counts correlated same-nucleus layers as independent validation."
)
writeLines(qa_lines, file.path(draft_root, "Figure6_redesign_draft_QA.txt"))

cat("Rendered Figure 6 draft to:", draft_root, "\n")
cat("Discovery programme axes:", n_discovery_pass, "/", n_discovery_total, "\n")
cat("Strict two-cohort replication axes:", n_strict_replication, "\n")
cat("Frozen diagnostic robust axes:", n_ruf_robust, "/", n_ruf_total, "\n")
cat("pTDP DARs:", n_dars, "(", n_odc_dars, "oligodendroglial )\n")
cat("Peak-gene associations:", n_peak_gene, "\n")
cat("Complete RNA-ATAC chains:", n_complete_chain, "\n")
cat("Glial perturbation directions opposite:", n_glial_opposite, "/", n_glial_tests, "\n")
cat("ALS GWAS programme sets passing P < 0.05:", n_gwas_set_hits, "/", n_gwas_sets, "\n")
cat("Bonferroni candidate-gene hits:", paste(gwas_gene_hits$GENE, collapse = ", "), "\n")

# Alternative radial composition. The rectangular draft above is intentionally
# retained so that both layouts can be compared before approval.
circle_xy <- function(cx, cy, radius, n = 240) {
  angle <- seq(0, 2 * pi, length.out = n)
  data.table(x = cx + radius * cos(angle), y = cy + radius * sin(angle))
}

annulus_xy <- function(cx, cy, r_inner, r_outer, start_deg, end_deg, n = 180) {
  angle <- seq(start_deg, end_deg, length.out = n) * pi / 180
  outer <- data.table(x = cx + r_outer * cos(angle), y = cy + r_outer * sin(angle))
  inner <- data.table(
    x = cx + r_inner * cos(rev(angle)),
    y = cy + r_inner * sin(rev(angle))
  )
  rbind(outer, inner)
}

add_circle <- function(plot, cx, cy, radius, fill, colour, linewidth = 0.45) {
  xy <- circle_xy(cx, cy, radius)
  plot + geom_polygon(
    data = xy, aes(x = x, y = y), inherit.aes = FALSE,
    fill = fill, colour = colour, linewidth = linewidth
  )
}

add_annulus <- function(plot, cx, cy, r_inner, r_outer, start_deg, end_deg,
                        fill, colour, linewidth = 0.45) {
  xy <- annulus_xy(cx, cy, r_inner, r_outer, start_deg, end_deg)
  plot + geom_polygon(
    data = xy, aes(x = x, y = y), inherit.aes = FALSE,
    fill = fill, colour = colour, linewidth = linewidth
  )
}

p_radial <- ggplot() +
  coord_fixed(xlim = c(0, 100), ylim = c(0, 108), ratio = 1,
              expand = FALSE, clip = "off") +
  theme_void(base_family = font_body) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(7, 9, 7, 9, unit = "pt")
  )

# Panel a: one circular evidence map. Spokes are deliberately undirected.
p_radial <- add_panel_heading(
  p_radial, "a", "A relational map of ALS–FTD programmes",
  "Concentric layers separate the recurrent programme model from the contexts and tests that constrain it",
  4, 105
)

map_cx <- 33.5
map_cy <- 72.5
orbit_r <- 22.7

# Quiet field and orbit guides.
p_radial <- add_circle(p_radial, map_cx, map_cy, 25.5, "#FAFBFC", pal[["line"]], 0.35)
p_radial <- add_circle(p_radial, map_cx, map_cy, orbit_r, NA, alpha(pal[["neutral"]], 0.65), 0.35) +
  annotate("segment", x = map_cx - 27.0, xend = map_cx + 27.0,
           y = map_cy, yend = map_cy, colour = alpha(pal[["line"]], 0.55), linewidth = 0.3) +
  annotate("segment", x = map_cx, xend = map_cx,
           y = map_cy - 27.0, yend = map_cy + 27.0,
           colour = alpha(pal[["line"]], 0.55), linewidth = 0.3)

nodes <- data.table(
  angle = c(145, 92, 47, 8, -38, -91, -137),
  label = c("PATHOLOGICAL\nCONTEXT", "DISCOVERY", "C9 TEST", "FROZEN\nTEST",
            "pTDP\nMULTIOME", "PERTURBATION", "ALS GWAS"),
  metric = c("C9 / TDP-43", "35/63", "13", "40/126", "409 · 7 · 0", "6/6", "0/9"),
  colour = c(pal[["violet"]], pal[["blue"]], pal[["teal"]], pal[["teal"]],
             pal[["amber"]], pal[["mauve"]], pal[["violet"]]),
  fill = c(pal[["violet_soft"]], pal[["blue_soft"]], pal[["teal_soft"]], pal[["teal_soft"]],
           pal[["amber_soft"]], pal[["mauve_soft"]], pal[["violet_soft"]]),
  radius = c(5.1, 4.6, 4.5, 4.9, 5.2, 4.7, 4.6)
)
nodes[, `:=`(
  x = map_cx + orbit_r * cos(angle * pi / 180),
  y = map_cy + orbit_r * sin(angle * pi / 180)
)]

for (i in seq_len(nrow(nodes))) {
  row <- nodes[i]
  vx <- row$x - map_cx
  vy <- row$y - map_cy
  distance <- sqrt(vx^2 + vy^2)
  p_radial <- p_radial + annotate(
    "segment",
    x = map_cx + 14.5 * vx / distance,
    y = map_cy + 14.5 * vy / distance,
    xend = row$x - row$radius * vx / distance,
    yend = row$y - row$radius * vy / distance,
    linewidth = 0.48, linetype = 2, colour = alpha(row$colour, 0.65)
  )
}

# Two programme-family arcs surround the central conditioned model.
p_radial <- add_annulus(
  p_radial, map_cx, map_cy, 9.2, 14.7, 18, 162,
  pal[["mauve_soft"]], pal[["mauve"]], 0.55
)
p_radial <- add_annulus(
  p_radial, map_cx, map_cy, 9.2, 14.7, 198, 342,
  pal[["blue_soft"]], pal[["blue"]], 0.55
)
p_radial <- add_circle(p_radial, map_cx, map_cy, 8.35, "#F8FAFC", pal[["ink"]], 0.65)

p_radial <- add_text(
  p_radial, map_cx, map_cy + 2.0,
  "CELL × REGION", 1.72, pal[["muted"]], font_display, "bold"
)
p_radial <- add_text(
  p_radial, map_cx, map_cy - 0.8,
  "CONDITIONED", 2.15, pal[["ink"]], font_display, "bold"
)
p_radial <- add_text(
  p_radial, map_cx, map_cy - 3.4,
  "PROGRAMMES", 2.15, pal[["ink"]], font_display, "bold"
)
p_radial <- add_text(
  p_radial, map_cx, map_cy + 11.6,
  "MITOCHONDRIAL /\nPROTEOSTASIS", 1.42, pal[["mauve"]], font_display, "bold",
  lineheight = 0.95
)
p_radial <- add_text(
  p_radial, map_cx, map_cy - 11.5,
  "SMALL-GTPASE / SYNAPTIC /\nAXONAL", 1.34, pal[["blue"]], font_display, "bold",
  lineheight = 0.95
)

for (i in seq_len(nrow(nodes))) {
  row <- nodes[i]
  p_radial <- add_circle(p_radial, row$x, row$y, row$radius, row$fill, row$colour, 0.55)
  p_radial <- add_text(
    p_radial, row$x, row$y + 1.25, row$label,
    1.22, row$colour, font_display, "bold", lineheight = 0.9
  )
  p_radial <- add_text(
    p_radial, row$x, row$y - 1.35, row$metric,
    ifelse(nchar(row$metric) > 8, 1.45, 2.05), pal[["ink"]], font_display, "bold"
  )
}

p_radial <- add_text(
  p_radial, map_cx, 43.5,
  "inner arcs = recurrent families     ·     outer nodes = independent constraints     ·     spokes = association, not causality",
  1.48, pal[["muted"]], hjust = 0.5
)

# Panel b: an open evidence ledger rather than another bank of cards.
p_radial <- add_panel_heading(
  p_radial, "b", "Quantitative anchors",
  "Exact counts decode the outer orbit",
  65.5, 101.5
)

anchor_rows <- data.table(
  y = c(92.5, 84.0, 75.5, 66.5, 57.5, 48.5),
  metric = c("35/63", "13", "40/126", "409 → 7 → 0", "6/6", "0/9"),
  detail = c(
    "global-FDR discovery axes\n16–17 donors",
    "strict dual-FDR + donor-LOO axes\ntwo C9 cohorts",
    "robust external-test axes\nin 79 donors",
    "DARs → peak–gene links → complete chains\n14 donors; 406/409 DARs in ODC",
    "glial directions opposite\nin two TDP-43-KD lines; P = 0.031",
    "ALS GWAS programme sets at P < 0.05\nMAPT alone passes 0.05/17"
  ),
  colour = c(pal[["blue"]], pal[["teal"]], pal[["teal"]], pal[["amber"]],
             pal[["mauve"]], pal[["violet"]])
)

p_radial <- p_radial + annotate(
  "segment", x = 67.5, xend = 67.5, y = 47.0, yend = 94.0,
  linewidth = 0.65, colour = pal[["line"]]
)
for (i in seq_len(nrow(anchor_rows))) {
  row <- anchor_rows[i]
  p_radial <- add_circle(p_radial, 67.5, row$y, 1.0, "white", row$colour, 0.65)
  p_radial <- add_text(
    p_radial, 70.0, row$y + 0.8, row$metric,
    ifelse(nchar(row$metric) > 7, 2.2, 2.75), row$colour, font_display, "bold", hjust = 0
  )
  p_radial <- add_text(
    p_radial, 79.5, row$y + 1.0, row$detail,
    1.45, pal[["ink"]], hjust = 0, vjust = 0.5, lineheight = 0.95
  )
}

# Panel c: three connected circular boundary nodes continue the visual grammar
# of the evidence orbit instead of appearing as a detached editorial strip.
p_radial <- add_panel_heading(
  p_radial, "c", "Interpretation boundary",
  "The synthesis is deliberately narrower than a universal signature, causal cascade or biomarker claim",
  4, 38.3
)

boundary_cols <- data.table(
  cx = c(19.5, 50.0, 80.5),
  cy = c(17.3, 17.3, 17.3),
  radius = c(13.7, 13.7, 13.7),
  heading = c("SUPPORTED", "NOT ESTABLISHED", "NEXT VALIDATION"),
  colour = c(pal[["teal"]], pal[["amber"]], pal[["violet"]]),
  fill = c(pal[["teal_soft"]], pal[["amber_soft"]], pal[["violet_soft"]]),
  body = c(
    "• Recurrent programme families\n• 13 strict C9 axes\n• 40 robust external-test axes\n• Bounded chromatin context",
    "• Universal ALS direction (0 axes)\n• Complete RNA–ATAC chain (0)\n• Enhancer-to-gene causality\n• Motif / TF-activity mechanism",
    "• Genetically matched C9 cohort\n• Functional BAG3–HSP90 perturbation\n• Cell-type enhancer validation\n• Junction-resolved RNA analysis"
  )
)

# Two quiet curves visually fold the analytical map and evidence ledger into
# the shared interpretation boundary; absence of arrows avoids causal meaning.
p_radial <- p_radial +
  annotate(
    "curve", x = 33.5, y = 43.2, xend = 24.0, yend = 30.4,
    curvature = 0.18, linetype = 2, linewidth = 0.45,
    colour = alpha(pal[["neutral"]], 0.7)
  ) +
  annotate(
    "curve", x = 67.5, y = 46.5, xend = 76.0, yend = 30.4,
    curvature = -0.18, linetype = 2, linewidth = 0.45,
    colour = alpha(pal[["neutral"]], 0.7)
  ) +
  annotate(
    "segment", x = 19.5, xend = 80.5, y = 17.3, yend = 17.3,
    linewidth = 0.55, linetype = 2, colour = alpha(pal[["line"]], 0.9)
  )

for (i in seq_len(nrow(boundary_cols))) {
  row <- boundary_cols[i]
  p_radial <- add_circle(
    p_radial, row$cx, row$cy, row$radius,
    alpha(row$fill, 0.82), row$colour, 0.65
  )
  p_radial <- add_circle(
    p_radial, row$cx, row$cy + 10.4, 1.15,
    "white", row$colour, 0.55
  )
  p_radial <- add_text(
    p_radial, row$cx, row$cy + 7.1, row$heading,
    1.95, row$colour, font_display, "bold", hjust = 0.5
  )
  p_radial <- p_radial + annotate(
    "segment", x = row$cx - 9.0, xend = row$cx + 9.0,
    y = row$cy + 4.7, yend = row$cy + 4.7,
    linewidth = 0.45, colour = alpha(row$colour, 0.55)
  )
  p_radial <- add_text(
    p_radial, row$cx - 9.4, row$cy + 2.7, row$body,
    1.55, pal[["ink"]], hjust = 0, vjust = 1, lineheight = 1.30
  )
}

radial_prefix <- file.path(draft_root, "Figure6_redesign_radial_draft")

svglite::svglite(paste0(radial_prefix, ".svg"), width = width_in, height = height_in, bg = "white")
print(p_radial)
dev.off()

grDevices::cairo_pdf(
  paste0(radial_prefix, ".pdf"), width = width_in, height = height_in,
  family = font_body, bg = "white"
)
print(p_radial)
dev.off()

ragg::agg_tiff(
  paste0(radial_prefix, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, background = "white", compression = "lzw"
)
print(p_radial)
dev.off()

ragg::agg_png(
  paste0(radial_prefix, "_preview.png"), width = width_in, height = height_in,
  units = "in", res = 300, background = "white"
)
print(p_radial)
dev.off()

radial_qa <- c(
  "Figure 6 radial alternative — relation-first synthesis",
  "The centre contains the cell- and region-conditioned model; inner arcs represent the two recurrent programme families; outer orbit nodes represent independent analytical or orthogonal constraints.",
  "All spokes are undirected and explicitly denote association rather than causality.",
  sprintf("Frozen anchors: discovery %d/%d; strict C9 replication %d; external test %d/%d; DAR-to-link-to-chain %d to %d to %d; TDP-43 perturbation %d/%d; ALS GWAS programme sets %d/%d.", n_discovery_pass, n_discovery_total, n_strict_replication, n_ruf_robust, n_ruf_total, n_dars, n_peak_gene, n_complete_chain, n_glial_opposite, n_glial_tests, n_gwas_set_hits, n_gwas_sets),
  "Times New Roman is used throughout. SVG text remains editable."
)
writeLines(radial_qa, file.path(draft_root, "Figure6_redesign_radial_draft_QA.txt"))

cat("Rendered Figure 6 radial alternative to:", draft_root, "\n")

# Fully integrated alternative: the former panel c is no longer a separate
# panel. Its three interpretive states become an outer halo around the model
# and the terminal entries of one continuous evidence-to-interpretation spine.
integrated_width_mm <- 190
integrated_height_mm <- 170
integrated_width_in <- integrated_width_mm / 25.4
integrated_height_in <- integrated_height_mm / 25.4

p_integrated <- ggplot() +
  coord_fixed(xlim = c(0, 100), ylim = c(0, 89.5), ratio = 1,
              expand = FALSE, clip = "off") +
  theme_void(base_family = font_body) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(8, 10, 8, 10, unit = "pt")
  )

p_integrated <- add_text(
  p_integrated, 4.0, 86.0,
  "Evidence-constrained synthesis of ALS–FTD programmes",
  3.7, pal[["ink"]], font_display, "bold", hjust = 0
)
p_integrated <- add_text(
  p_integrated, 4.0, 82.4,
  "One relational system links recurrent programme families to donor-aware tests, orthogonal anchors and explicit interpretive limits",
  1.85, pal[["muted"]], hjust = 0
)
p_integrated <- add_text(
  p_integrated, 5.5, 77.4, "RELATIONAL MODEL",
  1.55, pal[["muted"]], font_display, "bold", hjust = 0
)
p_integrated <- add_text(
  p_integrated, 67.2, 77.4, "EVIDENCE  →  INTERPRETATION",
  1.55, pal[["muted"]], font_display, "bold", hjust = 0
)

integrated_cx <- 34.5
integrated_cy <- 44.5
integrated_orbit <- 19.4

# Interpretation halo: three outer arcs now carry the former panel-c logic.
p_integrated <- add_annulus(
  p_integrated, integrated_cx, integrated_cy, 26.2, 28.1, 28, 142,
  alpha(pal[["teal_soft"]], 0.90), pal[["teal"]], 0.7
)
p_integrated <- add_annulus(
  p_integrated, integrated_cx, integrated_cy, 26.2, 28.1, 148, 262,
  alpha(pal[["amber_soft"]], 0.90), pal[["amber"]], 0.7
)
p_integrated <- add_annulus(
  p_integrated, integrated_cx, integrated_cy, 26.2, 28.1, 268, 382,
  alpha(pal[["violet_soft"]], 0.90), pal[["violet"]], 0.7
)

# A quiet model field and evidence orbit sit inside the interpretation halo.
p_integrated <- add_circle(
  p_integrated, integrated_cx, integrated_cy, 24.6,
  "#FBFCFD", pal[["line"]], 0.35
)
p_integrated <- add_circle(
  p_integrated, integrated_cx, integrated_cy, integrated_orbit,
  NA, alpha(pal[["neutral"]], 0.72), 0.35
)

integrated_nodes <- data.table(
  angle = c(150, 92, 34, -8, -70, -132),
  label = c("CONTEXT", "DISCOVERY", "C9 TEST", "FROZEN TEST", "pTDP MULTIOME", "ORTHOGONAL"),
  metric = c("C9 / TDP-43", "35/63", "13", "40/126", "409 · 7 · 0", "6/6 · 0/9"),
  colour = c(pal[["violet"]], pal[["blue"]], pal[["teal"]], pal[["teal"]],
             pal[["amber"]], pal[["mauve"]]),
  fill = c(pal[["violet_soft"]], pal[["blue_soft"]], pal[["teal_soft"]], pal[["teal_soft"]],
           pal[["amber_soft"]], pal[["mauve_soft"]]),
  radius = c(4.0, 3.9, 3.8, 4.1, 4.2, 4.1)
)
integrated_nodes[, `:=`(
  x = integrated_cx + integrated_orbit * cos(angle * pi / 180),
  y = integrated_cy + integrated_orbit * sin(angle * pi / 180)
)]

for (i in seq_len(nrow(integrated_nodes))) {
  row <- integrated_nodes[i]
  vx <- row$x - integrated_cx
  vy <- row$y - integrated_cy
  distance <- sqrt(vx^2 + vy^2)
  p_integrated <- p_integrated + annotate(
    "segment",
    x = integrated_cx + 14.3 * vx / distance,
    y = integrated_cy + 14.3 * vy / distance,
    xend = row$x - row$radius * vx / distance,
    yend = row$y - row$radius * vy / distance,
    linewidth = 0.42, linetype = 2, colour = alpha(row$colour, 0.62)
  )
}

# Recurrent programme-family ring and central conditioned model.
p_integrated <- add_annulus(
  p_integrated, integrated_cx, integrated_cy, 8.8, 13.9, 18, 162,
  pal[["mauve_soft"]], pal[["mauve"]], 0.55
)
p_integrated <- add_annulus(
  p_integrated, integrated_cx, integrated_cy, 8.8, 13.9, 198, 342,
  pal[["blue_soft"]], pal[["blue"]], 0.55
)
p_integrated <- add_circle(
  p_integrated, integrated_cx, integrated_cy, 8.0,
  "#F8FAFC", pal[["ink"]], 0.7
)

p_integrated <- add_text(
  p_integrated, integrated_cx, integrated_cy + 1.9,
  "CELL × REGION", 1.55, pal[["muted"]], font_display, "bold"
)
p_integrated <- add_text(
  p_integrated, integrated_cx, integrated_cy - 0.9,
  "CONDITIONED", 1.95, pal[["ink"]], font_display, "bold"
)
p_integrated <- add_text(
  p_integrated, integrated_cx, integrated_cy - 3.25,
  "PROGRAMMES", 1.95, pal[["ink"]], font_display, "bold"
)
p_integrated <- add_text(
  p_integrated, integrated_cx, integrated_cy + 11.1,
  "MITOCHONDRIAL /\nPROTEOSTASIS", 1.24, pal[["mauve"]], font_display, "bold",
  lineheight = 0.92
)
p_integrated <- add_text(
  p_integrated, integrated_cx, integrated_cy - 11.0,
  "SMALL-GTPASE / SYNAPTIC /\nAXONAL", 1.18, pal[["blue"]], font_display, "bold",
  lineheight = 0.92
)

for (i in seq_len(nrow(integrated_nodes))) {
  row <- integrated_nodes[i]
  p_integrated <- add_circle(
    p_integrated, row$x, row$y, row$radius,
    alpha(row$fill, 0.96), row$colour, 0.58
  )
  p_integrated <- add_text(
    p_integrated, row$x, row$y + 1.05, row$label,
    ifelse(nchar(row$label) > 11, 0.98, 1.12), row$colour,
    font_display, "bold"
  )
  p_integrated <- add_text(
    p_integrated, row$x, row$y - 1.15, row$metric,
    ifelse(nchar(row$metric) > 8, 1.28, 1.70), pal[["ink"]],
    font_display, "bold"
  )
}

# Direct halo labels make the integrated interpretation layer explicit.
p_integrated <- add_text(
  p_integrated, integrated_cx, 73.6, "SUPPORTED",
  1.50, pal[["teal"]], font_display, "bold"
)
p_integrated <- add_text(
  p_integrated, 7.0, 23.0, "NOT ESTABLISHED",
  1.42, pal[["amber"]], font_display, "bold", hjust = 0
)
p_integrated <- add_text(
  p_integrated, 50.0, 18.9, "NEXT VALIDATION",
  1.42, pal[["violet"]], font_display, "bold", hjust = 0
)
p_integrated <- add_text(
  p_integrated, integrated_cx, 12.5,
  "outer halo = interpretation     ·     nodes = analytical constraints     ·     spokes = association, not causality",
  1.28, pal[["muted"]]
)

# One continuous evidence-to-interpretation spine replaces separate panels b/c.
spine_x <- 68.5
p_integrated <- p_integrated +
  annotate(
    "segment", x = 62.7, xend = spine_x, y = integrated_cy, yend = integrated_cy,
    linewidth = 0.55, linetype = 2, colour = pal[["line"]]
  ) +
  annotate(
    "segment", x = spine_x, xend = spine_x, y = 9.0, yend = 72.5,
    linewidth = 0.7, colour = pal[["line"]]
  )

ledger_evidence <- data.table(
  y = c(69.0, 61.5, 54.0, 46.5, 39.0),
  metric = c("35/63", "13", "40/126", "409 → 7 → 0", "6/6 · 0/9"),
  detail = c(
    "global-FDR discovery axes · 16–17 donors",
    "strict dual-FDR + donor-LOO axes · two C9 cohorts",
    "robust external-test axes · 79 donors",
    "DARs → peak–gene links → complete chains · 14 donors",
    "opposite glial perturbation directions · ALS GWAS sets"
  ),
  foot = c("", "", "0 universal three-cohort ALS axes",
           "406/409 DARs in ODC", "P = 0.031 · MAPT alone passes 0.05/17"),
  colour = c(pal[["blue"]], pal[["teal"]], pal[["teal"]], pal[["amber"]], pal[["violet"]])
)

for (i in seq_len(nrow(ledger_evidence))) {
  row <- ledger_evidence[i]
  p_integrated <- add_circle(p_integrated, spine_x, row$y, 0.92, "white", row$colour, 0.62)
  p_integrated <- add_text(
    p_integrated, 71.0, row$y + 0.65, row$metric,
    ifelse(nchar(row$metric) > 8, 1.90, 2.35), row$colour,
    font_display, "bold", hjust = 0
  )
  p_integrated <- add_text(
    p_integrated, 79.4, row$y + 0.65, row$detail,
    1.23, pal[["ink"]], hjust = 0
  )
  if (nzchar(row$foot)) {
    p_integrated <- add_text(
      p_integrated, 79.4, row$y - 1.45, row$foot,
      1.10, pal[["muted"]], hjust = 0
    )
  }
}

p_integrated <- p_integrated + annotate(
  "segment", x = spine_x, xend = 97.0, y = 34.4, yend = 34.4,
  linewidth = 0.35, colour = pal[["line"]]
)
p_integrated <- add_text(
  p_integrated, 71.0, 32.4, "INTERPRETATION HALO",
  1.30, pal[["muted"]], font_display, "bold", hjust = 0
)

ledger_status <- data.table(
  y = c(27.5, 19.5, 11.5),
  heading = c("SUPPORTED", "NOT ESTABLISHED", "NEXT VALIDATION"),
  detail = c(
    "recurrent families · 13 strict C9 axes\n40/126 robust external-test axes",
    "0 universal ALS axes · 0 complete chains\nno enhancer causality or TF mechanism",
    "matched C9 cohort · BAG3–HSP90 function\nenhancer validation · junction-resolved RNA"
  ),
  colour = c(pal[["teal"]], pal[["amber"]], pal[["violet"]]),
  fill = c(pal[["teal_soft"]], pal[["amber_soft"]], pal[["violet_soft"]])
)

for (i in seq_len(nrow(ledger_status))) {
  row <- ledger_status[i]
  p_integrated <- add_circle(p_integrated, spine_x, row$y, 1.15, row$fill, row$colour, 0.68)
  p_integrated <- add_text(
    p_integrated, 71.0, row$y + 0.9, row$heading,
    1.42, row$colour, font_display, "bold", hjust = 0
  )
  p_integrated <- add_text(
    p_integrated, 71.0, row$y - 1.1, row$detail,
    1.14, pal[["ink"]], hjust = 0, lineheight = 0.95
  )
}

integrated_prefix <- file.path(draft_root, "Figure6_redesign_integrated_draft")

svglite::svglite(
  paste0(integrated_prefix, ".svg"),
  width = integrated_width_in, height = integrated_height_in, bg = "white"
)
print(p_integrated)
dev.off()

grDevices::cairo_pdf(
  paste0(integrated_prefix, ".pdf"),
  width = integrated_width_in, height = integrated_height_in,
  family = font_body, bg = "white"
)
print(p_integrated)
dev.off()

ragg::agg_tiff(
  paste0(integrated_prefix, ".tiff"),
  width = integrated_width_in, height = integrated_height_in,
  units = "in", res = 600, background = "white", compression = "lzw"
)
print(p_integrated)
dev.off()

ragg::agg_png(
  paste0(integrated_prefix, "_preview.png"),
  width = integrated_width_in, height = integrated_height_in,
  units = "in", res = 300, background = "white"
)
print(p_integrated)
dev.off()

integrated_qa <- c(
  "Figure 6 fully integrated draft — evidence and interpretation in one relational system",
  "The former panel c has been removed as an independent panel.",
  "The outer halo encodes supported evidence, unestablished claims and next validation, while the right-hand spine continues from quantitative anchors into the same three interpretive states.",
  "No arrows are used; dashed spokes indicate association rather than causality.",
  sprintf("Frozen anchors: discovery %d/%d; strict C9 replication %d; external test %d/%d; DAR-to-link-to-chain %d to %d to %d; TDP-43 perturbation %d/%d; ALS GWAS programme sets %d/%d.", n_discovery_pass, n_discovery_total, n_strict_replication, n_ruf_robust, n_ruf_total, n_dars, n_peak_gene, n_complete_chain, n_glial_opposite, n_glial_tests, n_gwas_set_hits, n_gwas_sets),
  "Times New Roman is used throughout. SVG text remains editable."
)
writeLines(
  integrated_qa,
  file.path(draft_root, "Figure6_redesign_integrated_draft_QA.txt")
)

cat("Rendered fully integrated Figure 6 alternative to:", draft_root, "\n")

# Fused circular alternative: no in-figure title and no detached right-hand
# ledger. Every frozen quantitative anchor is placed inside its evidence node;
# interpretation remains the outermost halo of the same relational system.
fused_width_mm <- 190
fused_height_mm <- 182
fused_width_in <- fused_width_mm / 25.4
fused_height_in <- fused_height_mm / 25.4

p_fused <- ggplot() +
  coord_fixed(xlim = c(0, 100), ylim = c(0, 84.2), ratio = 1,
              expand = FALSE, clip = "off") +
  theme_void(base_family = font_body) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(7, 9, 7, 9, unit = "pt")
  )

fused_cx <- 50.0
fused_cy <- 45.0
fused_orbit <- 22.0

# The former interpretation panel becomes a three-part halo.
p_fused <- add_annulus(
  p_fused, fused_cx, fused_cy, 31.0, 33.0, 30, 142,
  alpha(pal[["teal_soft"]], 0.92), pal[["teal"]], 0.72
)
p_fused <- add_annulus(
  p_fused, fused_cx, fused_cy, 31.0, 33.0, 148, 262,
  alpha(pal[["amber_soft"]], 0.92), pal[["amber"]], 0.72
)
p_fused <- add_annulus(
  p_fused, fused_cx, fused_cy, 31.0, 33.0, 268, 382,
  alpha(pal[["violet_soft"]], 0.92), pal[["violet"]], 0.72
)

# Interpretation labels are attached directly to their halo segments.
p_fused <- add_text(
  p_fused, fused_cx, 81.0, "SUPPORTED",
  1.62, pal[["teal"]], font_display, "bold"
)
p_fused <- add_text(
  p_fused, fused_cx, 78.2,
  "13 strict C9 axes · 40/126 robust external-test axes",
  1.18, pal[["teal"]]
)
p_fused <- add_text(
  p_fused, 17.0, 9.8, "NOT ESTABLISHED",
  1.52, pal[["amber"]], font_display, "bold", hjust = 0
)
p_fused <- add_text(
  p_fused, 17.0, 7.1,
  "0 universal ALS axes · 0 complete RNA–ATAC chains",
  1.12, pal[["amber"]], hjust = 0
)
p_fused <- add_text(
  p_fused, 59.2, 9.8, "NEXT VALIDATION",
  1.52, pal[["violet"]], font_display, "bold", hjust = 0
)
p_fused <- add_text(
  p_fused, 59.2, 7.1,
  "matched C9 · functional perturbation · enhancer + junction RNA",
  1.08, pal[["violet"]], hjust = 0
)

# Quiet model field and evidence orbit.
p_fused <- add_circle(
  p_fused, fused_cx, fused_cy, 29.3,
  "#FBFCFD", pal[["line"]], 0.35
)
p_fused <- add_circle(
  p_fused, fused_cx, fused_cy, fused_orbit,
  NA, alpha(pal[["neutral"]], 0.72), 0.38
)

fused_nodes <- data.table(
  angle = c(150, 105, 55, -8, -70, -100, -145),
  label = c(
    "PATHOLOGICAL\nCONTEXT",
    "DISCOVERY\nGLOBAL-FDR",
    "STRICT C9\nDUAL-FDR + LOO",
    "FROZEN TEST\n79 DONORS",
    "pTDP MULTIOME\n14 DONORS",
    "TDP-43 KD\nOPPOSITE",
    "ALS GWAS\nMAPT ONLY"
  ),
  metric = c("C9 / TDP-43", "35/63", "13", "40/126", "409 → 7 → 0", "6/6", "0/9"),
  colour = c(pal[["violet"]], pal[["blue"]], pal[["teal"]], pal[["teal"]],
             pal[["amber"]], pal[["mauve"]], pal[["violet"]]),
  fill = c(pal[["violet_soft"]], pal[["blue_soft"]], pal[["teal_soft"]], pal[["teal_soft"]],
           pal[["amber_soft"]], pal[["mauve_soft"]], pal[["violet_soft"]]),
  radius = c(4.7, 4.5, 4.5, 4.7, 5.0, 4.5, 4.5),
  label_size = c(1.08, 1.00, 0.94, 1.00, 0.96, 1.00, 1.00),
  metric_size = c(1.38, 1.82, 1.85, 1.78, 1.34, 1.78, 1.78)
)
fused_nodes[, `:=`(
  x = fused_cx + fused_orbit * cos(angle * pi / 180),
  y = fused_cy + fused_orbit * sin(angle * pi / 180)
)]

for (i in seq_len(nrow(fused_nodes))) {
  row <- fused_nodes[i]
  vx <- row$x - fused_cx
  vy <- row$y - fused_cy
  distance <- sqrt(vx^2 + vy^2)
  p_fused <- p_fused + annotate(
    "segment",
    x = fused_cx + 15.3 * vx / distance,
    y = fused_cy + 15.3 * vy / distance,
    xend = row$x - row$radius * vx / distance,
    yend = row$y - row$radius * vy / distance,
    linewidth = 0.48, linetype = 2, colour = alpha(row$colour, 0.62)
  )
}

# Recurrent programme families and the conditioned core.
p_fused <- add_annulus(
  p_fused, fused_cx, fused_cy, 9.4, 14.9, 18, 162,
  pal[["mauve_soft"]], pal[["mauve"]], 0.58
)
p_fused <- add_annulus(
  p_fused, fused_cx, fused_cy, 9.4, 14.9, 198, 342,
  pal[["blue_soft"]], pal[["blue"]], 0.58
)
p_fused <- add_circle(
  p_fused, fused_cx, fused_cy, 8.5,
  "#F8FAFC", pal[["ink"]], 0.72
)

p_fused <- add_text(
  p_fused, fused_cx, fused_cy + 2.1,
  "CELL × REGION", 1.62, pal[["muted"]], font_display, "bold"
)
p_fused <- add_text(
  p_fused, fused_cx, fused_cy - 0.9,
  "CONDITIONED", 2.10, pal[["ink"]], font_display, "bold"
)
p_fused <- add_text(
  p_fused, fused_cx, fused_cy - 3.5,
  "PROGRAMMES", 2.10, pal[["ink"]], font_display, "bold"
)
p_fused <- add_text(
  p_fused, fused_cx, fused_cy + 11.9,
  "MITOCHONDRIAL /\nPROTEOSTASIS", 1.35, pal[["mauve"]], font_display, "bold",
  lineheight = 0.92
)
p_fused <- add_text(
  p_fused, fused_cx, fused_cy - 11.8,
  "SMALL-GTPASE / SYNAPTIC /\nAXONAL", 1.27, pal[["blue"]], font_display, "bold",
  lineheight = 0.92
)

for (i in seq_len(nrow(fused_nodes))) {
  row <- fused_nodes[i]
  p_fused <- add_circle(
    p_fused, row$x, row$y, row$radius,
    alpha(row$fill, 0.97), row$colour, 0.60
  )
  p_fused <- add_text(
    p_fused, row$x, row$y + 1.45, row$label,
    row$label_size, row$colour, font_display, "bold", lineheight = 0.90
  )
  p_fused <- add_text(
    p_fused, row$x, row$y - 1.35, row$metric,
    row$metric_size, pal[["ink"]], font_display, "bold"
  )
}

p_fused <- add_text(
  p_fused, fused_cx, 2.7,
  "inner arcs = recurrent programme families     ·     dashed spokes = association, not causality",
  1.22, pal[["muted"]]
)

fused_prefix <- file.path(draft_root, "Figure6_redesign_fused_draft")

svglite::svglite(
  paste0(fused_prefix, ".svg"),
  width = fused_width_in, height = fused_height_in, bg = "white"
)
print(p_fused)
dev.off()

grDevices::cairo_pdf(
  paste0(fused_prefix, ".pdf"),
  width = fused_width_in, height = fused_height_in,
  family = font_body, bg = "white"
)
print(p_fused)
dev.off()

ragg::agg_tiff(
  paste0(fused_prefix, ".tiff"),
  width = fused_width_in, height = fused_height_in,
  units = "in", res = 600, background = "white", compression = "lzw"
)
print(p_fused)
dev.off()

ragg::agg_png(
  paste0(fused_prefix, "_preview.png"),
  width = fused_width_in, height = fused_height_in,
  units = "in", res = 300, background = "white"
)
print(p_fused)
dev.off()

fused_qa <- c(
  "Figure 6 fused circular draft — no in-figure title and no detached ledger",
  "All frozen numbers are embedded in their corresponding evidence nodes.",
  "The outer halo integrates supported evidence, unestablished claims and next validation into the same relational system.",
  "No arrows are used; dashed spokes indicate association rather than causality.",
  sprintf("Frozen anchors: discovery %d/%d; strict C9 replication %d; external test %d/%d; DAR-to-link-to-chain %d to %d to %d; TDP-43 perturbation %d/%d; ALS GWAS programme sets %d/%d.", n_discovery_pass, n_discovery_total, n_strict_replication, n_ruf_robust, n_ruf_total, n_dars, n_peak_gene, n_complete_chain, n_glial_opposite, n_glial_tests, n_gwas_set_hits, n_gwas_sets),
  "Times New Roman is used throughout. SVG text remains editable."
)
writeLines(fused_qa, file.path(draft_root, "Figure6_redesign_fused_draft_QA.txt"))

cat("Rendered title-free fused Figure 6 alternative to:", draft_root, "\n")

# Enriched fused revision: exact seven-fold angular spacing, complete result
# context inside each node, and editable character-by-character curved labels.
curved_text_df <- function(label, cx, cy, radius, start_deg, end_deg,
                           rotation_offset) {
  tokens <- strsplit(label, "[[:space:]]+", perl = TRUE)[[1]]
  token_width <- pmax(nchar(tokens, type = "width"), 1) + 0.9
  token_position <- (cumsum(token_width) - token_width / 2) / sum(token_width)
  angles <- start_deg + (end_deg - start_deg) * token_position
  radians <- angles * pi / 180
  data.table(
    x = cx + radius * cos(radians),
    y = cy + radius * sin(radians),
    token = tokens,
    rotation = angles + rotation_offset
  )
}

add_curved_text <- function(plot, label, cx, cy, radius, start_deg, end_deg,
                            rotation_offset, colour, size = 1.0,
                            family = font_display, face = "bold") {
  letters <- curved_text_df(
    label, cx, cy, radius, start_deg, end_deg, rotation_offset
  )
  plot + geom_text(
    data = letters,
    aes(x = x, y = y, label = token, angle = rotation),
    inherit.aes = FALSE, family = family, fontface = face,
    colour = colour, size = size, hjust = 0.5, vjust = 0.5
  )
}

p_enriched <- ggplot() +
  coord_fixed(xlim = c(0, 100), ylim = c(0, 95.8), ratio = 1,
              expand = FALSE, clip = "off") +
  theme_void(base_family = font_body) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(7, 9, 7, 9, unit = "pt")
  )

enriched_cx <- 50.0
enriched_cy <- 48.0
enriched_orbit <- 25.4

# Three interpretation arcs; their explanations are rendered along the arc.
p_enriched <- add_annulus(
  p_enriched, enriched_cx, enriched_cy, 34.5, 36.8, 30, 142,
  alpha(pal[["teal_soft"]], 0.92), pal[["teal"]], 0.72
)
p_enriched <- add_annulus(
  p_enriched, enriched_cx, enriched_cy, 34.5, 36.8, 148, 262,
  alpha(pal[["amber_soft"]], 0.92), pal[["amber"]], 0.72
)
p_enriched <- add_annulus(
  p_enriched, enriched_cx, enriched_cy, 34.5, 36.8, 268, 382,
  alpha(pal[["violet_soft"]], 0.92), pal[["violet"]], 0.72
)

p_enriched <- add_curved_text(
  p_enriched,
  "SUPPORTED",
  enriched_cx, enriched_cy, 38.3, 101, 79, -90,
  pal[["teal"]], size = 1.62
)
p_enriched <- add_curved_text(
  p_enriched,
  "13 STRICT C9 AXES · 40/126 ROBUST EXTERNAL-TEST AXES",
  enriched_cx, enriched_cy, 35.6, 122, 58, -90,
  pal[["teal"]], size = 1.16, face = "plain"
)
p_enriched <- add_curved_text(
  p_enriched,
  "NOT ESTABLISHED",
  enriched_cx, enriched_cy, 38.3, 216, 240, -270,
  pal[["amber"]], size = 1.48
)
p_enriched <- add_curved_text(
  p_enriched,
  "0 UNIVERSAL ALS AXES · 0 COMPLETE RNA-ATAC CHAINS",
  enriched_cx, enriched_cy, 35.6, 205, 245, -270,
  pal[["amber"]], size = 1.10, face = "plain"
)
p_enriched <- add_curved_text(
  p_enriched,
  "NEXT VALIDATION",
  enriched_cx, enriched_cy, 38.3, 324, 300, -270,
  pal[["violet"]], size = 1.48
)
p_enriched <- add_curved_text(
  p_enriched,
  "MATCHED C9 · FUNCTIONAL PERTURBATION · ENHANCER + JUNCTION RNA",
  enriched_cx, enriched_cy, 35.6, 340, 296, -270,
  pal[["violet"]], size = 1.04, face = "plain"
)

# Model field and a seven-position evidence orbit.
p_enriched <- add_circle(
  p_enriched, enriched_cx, enriched_cy, 33.1,
  "#FBFCFD", pal[["line"]], 0.35
)
p_enriched <- add_circle(
  p_enriched, enriched_cx, enriched_cy, enriched_orbit,
  NA, alpha(pal[["neutral"]], 0.72), 0.38
)

equal_angles <- 90 - (0:6) * (360 / 7)
enriched_nodes <- data.table(
  angle = equal_angles,
  label = c(
    "DISCOVERY",
    "C9 REPLICATION",
    "FROZEN TEST",
    "pTDP MULTIOME",
    "TDP-43 PERTURBATION",
    "ALS GWAS",
    "PATHOLOGICAL CONTEXT"
  ),
  metric = c("35/63", "13", "40/126", "409 → 7 → 0", "6/6", "0/9", "C9 / TDP-43"),
  detail = c(
    "global-FDR programme axes\n16–17 donors",
    "dual-FDR + donor-LOO\ntwo independent C9 cohorts",
    "robust axes · 79 donors\n0 universal ALS axes",
    "409 DARs (406 ODC)\n7 links · 0 complete chains",
    "opposite glial directions\ntwo lines · exact P = 0.031",
    "0/9 sets at P < 0.05\nMAPT P = 4.8 × 10⁻⁴",
    "C9orf72 and TDP-43\ncell/region conditioned"
  ),
  interpretation = c(
    "Recurrent systems emerge\nunder donor-aware testing",
    "Replication is strong but\nrestricted to C9 context",
    "Robust across testing,\nnot universally directed",
    "Chromatin context is present;\ncomplete chain is absent",
    "Direction is model-dependent,\nnot a simple causal match",
    "Limited common-variant support;\nMAPT is the exception",
    "Cell, region and cohort\ncondition programme direction"
  ),
  colour = c(
    pal[["blue"]], pal[["teal"]], pal[["teal"]], pal[["amber"]],
    pal[["mauve"]], pal[["violet"]], pal[["violet"]]
  ),
  fill = c(
    pal[["blue_soft"]], pal[["teal_soft"]], pal[["teal_soft"]], pal[["amber_soft"]],
    pal[["mauve_soft"]], pal[["violet_soft"]], pal[["violet_soft"]]
  ),
  metric_size = c(2.38, 2.48, 2.32, 1.82, 2.36, 2.36, 1.76)
)
enriched_nodes[, `:=`(
  x = enriched_cx + enriched_orbit * cos(angle * pi / 180),
  y = enriched_cy + enriched_orbit * sin(angle * pi / 180)
)]

node_radius <- 7.9
for (i in seq_len(nrow(enriched_nodes))) {
  row <- enriched_nodes[i]
  vx <- row$x - enriched_cx
  vy <- row$y - enriched_cy
  distance <- sqrt(vx^2 + vy^2)
  p_enriched <- p_enriched + annotate(
    "segment",
    x = enriched_cx + 17.15 * vx / distance,
    y = enriched_cy + 17.15 * vy / distance,
    xend = row$x - node_radius * vx / distance,
    yend = row$y - node_radius * vy / distance,
    linewidth = 0.48, linetype = 2, colour = alpha(row$colour, 0.62)
  )
}

# Recurrent systems and their context-dependent central interpretation.
p_enriched <- add_annulus(
  p_enriched, enriched_cx, enriched_cy, 10.5, 16.8, 18, 162,
  pal[["mauve_soft"]], pal[["mauve"]], 0.58
)
p_enriched <- add_annulus(
  p_enriched, enriched_cx, enriched_cy, 10.5, 16.8, 198, 342,
  pal[["blue_soft"]], pal[["blue"]], 0.58
)
p_enriched <- add_circle(
  p_enriched, enriched_cx, enriched_cy, 9.6,
  "#F8FAFC", pal[["ink"]], 0.72
)

p_enriched <- add_text(
  p_enriched, enriched_cx, enriched_cy + 5.65,
  "CORE CONCLUSION", 1.18, pal[["teal"]], font_display, "bold"
)
p_enriched <- add_text(
  p_enriched, enriched_cx, enriched_cy + 3.15,
  "ALS–FTD PROGRAMMES", 1.90, pal[["ink"]], font_display, "bold"
)
p_enriched <- add_text(
  p_enriched, enriched_cx, enriched_cy + 1.05,
  "TRANSFER BY FAMILY", 2.08, pal[["ink"]], font_display, "bold"
)
p_enriched <- add_text(
  p_enriched, enriched_cx, enriched_cy - 1.15,
  "but depend on cellular", 1.22, pal[["muted"]]
)
p_enriched <- add_text(
  p_enriched, enriched_cx, enriched_cy - 2.65,
  "and cohort context", 1.22, pal[["muted"]]
)
p_enriched <- add_text(
  p_enriched, enriched_cx, enriched_cy - 4.45,
  "External, perturbational and genetic audit", 1.04, pal[["blue"]]
)
p_enriched <- add_text(
  p_enriched, enriched_cx, enriched_cy - 5.95,
  "must precede biomarker or mechanism claims", 1.04, pal[["blue"]]
)

p_enriched <- add_text(
  p_enriched, enriched_cx, enriched_cy + 14.0,
  "MITOCHONDRIAL / PROTEOSTASIS", 1.48,
  pal[["mauve"]], font_display, "bold"
)
p_enriched <- add_text(
  p_enriched, enriched_cx, enriched_cy + 12.0,
  "recurrent across C9 cohorts", 1.18, pal[["ink"]]
)
p_enriched <- add_text(
  p_enriched, enriched_cx, enriched_cy - 11.6,
  "SMALL-GTPASE / SYNAPTIC / AXONAL", 1.40,
  pal[["blue"]], font_display, "bold"
)
p_enriched <- add_text(
  p_enriched, enriched_cx, enriched_cy - 13.7,
  "contextual in external testing", 1.16, pal[["ink"]]
)

for (i in seq_len(nrow(enriched_nodes))) {
  row <- enriched_nodes[i]
  p_enriched <- add_circle(
    p_enriched, row$x, row$y, node_radius,
    alpha(row$fill, 0.97), row$colour, 0.62
  )
  p_enriched <- add_text(
    p_enriched, row$x, row$y + 5.55, row$label,
    ifelse(nchar(row$label) > 18, 1.28, 1.43),
    row$colour, font_display, "bold"
  )
  p_enriched <- add_text(
    p_enriched, row$x, row$y + 3.05, row$metric,
    row$metric_size, pal[["ink"]], font_display, "bold"
  )
  p_enriched <- add_text(
    p_enriched, row$x, row$y + 0.35, row$detail,
    1.22, pal[["ink"]], hjust = 0.5, lineheight = 0.93
  )
  p_enriched <- p_enriched + annotate(
    "segment",
    x = row$x - 4.75, xend = row$x + 4.75,
    y = row$y - 1.82, yend = row$y - 1.82,
    linewidth = 0.30, colour = alpha(row$colour, 0.45)
  )
  p_enriched <- add_text(
    p_enriched, row$x, row$y - 2.72, "INTERPRETATION",
    1.00, row$colour, font_display, "bold"
  )
  p_enriched <- add_text(
    p_enriched, row$x, row$y - 4.85, row$interpretation,
    1.14, pal[["ink"]], hjust = 0.5, lineheight = 0.91
  )
}

p_enriched <- add_text(
  p_enriched, enriched_cx, 3.3,
  "Dashed spokes denote evidence relationships only; no causal direction is implied",
  1.38, pal[["muted"]]
)

# Overwrite the fused draft outputs with this enriched, evenly spaced revision.
svglite::svglite(
  paste0(fused_prefix, ".svg"),
  width = fused_width_in, height = fused_height_in, bg = "white"
)
print(p_enriched)
dev.off()

grDevices::cairo_pdf(
  paste0(fused_prefix, ".pdf"),
  width = fused_width_in, height = fused_height_in,
  family = font_body, bg = "white"
)
print(p_enriched)
dev.off()

ragg::agg_tiff(
  paste0(fused_prefix, ".tiff"),
  width = fused_width_in, height = fused_height_in,
  units = "in", res = 600, background = "white", compression = "lzw"
)
print(p_enriched)
dev.off()

ragg::agg_png(
  paste0(fused_prefix, "_preview.png"),
  width = fused_width_in, height = fused_height_in,
  units = "in", res = 300, background = "white"
)
print(p_enriched)
dev.off()

fused_qa <- c(
  "Figure 6 enriched fused draft — title free, self-contained and evenly spaced",
  "The seven evidence nodes use exact 360/7 angular spacing and a common radius.",
  "Each node carries its frozen metric plus cohort/test context, including sample counts, zero results and exact inferential anchors.",
  "Each node also includes a bounded one-sentence interpretation that does not exceed the observed evidence.",
  "Supported, not established and next validation explanations follow the outer arcs as editable Times New Roman text.",
  "The centre states the manuscript-level conclusion: ALS-FTD programmes transfer by family but remain context dependent, and external, perturbational and genetic audit must precede biomarker or mechanism claims.",
  "No arrows are used; dashed spokes indicate association rather than causality.",
  sprintf("Frozen anchors: discovery %d/%d; strict C9 replication %d; external test %d/%d; DAR-to-link-to-chain %d to %d to %d; TDP-43 perturbation %d/%d; ALS GWAS programme sets %d/%d.", n_discovery_pass, n_discovery_total, n_strict_replication, n_ruf_robust, n_ruf_total, n_dars, n_peak_gene, n_complete_chain, n_glial_opposite, n_glial_tests, n_gwas_set_hits, n_gwas_sets),
  "SVG text remains editable."
)
writeLines(fused_qa, file.path(draft_root, "Figure6_redesign_fused_draft_QA.txt"))

cat("Rendered enriched evenly-spaced Figure 6 revision to:", draft_root, "\n")
