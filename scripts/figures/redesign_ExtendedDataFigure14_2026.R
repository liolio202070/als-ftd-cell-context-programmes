#!/usr/bin/env Rscript

# Submission-ready redesign of Extended Data Figure 14.
# Figure contract:
# Core conclusion: the absent three-cohort ALS programme axis reflects both gate
# attrition and genuine direction reversal, while perturbation and GWAS analyses
# provide bounded orthogonal context rather than universal validation.
# Archetype: quantitative grid with six non-redundant panels.
# Backend: R only. Output: editable SVG/PDF, 600-dpi TIFF and PNG preview.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
project_lib <- Sys.getenv("ALS_FTD_R_LIBRARY", unset = "")
if (nzchar(project_lib) && dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
})

extended_source <- file.path(project_root, "data", "source_data", "extended_figures")
extended_inputs <- file.path(project_root, "data", "figure_inputs", "ExtendedDataFigure14")
output_root <- file.path(project_root, "outputs", "figures", "ExtendedDataFigure14")
source_output <- file.path(output_root, "source_data")
dir.create(source_output, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  waterfall = file.path(extended_source, "ExtendedDataFigure14_panel_a_source_data.tsv"),
  failures = file.path(extended_inputs, "ruf_gate_first_failure_attribution.tsv"),
  attainability = file.path(extended_source, "ExtendedDataFigure14_panel_b_source_data.tsv"),
  perturbation = file.path(extended_source, "ExtendedDataFigure14_panel_c_source_data.tsv"),
  chaperone = file.path(extended_source, "ExtendedDataFigure14_panel_d_source_data.tsv"),
  gwas_sets = file.path(extended_source, "ExtendedDataFigure14_panel_e_source_data.tsv"),
  gwas_genes = file.path(extended_source, "ExtendedDataFigure14_panel_f_source_data.tsv")
)
if (!all(file.exists(unlist(paths)))) stop("One or more frozen Extended Data Figure 14 inputs are missing.")

pal <- c(
  ink = "#172033",
  muted = "#667085",
  grid = "#D7DEE8",
  soft = "#F5F7FA",
  blue = "#2D6FB3",
  blue_soft = "#DDE9F5",
  teal = "#2B948B",
  teal_soft = "#DCEFEB",
  mauve = "#B15677",
  mauve_soft = "#F2DFE7",
  amber = "#D48725",
  amber_soft = "#F7E8D2",
  red = "#D6574B",
  neutral = "#9AA7B7",
  neutral_soft = "#E9EEF4"
)
font_body <- "Arial"
font_display <- "Avenir Next"

theme_ext <- function(base_size = 6.5) {
  theme_classic(base_size = base_size, base_family = font_body) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = pal[["ink"]]),
      axis.ticks = element_line(linewidth = 0.35, colour = pal[["ink"]]),
      axis.title = element_text(size = base_size, colour = pal[["ink"]]),
      axis.text = element_text(size = base_size - 0.35, colour = pal[["ink"]]),
      plot.title = element_text(
        family = font_display, face = "bold", size = 8.2,
        colour = pal[["ink"]], margin = margin(b = 2)
      ),
      plot.subtitle = element_text(
        size = 5.9, colour = pal[["muted"]], lineheight = 1.05,
        margin = margin(b = 4)
      ),
      plot.caption = element_text(size = 5.2, colour = pal[["muted"]], hjust = 0),
      legend.title = element_text(size = 5.6, face = "bold"),
      legend.text = element_text(size = 5.25),
      legend.key.height = grid::unit(2.2, "mm"),
      legend.key.width = grid::unit(3.0, "mm"),
      panel.grid = element_blank(),
      plot.margin = margin(5, 8, 5, 8, "pt"),
      plot.tag = element_text(
        family = font_display, face = "bold", size = 9.2,
        colour = pal[["ink"]]
      )
    )
}
theme_set(theme_ext())

write_clean_tsv <- function(x, path) {
  y <- copy(as.data.table(x))
  text_columns <- names(y)[vapply(y, function(z) is.character(z) || is.factor(z), logical(1))]
  for (column in text_columns) {
    value <- gsub("[\r\n\t]", " ", as.character(y[[column]]), perl = TRUE)
    set(y, j = column, value = value)
  }
  fwrite(y, path, sep = "\t", quote = FALSE, na = "")
}

# a. Five-gate attrition -----------------------------------------------------
waterfall <- fread(paths$waterfall)
waterfall <- waterfall[gate != "strict_base字段对照"]
waterfall[, gate_label := c(
  "All frozen axes", "G1  Family FDR", "G2  Donor LOO",
  "G3  Technical + age", "G4  SVA12 direction", "G5  SVA12 FDR"
)]
waterfall[, gate_label := factor(gate_label, levels = rev(gate_label))]
waterfall[, final_gate := seq_len(.N) == .N]

failures <- fread(paths$failures)[subset == "全部126轴"]
failure_lookup <- setNames(failures$n_axes, failures$first_failed_gate)
failure_caption <- sprintf(
  "First failures: G1 %d; G3 %d; G4 %d; G5 %d.",
  failure_lookup[["G1_全局家族FDR"]], failure_lookup[["G3_技术+年龄方向"]],
  failure_lookup[["G4_SVA12方向"]], failure_lookup[["G5_SVA12全局FDR"]]
)

p_a <- ggplot(waterfall, aes(x = pass_n, y = gate_label)) +
  geom_segment(
    aes(x = 0, xend = pass_n, yend = gate_label),
    linewidth = 3.8, colour = pal[["neutral_soft"]], lineend = "round"
  ) +
  geom_point(aes(colour = final_gate), size = 2.9) +
  geom_text(aes(label = pass_n), hjust = -0.55, size = 2.25, fontface = "bold") +
  scale_colour_manual(values = c(`FALSE` = pal[["blue"]], `TRUE` = pal[["teal"]]), guide = "none") +
  scale_x_continuous(limits = c(0, 150), breaks = c(0, 50, 100, 150), expand = c(0, 0)) +
  labs(
    title = "Five-gate attrition",
    subtitle = "126 axes; donor leave-one-out removed none",
    caption = failure_caption,
    x = "Axes retained (cumulative)", y = NULL
  ) +
  theme(axis.text.y = element_text(size = 5.65), plot.caption = element_text(margin = margin(t = 3)))

# b. Attainability and direction reversal ----------------------------------
attainability <- fread(paths$attainability)
attainability[, direction_class := fcase(
  sign(z_ruf) != sign(z_disc) & sign(z_ruf) != sign(z_git), "Flipped in both C9 cohorts",
  sign(z_ruf) != sign(z_disc) | sign(z_ruf) != sign(z_git), "Flipped in one C9 cohort",
  default = "Concordant"
)]
attainability[, point_label := fifelse(
  powerA > 0.5,
  paste(
    fcase(
      cell_type_label == "Astrocyte", "Astro",
      cell_type_label == "Oligodendrocyte", "Oligo",
      default = cell_type_label
    ),
    gsub("_", " ", family_id), sep = " · "
  ),
  NA_character_
)]
direction_values <- c(
  "Concordant" = pal[["blue"]],
  "One C9 flip" = pal[["amber"]],
  "Both C9 flip" = pal[["red"]]
)
attainability[, direction_class := fcase(
  direction_class == "Flipped in both C9 cohorts", "Both C9 flip",
  direction_class == "Flipped in one C9 cohort", "One C9 flip",
  default = "Concordant"
)]

p_b <- ggplot(attainability, aes(abs(z_ruf), powerA, colour = direction_class)) +
  geom_hline(yintercept = 0.5, linetype = "22", linewidth = 0.35, colour = pal[["muted"]]) +
  geom_point(size = 2.0, alpha = 0.95) +
  ggrepel::geom_text_repel(
    data = attainability[!is.na(point_label)], aes(label = point_label),
    size = 1.65, box.padding = 0.25, point.padding = 0.12,
    min.segment.length = 0, segment.colour = pal[["grid"]],
    show.legend = FALSE, max.overlaps = Inf
  ) +
  scale_colour_manual(values = direction_values, name = NULL) +
  scale_y_continuous(limits = c(0, 1.04), breaks = c(0, 0.5, 1), labels = percent_format(accuracy = 1)) +
  scale_x_continuous(limits = c(0, 6.45), expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title = "Three-cohort gate attainability",
    subtitle = "Median 1.2%; 3/21 axes exceed 50%, all direction-flipped",
    x = "|z| in Ruf (ALS versus control)", y = "Attainable probability"
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    legend.position = "bottom", legend.box.margin = margin(t = -2),
    legend.spacing.x = grid::unit(1.2, "mm")
  )

# c. Perturbation-to-pathology alignment -----------------------------------
perturbation <- fread(paths$perturbation)
perturbation[, dataset_label := factor(
  fifelse(dataset == "GSE296710_zanovello", "GSE296710\nZanovello", "GSE296714\nHumphrey"),
  levels = c("GSE296710\nZanovello", "GSE296714\nHumphrey")
)]
cell_levels <- c("EX", "MG", "IN", "OPC", "ODC", "ASC")
cell_labels <- c(EX = "Excitatory N", MG = "Microglia", IN = "Inhibitory N", OPC = "OPC", ODC = "Oligodendrocyte", ASC = "Astrocyte")
perturbation[, cell_type := factor(cell_type, levels = rev(cell_levels))]

p_c <- ggplot(perturbation, aes(dataset_label, cell_type, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.65) +
  geom_text(
    aes(label = sprintf("%.2f%s", rho, ifelse(p < 0.05, "*", ""))),
    size = 2.05, colour = pal[["ink"]]
  ) +
  scale_fill_gradient2(
    low = pal[["red"]], mid = "white", high = pal[["blue"]], midpoint = 0,
    limits = c(-0.75, 0.75), name = "Spearman ρ"
  ) +
  scale_y_discrete(labels = cell_labels) +
  labs(
    title = "Perturbation–pathology\nalignment",
    subtitle = "Two TDP-43-loss lines; 6/6 glial correlations negative",
    x = NULL, y = NULL
  ) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    legend.position = "bottom", legend.key.width = grid::unit(10, "mm")
  )

# d. Chaperone-arm response -------------------------------------------------
chaperone <- fread(paths$chaperone)
chaperone[, dataset_label := factor(
  fifelse(dataset == "GSE296710_zanovello", "GSE296710", "GSE296714"),
  levels = c("GSE296710", "GSE296714")
)]
chaperone[, significant := factor(FDR < 0.05, levels = c(FALSE, TRUE), labels = c("not significant", "FDR < 0.05"))]
gene_order <- chaperone[, .(mean_effect = mean(KD_log2FC)), by = symbol][order(mean_effect), symbol]
chaperone[, symbol := factor(symbol, levels = gene_order)]

p_d <- ggplot(chaperone, aes(KD_log2FC, symbol, fill = dataset_label, alpha = significant)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = pal[["muted"]]) +
  geom_col(position = position_dodge(width = 0.78), width = 0.68, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("GSE296710" = pal[["blue"]], "GSE296714" = pal[["amber"]]), name = NULL) +
  scale_alpha_manual(values = c("not significant" = 0.30, "FDR < 0.05" = 1.0), name = NULL) +
  scale_x_continuous(limits = c(-1.8, 0.75), breaks = c(-1.5, -1, -0.5, 0, 0.5)) +
  labs(
    title = "Chaperone-arm response to TDP-43 loss",
    subtitle = "Knockdown versus control; colour denotes dataset",
    x = "RNA log2 fold change", y = NULL
  ) +
  guides(alpha = guide_legend(order = 2), fill = guide_legend(order = 1)) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 5.7))

# e. Competitive programme-set GWAS tests ---------------------------------
gwas_sets <- fread(paths$gwas_sets)
set_name <- function(x) {
  x <- sub("^FAM_pool3_transferable$", "Transferable three-family pool", x)
  x <- sub("^CAND17_frozen$", "17 frozen candidates", x)
  x <- sub("^LE_strict13_pooled$", "Strict-13 leading edge", x)
  x <- sub("^FAM_", "", x)
  x <- gsub("_", " ", x)
  vapply(x, function(value) paste(strwrap(value, width = 17), collapse = "\n"), character(1))
}
gwas_sets[, `:=`(set_label = set_name(VARIABLE), neglogP = -log10(P))]
setorder(gwas_sets, neglogP)
gwas_sets[, set_label := factor(set_label, levels = set_label)]

p_e <- ggplot(gwas_sets, aes(neglogP, set_label)) +
  geom_vline(xintercept = -log10(0.05), linetype = "22", linewidth = 0.4, colour = pal[["muted"]]) +
  geom_col(width = 0.60, fill = pal[["blue"]]) +
  geom_text(aes(label = sprintf("P = %.2f", P)), hjust = -0.12, size = 1.72, colour = pal[["ink"]]) +
  scale_x_continuous(limits = c(0, 1.48), breaks = c(0, 0.5, 1.0), expand = c(0, 0)) +
  labs(
    title = "ALS GWAS programme-set tests",
    subtitle = "Nine prespecified competitive tests; none passes P < 0.05",
    x = "−log10 P", y = NULL
  ) +
  theme(axis.text.y = element_text(size = 5.35, lineheight = 0.92))

# f. Candidate-gene ALS association ----------------------------------------
gwas_genes <- fread(paths$gwas_genes)
gwas_genes[, neglogP := -log10(P)]
setorder(gwas_genes, neglogP)
gwas_genes[, GENE := factor(GENE, levels = GENE)]
bonferroni_x <- -log10(0.05 / nrow(gwas_genes))

p_f <- ggplot(gwas_genes, aes(neglogP, GENE)) +
  geom_vline(xintercept = bonferroni_x, linetype = "22", linewidth = 0.45, colour = pal[["red"]]) +
  geom_point(aes(colour = GENE == "MAPT"), size = 2.15) +
  geom_text(
    data = gwas_genes[GENE == "MAPT"],
    aes(label = "P = 4.8 × 10⁻⁴"), hjust = 1.15,
    size = 1.8, colour = pal[["red"]], fontface = "bold"
  ) +
  scale_colour_manual(values = c(`FALSE` = "#737B87", `TRUE` = pal[["red"]]), guide = "none") +
  scale_x_continuous(limits = c(0, 4.15), breaks = 0:4, expand = c(0, 0)) +
  annotate(
    "text", x = bonferroni_x + 0.04, y = 1.1, label = "0.05/17",
    angle = 90, hjust = 0, size = 1.65, colour = pal[["red"]]
  ) +
  labs(
    title = "Candidate-gene ALS\nassociation",
    subtitle = "17 candidates; Bonferroni: MAPT only",
    x = "−log10 P", y = NULL
  ) +
  theme(axis.text.y = element_text(size = 5.45))

top_row <- p_a + p_b + p_c + plot_layout(widths = c(0.95, 1.18, 0.92))
bottom_row <- p_d + p_e + p_f + plot_layout(widths = c(1.05, 1.05, 1.0))
figure <- top_row / bottom_row +
  plot_layout(heights = c(1.0, 1.12)) +
  plot_annotation(tag_levels = "a")

width_mm <- 190
height_mm <- 172
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
prefix <- file.path(output_root, "ExtendedDataFigure14")

svglite::svglite(paste0(prefix, ".svg"), width = width_in, height = height_in, bg = "white")
print(figure)
dev.off()

grDevices::cairo_pdf(
  paste0(prefix, ".pdf"), width = width_in, height = height_in,
  family = font_body, bg = "white"
)
print(figure)
dev.off()

ragg::agg_tiff(
  paste0(prefix, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, background = "white", compression = "lzw"
)
print(figure)
dev.off()

ragg::agg_png(
  paste0(prefix, "_preview.png"), width = width_in, height = height_in,
  units = "in", res = 300, background = "white"
)
print(figure)
dev.off()

write_clean_tsv(waterfall, file.path(source_output, "ExtendedDataFigure14_panel_a_source_data.tsv"))
write_clean_tsv(attainability, file.path(source_output, "ExtendedDataFigure14_panel_b_source_data.tsv"))
write_clean_tsv(perturbation, file.path(source_output, "ExtendedDataFigure14_panel_c_source_data.tsv"))
write_clean_tsv(chaperone, file.path(source_output, "ExtendedDataFigure14_panel_d_source_data.tsv"))
write_clean_tsv(gwas_sets, file.path(source_output, "ExtendedDataFigure14_panel_e_source_data.tsv"))
write_clean_tsv(gwas_genes, file.path(source_output, "ExtendedDataFigure14_panel_f_source_data.tsv"))

qa_lines <- c(
  "Extended Data Figure 14 — final figure contract and QA",
  "Core conclusion: three-cohort failure reflects both gate attrition and direction reversal; perturbation and GWAS are bounded orthogonal anchors.",
  "Backend: R only.",
  "Dimensions: 190 × 172 mm.",
  "Exports: editable SVG/PDF, 600-dpi LZW TIFF and 300-dpi PNG preview.",
  "Panel a uses horizontal gate labels to prevent overlap.",
  "Panels b and f use expanded limits and direct labels to prevent clipping.",
  "All quantitative values are copied from frozen source tables; no new inferential test is introduced."
)
writeLines(qa_lines, file.path(output_root, "ExtendedDataFigure14_QA.txt"))

cat("Rendered submission-ready Extended Data Figure 14 to:", output_root, "\n")
