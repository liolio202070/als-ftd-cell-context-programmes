script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
project_lib <- Sys.getenv("ALS_FTD_R_LIBRARY", unset = "")
if (nzchar(project_lib) && dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
})

source_root <- file.path(project_root, "data", "source_data", "main_figures")
draft_root <- file.path(project_root, "outputs", "figures", "Figure5")
draft_source <- file.path(draft_root, "source_data")
dir.create(draft_source, recursive = TRUE, showWarnings = FALSE)

source_paths <- file.path(
  source_root, sprintf("Figure5_panel_%s_source_data.tsv", letters[1:5])
)
if (!all(file.exists(source_paths))) {
  stop("One or more frozen Figure 5 source files are missing.")
}
file.copy(source_paths, draft_source, overwrite = TRUE)

panel_a <- fread(source_paths[1])
panel_b <- fread(source_paths[2])
panel_c <- fread(source_paths[3])
panel_d <- fread(source_paths[4])
panel_e <- fread(source_paths[5])

pal <- c(
  ink = "#172033",
  muted = "#667085",
  grid = "#D9E0E8",
  paper = "#FFFFFF",
  soft = "#F4F7FA",
  blue = "#2867A8",
  blue_soft = "#DDE9F5",
  teal = "#278F88",
  teal_soft = "#DCEFEB",
  mauve = "#A85674",
  mauve_soft = "#F1DEE6",
  amber = "#CF872C",
  amber_soft = "#F7E9D5",
  violet = "#6955A5",
  violet_soft = "#E9E4F3",
  neutral = "#AAB4C2",
  neutral_soft = "#EEF2F6",
  negative = "#3D6FB1",
  positive = "#B75E68"
)

font_body <- "Arial"
font_display <- "Avenir Next"

theme_fig5 <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = font_body) +
    theme(
      text = element_text(colour = pal[["ink"]]),
      axis.line = element_line(linewidth = 0.35, colour = pal[["ink"]]),
      axis.ticks = element_line(linewidth = 0.35, colour = pal[["ink"]]),
      axis.ticks.length = grid::unit(1.2, "mm"),
      axis.title = element_text(size = 7.0, colour = pal[["ink"]]),
      axis.text = element_text(size = 6.4, colour = pal[["ink"]]),
      strip.background = element_blank(),
      strip.text = element_text(
        family = font_display, size = 7.0, face = "bold", colour = pal[["ink"]]
      ),
      legend.title = element_text(size = 6.3, face = "bold"),
      legend.text = element_text(size = 6.0),
      legend.key.height = grid::unit(3.0, "mm"),
      plot.title = element_text(
        family = font_display, face = "bold", size = 8.5,
        colour = pal[["ink"]], margin = margin(b = 2.5)
      ),
      plot.subtitle = element_text(
        size = 6.25, colour = pal[["muted"]], lineheight = 1.05,
        margin = margin(b = 4.5)
      ),
      plot.tag = element_text(
        family = font_display, face = "bold", size = 9.5,
        colour = pal[["ink"]]
      ),
      plot.tag.position = c(0, 1),
      plot.margin = margin(6, 7, 6, 7, unit = "pt")
    )
}

theme_set(theme_fig5())

write_clean_tsv <- function(x, path) {
  y <- copy(as.data.table(x))
  text_columns <- names(y)[vapply(y, function(z) is.character(z) || is.factor(z), logical(1))]
  for (column in text_columns) {
    value <- gsub("[\r\n\t]", " ", as.character(y[[column]]), perl = TRUE)
    set(y, j = column, value = value)
  }
  fwrite(y, path, sep = "\t", quote = FALSE, na = "")
}

# Figure contract:
# Core conclusion: ordinal pTDP burden is associated with an astrocyte proteostasis
# decline and predominantly oligodendroglial accessibility differences, but no
# candidate passes a complete RNA–ATAC disease-trend chain.
# Archetype: quantitative grid with two RNA panels, two chromatin-context panels,
# and a full-width evidence-boundary row.
# Backend: R only. Cell-resolved inferential unit: 14 annotated Emory donors.

cell_colours <- c(
  ASC = pal[["mauve"]],
  MG = pal[["teal"]],
  OPC = pal[["blue"]],
  ODC = pal[["violet"]]
)

# Panel a — programme-level primary signal; only this frozen panel has a defensible CI.
family_levels <- c(
  "Mitochondrial energy", "Protein stabilization", "Protein folding",
  "Heat response", "Small GTPase", "Synapse / junction assembly", "Axon development"
)
a_plot <- panel_a[cell_type == "ASC"]
a_plot[, family := factor(family, levels = family_levels)]
a_plot[, highlighted := family == "Protein stabilization"]
a_signal <- a_plot[highlighted == TRUE]
a_signal[, signal_label := sprintf(
  "β %.2f  [%.2f, %.2f]\nglobal BH q %.3f",
  beta_per_pTDP_level, ci_low, ci_high, global_FDR
)]

p_a <- ggplot(a_plot, aes(beta_per_pTDP_level, family)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45, colour = pal[["neutral"]]) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high, colour = highlighted),
    orientation = "y", width = 0.15, linewidth = 0.55
  ) +
  geom_point(aes(fill = highlighted), shape = 21, size = 2.5, stroke = 0.55, colour = "white") +
  geom_label(
    data = a_signal,
    aes(x = -3.28, y = family, label = signal_label),
    inherit.aes = FALSE, hjust = 0, size = 1.82, lineheight = 1.02,
    linewidth = 0, fill = pal[["mauve_soft"]], colour = pal[["mauve"]]
  ) +
  scale_colour_manual(values = c(`FALSE` = pal[["neutral"]], `TRUE` = pal[["mauve"]]), guide = "none") +
  scale_fill_manual(values = c(`FALSE` = pal[["neutral"]], `TRUE` = pal[["mauve"]]), guide = "none") +
  scale_x_continuous(breaks = c(-3, -2, -1, 0, 1)) +
  coord_cartesian(xlim = c(-3.35, 1.12), clip = "off") +
  labs(
    x = "RNA programme β per pTDP level (normal-approximation 95% CI)", y = NULL,
    title = "Astrocyte proteostasis declines with pTDP burden",
    subtitle = "14 annotated Emory donors; global BH across 42 programme-by-cell-type tests",
    tag = "a"
  ) +
  theme(
    axis.text.y = element_text(size = 6.25),
    plot.tag.position = c(0.20, 1),
    plot.margin = margin(6, 12, 6, 18, unit = "pt")
  )

# Panel b — eight candidate-set RNA signals, all direction-stable in 14/14 donor omissions.
b_hits <- panel_b[cell_type != "whole_tissue" & target_global_FDR < 0.05]
b_hits[, display_cell := fcase(
  cell_type == "ASC", "Astrocyte",
  cell_type == "MG", "Microglia",
  cell_type == "OPC", "OPC",
  default = cell_type
)]
b_hits[, axis_label := paste(display_cell, feature, sep = "  ·  ")]
b_hits[, transcriptome_wide := FDR < 0.05]
b_hits[, axis_label := factor(axis_label, levels = axis_label[order(log2FC_per_pTDP_level)])]

p_b <- ggplot(b_hits, aes(log2FC_per_pTDP_level, axis_label)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45, colour = pal[["neutral"]]) +
  geom_segment(
    aes(x = 0, xend = log2FC_per_pTDP_level, yend = axis_label),
    colour = pal[["grid"]], linewidth = 1.1
  ) +
  geom_point(
    aes(fill = cell_type, shape = transcriptome_wide),
    size = 2.8, stroke = 0.75, colour = "white"
  ) +
  scale_fill_manual(
    values = cell_colours[c("ASC", "MG", "OPC")],
    labels = c(ASC = "Astrocyte", MG = "Microglia", OPC = "OPC"),
    guide = guide_legend(title = NULL)
  ) +
  scale_shape_manual(
    values = c(`FALSE` = 21, `TRUE` = 23),
    labels = c(`FALSE` = "candidate-set only", `TRUE` = "also transcriptome-wide"),
    name = "FDR layer"
  ) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1)) +
  coord_cartesian(xlim = c(-2.45, 1.10), clip = "off") +
  labs(
    x = "RNA log2FC per pTDP level", y = NULL,
    title = "Candidate RNA trends",
    subtitle = "14 donors; 8 candidate-set signals, all 14/14 LOO;\nBAG3 alone passes transcriptome-wide FDR",
    tag = "b"
  ) +
  theme(
    axis.text.y = element_text(size = 6.1),
    legend.position = "bottom", legend.box = "vertical",
    legend.margin = margin(t = -2), legend.key.width = grid::unit(4, "mm"),
    plot.tag.position = c(0.24, 1),
    plot.margin = margin(6, 10, 6, 10, unit = "pt")
  ) +
  guides(
    fill = guide_legend(
      title = "Cell type", order = 1,
      override.aes = list(shape = 21, size = 2.8, stroke = 0.6, colour = "white")
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(fill = pal[["neutral"]], size = 2.8, stroke = 0.6, colour = "white")
    )
  )

# Panel c — detected DAR count summary; counts do not imply selective causality.
c_summary <- panel_c[, .(
  total = .N,
  loss = sum(log2FC_per_pTDP_level < 0),
  gain = sum(log2FC_per_pTDP_level > 0),
  n_donors = unique(n_donors),
  peak_test_universe = unique(peak_test_universe)
), by = cell_type]
c_summary[, display_cell := fifelse(cell_type == "ODC", "Oligodendrocyte", "Microglia")]
c_long <- melt(
  c_summary,
  id.vars = c("cell_type", "display_cell", "total", "n_donors", "peak_test_universe"),
  measure.vars = c("loss", "gain"),
  variable.name = "direction", value.name = "n"
)
c_long[, direction := factor(direction, levels = c("loss", "gain"), labels = c("accessibility loss", "accessibility gain"))]
c_long[, display_cell := factor(display_cell, levels = c("Microglia", "Oligodendrocyte"))]

p_c <- ggplot(c_long, aes(n, display_cell, fill = direction)) +
  geom_col(width = 0.58, position = position_stack(reverse = TRUE)) +
  geom_text(
    data = c_summary,
    aes(x = total + 7, y = display_cell, label = total),
    inherit.aes = FALSE, hjust = 0, size = 2.35, fontface = "bold", colour = pal[["ink"]]
  ) +
  scale_fill_manual(values = c("accessibility loss" = pal[["blue"]], "accessibility gain" = pal[["teal_soft"]])) +
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(xlim = c(0, 455), clip = "off") +
  labs(
    x = "Genome-wide pTDP-associated DARs", y = NULL, fill = NULL,
    title = "Detected accessibility differences are oligodendroglial",
    subtitle = "14 donors; 406 of 409 DARs occur in oligodendrocytes under the fitted model",
    tag = "c"
  ) +
  theme(
    axis.text.y = element_text(size = 6.35),
    legend.position = "bottom", legend.direction = "horizontal",
    plot.tag.position = c(0.20, 1),
    plot.margin = margin(6, 12, 6, 18, unit = "pt")
  )

# Panel d — adjusted donor-level residual correlations, shown without invented CI.
d_plot <- copy(panel_d)
d_plot[, distance_label := fifelse(
  distance_to_gene == 0,
  "promoter",
  fifelse(distance_to_gene < 1000, sprintf("%.1f kb", distance_to_gene / 1000), sprintf("%.0f kb", distance_to_gene / 1000))
)]
d_plot[, axis_label := paste(gene, distance_label, sep = "  ·  ")]
d_plot[, axis_label := factor(axis_label, levels = axis_label[order(rho)])]
d_plot[, direction := ifelse(rho < 0, "negative", "positive")]
d_plot[, q_label := sprintf("q %.3f", FDR)]
d_plot[, label_x := rho + ifelse(rho < 0, -0.045, 0.045)]
d_plot[, label_hjust := ifelse(rho < 0, 1, 0)]

p_d <- ggplot(d_plot, aes(rho, axis_label, colour = direction)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45, colour = pal[["neutral"]]) +
  geom_segment(aes(x = 0, xend = rho, yend = axis_label), linewidth = 1.2, alpha = 0.22) +
  geom_point(size = 2.7) +
  geom_text(
    aes(x = label_x, label = q_label, hjust = label_hjust),
    size = 1.8, fontface = "bold", show.legend = FALSE
  ) +
  scale_colour_manual(values = c(negative = pal[["blue"]], positive = pal[["amber"]]), guide = "none") +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) +
  coord_cartesian(xlim = c(-1.18, 1.18), clip = "off") +
  labs(
    x = "Adjusted donor-level peak–gene Spearman correlation (ρ)", y = NULL,
    title = "Peak–gene associations",
    subtitle = "14 donors; 7 within-cell-type BH hits;\npoint estimates only (no frozen CI)",
    tag = "d"
  ) +
  theme(
    axis.text.y = element_text(size = 6.0),
    plot.tag.position = c(0.24, 1),
    plot.margin = margin(6, 16, 6, 10, unit = "pt")
  )

# Panel e — gene-level evidence components and the required null complete-chain result.
criterion_levels <- c(
  "RNA pTDP trend", "Peak–gene correlation",
  "Nearest-peak pTDP trend", "Full three-part chain"
)
criterion_labels <- c(
  "RNA pTDP\ntrend", "Peak–gene\nassociation",
  "Nearest-peak\npTDP trend", "Complete\nthree-part chain"
)
panel_e[, criterion := factor(criterion, levels = criterion_levels, labels = criterion_labels)]
gene_order_data <- panel_e[, .(
  evidence_count = sum(passed),
  has_rna = any(passed & as.character(criterion) == "RNA pTDP\ntrend"),
  has_link = any(passed & as.character(criterion) == "Peak–gene\nassociation")
), by = gene]
setorder(gene_order_data, -evidence_count, -has_rna, -has_link, gene)
panel_e[, gene := factor(gene, levels = rev(gene_order_data$gene))]
panel_e[, pass_class := fifelse(passed, as.character(criterion), "not supported")]

e_colours <- c(
  "not supported" = pal[["neutral_soft"]],
  "RNA pTDP\ntrend" = pal[["mauve"]],
  "Peak–gene\nassociation" = pal[["teal"]],
  "Nearest-peak\npTDP trend" = pal[["blue"]],
  "Complete\nthree-part chain" = pal[["ink"]]
)

p_e <- ggplot(panel_e, aes(criterion, gene)) +
  geom_tile(aes(fill = pass_class), colour = "white", linewidth = 0.55) +
  geom_text(aes(label = ifelse(passed, "•", "")), colour = "white", size = 3.2) +
  scale_fill_manual(values = e_colours, guide = "none") +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL, y = NULL, title = "Joint-evidence gate",
    subtitle = "14 donors; RNA trends: 6 genes · peak–gene associations: 5 genes\nNearest-peak trends: 0 · complete chains: 0",
    tag = "e"
  ) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 6.05, face = "bold", lineheight = 0.94, margin = margin(b = 3)),
    axis.text.y = element_text(size = 5.9, face = "bold", margin = margin(r = 3)),
    plot.subtitle = element_text(
      size = 6.0, colour = pal[["muted"]], lineheight = 1.05,
      margin = margin(b = 4.5)
    ),
    plot.tag.position = c(0.14, 1),
    plot.margin = margin(6, 18, 6, 18, unit = "pt")
  )

top_row <- p_a | p_b
top_row <- top_row + plot_layout(widths = c(1, 1))
middle_row <- p_c | p_d
middle_row <- middle_row + plot_layout(widths = c(1, 1))
figure5 <- top_row / middle_row / p_e +
  plot_layout(heights = c(1.0, 0.84, 1.16))

width_mm <- 190
height_mm <- 205
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
prefix <- file.path(draft_root, "Figure5_redesign_draft")

svglite::svglite(paste0(prefix, ".svg"), width = width_in, height = height_in, bg = "white")
print(figure5)
dev.off()

grDevices::cairo_pdf(
  paste0(prefix, ".pdf"), width = width_in, height = height_in,
  family = font_body, bg = "white"
)
print(figure5)
dev.off()

ragg::agg_tiff(
  paste0(prefix, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, background = "white", compression = "lzw"
)
print(figure5)
dev.off()

ragg::agg_png(
  paste0(prefix, "_preview.png"), width = width_in, height = height_in,
  units = "in", res = 300, background = "white"
)
print(figure5)
dev.off()

write_clean_tsv(b_hits, file.path(draft_source, "Figure5_panel_b_eight_candidate_RNA_hits.tsv"))
write_clean_tsv(c_summary, file.path(draft_source, "Figure5_panel_c_DAR_summary.tsv"))
write_clean_tsv(d_plot, file.path(draft_source, "Figure5_panel_d_seven_peak_gene_associations.tsv"))

qa_lines <- c(
  "Figure 5 redesign draft — statistical and interpretive boundaries",
  "Core conclusion: ordinal pTDP burden is associated with astrocyte proteostasis decline and predominantly oligodendroglial accessibility differences, but no candidate passes a complete RNA–ATAC disease-trend chain.",
  "All displayed cell-resolved inferential panels use 14 annotated Emory biological donors; the 26-donor full GEO series is not the cell-resolved inferential denominator.",
  "Panel a: OLS donor programme-score trend with technical covariates; normal-approximation 95% CI is available only for this panel; global BH spans 42 programme-by-cell-type tests.",
  "Panel b: eight candidate RNA signals pass BH across 119 frozen candidate tests including 17 whole-tissue rows; all retain direction in 14/14 donor omissions. BAG3 microglia alone also passes within-model transcriptome-wide FDR. No CI is available in the frozen edgeR export.",
  "Panel c: 409 pTDP-associated DARs pass BH separately within donor-cell-type peak models; 406 occur in oligodendrocytes and 3 in microglia. This detected distribution does not establish selective causal remodelling.",
  "Panel d: seven residualized donor-level Spearman associations pass BH across 495 tested target correlations within astrocytes. Points are shown without fabricated CIs and are not enhancer-target proof.",
  "Panel e: six genes have an RNA trend, five have a peak-gene association, no nearest candidate peak passes the pTDP trend gate, and no complete three-part chain passes.",
  "Motif/chromVAR interpretation is reported as a methodological boundary in the manuscript rather than as a data panel because a prespecified validated donor-level model with technical adjustment and LOO stability is unavailable.",
  "No ATAC-to-RNA causal, enhancer-target or TF-activity claim is made."
)
writeLines(qa_lines, file.path(draft_root, "Figure5_redesign_draft_QA.txt"))

cat("Rendered Figure 5 draft to:", draft_root, "\n")
cat("Programme tests:", nrow(panel_a), "\n")
cat("Candidate RNA hits:", nrow(b_hits), "\n")
cat("Transcriptome-wide candidate RNA hits:", sum(b_hits$transcriptome_wide), "\n")
cat("DARs:", nrow(panel_c), "\n")
cat("Peak-gene associations:", nrow(panel_d), "\n")
cat("Complete three-part chains:", sum(panel_e$passed[as.character(panel_e$criterion) == "Complete\nthree-part chain"]), "\n")
