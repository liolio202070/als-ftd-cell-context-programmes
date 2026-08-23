script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
project_lib <- Sys.getenv("ALS_FTD_R_LIBRARY", unset = "")
if (nzchar(project_lib) && dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(svglite)
  library(ragg)
})

package_source <- file.path(project_root, "data", "source_data", "main_figures")
detailed_exact_path <- file.path(project_root, "data", "figure_inputs", "Figure3", "exact_estimand_program_replication.tsv")
table3_path <- file.path(project_root, "data", "figure_inputs", "Figure3", "Table3_strict_replicated_and_discordant_programme_axes.tsv")
draft_root <- file.path(project_root, "outputs", "figures", "Figure3")
draft_source <- file.path(draft_root, "source_data")
dir.create(draft_source, recursive = TRUE, showWarnings = FALSE)

frozen_paths <- c(
  detailed_exact_path,
  table3_path,
  file.path(package_source, sprintf("Figure3_panel_%s_source_data.tsv", c("a", "b", "c", "d", "e")))
)
if (!all(file.exists(frozen_paths))) {
  stop("One or more frozen Figure 3 source files are missing.")
}
file.copy(frozen_paths, draft_source, overwrite = TRUE)

exact <- fread(detailed_exact_path)
table3 <- fread(table3_path)
panel_b <- fread(file.path(package_source, "Figure3_panel_b_source_data.tsv"))
panel_c <- fread(file.path(package_source, "Figure3_panel_c_source_data.tsv"))
panel_d <- fread(file.path(package_source, "Figure3_panel_d_source_data.tsv"))
panel_e <- fread(file.path(package_source, "Figure3_panel_e_source_data.tsv"))

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
  neutral = "#AAB4C2",
  neutral_soft = "#EEF2F6",
  negative = "#3D6FB1",
  positive = "#B75E68"
)

font_body <- "Arial"
font_display <- "Avenir Next"

theme_fig3 <- function(base_size = 7.2) {
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
      legend.title = element_text(size = 6.4, face = "bold"),
      legend.text = element_text(size = 6.0),
      legend.key.height = grid::unit(3.2, "mm"),
      plot.title = element_text(
        family = font_display, face = "bold", size = 8.5,
        colour = pal[["ink"]], margin = margin(b = 2.5)
      ),
      plot.subtitle = element_text(
        size = 6.35, colour = pal[["muted"]], lineheight = 1.05,
        margin = margin(b = 4.5)
      ),
      plot.tag = element_text(
        family = font_display, face = "bold", size = 9.5,
        colour = pal[["ink"]]
      ),
      plot.tag.position = c(0, 1),
      plot.margin = margin(5, 6, 5, 6, unit = "pt")
    )
}

theme_set(theme_fig3())

write_clean_tsv <- function(x, path) {
  y <- copy(as.data.table(x))
  text_columns <- names(y)[vapply(y, function(z) is.character(z) || is.factor(z), logical(1))]
  for (column in text_columns) {
    value <- gsub("[\r\n\t]", " ", as.character(y[[column]]), perl = TRUE)
    set(y, j = column, value = value)
  }
  fwrite(y, path, sep = "\t", quote = FALSE, na = "")
}

estimand_levels <- c("ALS_vs_Control", "FTD_vs_Control", "ALS_vs_FTD")
estimand_labels <- c(
  "ALS_vs_Control" = "ALS vs control",
  "FTD_vs_Control" = "FTD vs control",
  "ALS_vs_FTD" = "ALS vs FTD\n(partition-bridged DiD)"
)
cell_shapes <- c(Astrocyte = 21, OPC = 22, Oligodendrocyte = 24)

# Panel a — all 63 matched axes are visible. Strict and stable-discordant axes are highlighted.
exact[, strict_replicated :=
        same_direction & both_global_FDR_0_05 &
        loo_direction_fraction_discovery >= 0.90 &
        loo_direction_fraction_replication >= 0.90]
exact[, stable_discordant :=
        !same_direction & both_global_FDR_0_05 &
        loo_direction_fraction_discovery >= 0.90 &
        loo_direction_fraction_replication >= 0.90]
exact[, evidence_class := fifelse(
  strict_replicated, "Strict replicated",
  fifelse(stable_discordant, "Stable discordant", "Other matched axis")
)]
exact[, evidence_class := factor(
  evidence_class,
  levels = c("Other matched axis", "Strict replicated", "Stable discordant")
)]
exact[, estimand_label := factor(
  unname(estimand_labels[estimand]), levels = unname(estimand_labels[estimand_levels])
)]
exact[, cell_type := factor(cell_type, levels = names(cell_shapes))]
exact[, discord_label := fifelse(
  stable_discordant,
  paste0(
    fifelse(
      cell_type == "Astrocyte", "Astro",
      fifelse(cell_type == "Oligodendrocyte", "Oligo", as.character(cell_type))
    ),
    " · ", family_label
  ),
  NA_character_
)]
setorder(exact, evidence_class)

evidence_fill <- c(
  "Other matched axis" = pal[["neutral_soft"]],
  "Strict replicated" = pal[["teal"]],
  "Stable discordant" = pal[["amber"]]
)
evidence_colour <- c(
  "Other matched axis" = pal[["neutral"]],
  "Strict replicated" = "#176C67",
  "Stable discordant" = "#9B5D13"
)

p_a <- ggplot(exact, aes(effect_discovery, effect_replication)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = pal[["teal_soft"]], alpha = 0.40) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = pal[["teal_soft"]], alpha = 0.40) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = pal[["amber_soft"]], alpha = 0.34) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = pal[["amber_soft"]], alpha = 0.34) +
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.34, colour = pal[["muted"]]) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.34, colour = pal[["muted"]]) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.35, colour = "white") +
  geom_point(
    aes(shape = cell_type, fill = evidence_class, colour = evidence_class),
    size = 2.5, stroke = 0.62
  ) +
  ggrepel::geom_text_repel(
    data = exact[stable_discordant == TRUE], aes(label = discord_label),
    family = font_body, size = 2.0, colour = pal[["ink"]],
    box.padding = 0.30, point.padding = 0.25,
    min.segment.length = 0, segment.size = 0.28,
    segment.colour = pal[["muted"]], max.overlaps = Inf,
    seed = 31415, show.legend = FALSE
  ) +
  facet_wrap(~estimand_label, nrow = 1) +
  scale_shape_manual(values = cell_shapes, name = "Cell type") +
  scale_fill_manual(values = evidence_fill, name = "Evidence class") +
  scale_colour_manual(values = evidence_colour, guide = "none") +
  scale_x_continuous(breaks = -2:3, limits = c(-2.70, 3.38), expand = c(0, 0)) +
  scale_y_continuous(breaks = -2:3, limits = c(-2.70, 3.38), expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title = "Exact-estimand programme replication",
    subtitle = paste0(
      "All 63 matched axes: 39 share direction, 13 meet the strict dual-cohort gate, and 3 are significant and LOO-stable but discordant.\n",
      "Pale teal quadrants indicate direction agreement; pale amber quadrants indicate disagreement."
    ),
    x = "Discovery cohort NES", y = "External C9 cohort NES", tag = "a"
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(shape = 21, size = 2.8)),
    shape = guide_legend(order = 2, override.aes = list(fill = "white", colour = pal[["ink"]], size = 2.6))
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = pal[["grid"]], linewidth = 0.3),
    panel.grid = element_blank(),
    axis.text = element_text(size = 6.0),
    strip.text = element_text(size = 6.9, lineheight = 0.95, margin = margin(b = 3)),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 5.7),
    legend.title = element_text(size = 6.1, face = "bold"),
    legend.spacing.x = grid::unit(1.2, "mm"),
    legend.margin = margin(t = -1),
    plot.tag.position = c(0.18, 1),
    plot.margin = margin(6, 80, 0, 8, unit = "pt")
  )

# Panel b — broader external-cohort programme landscape across all profiled cell types.
programme_levels <- c(
  "Axon development", "Heat response", "Mitochondrial energy",
  "Protein folding", "Protein stabilization", "Small GTPase", "Synapse/junction"
)
gittings_contrasts <- c(
  "C9_ALS_vs_Ctrl_Frontal", "C9_FTD_vs_Ctrl_Frontal",
  "C9_ALS_vs_C9_FTD_Frontal_batch_bridged"
)
gittings_labels <- c(
  "C9_ALS_vs_Ctrl_Frontal" = "C9-ALS\nvs control",
  "C9_FTD_vs_Ctrl_Frontal" = "C9-FTD\nvs control",
  "C9_ALS_vs_C9_FTD_Frontal_batch_bridged" = "ALS–FTD\nbridged DiD"
)
b_grid <- CJ(contrast = gittings_contrasts, set_label = programme_levels, unique = TRUE)
b_plot <- merge(
  b_grid,
  panel_b[, .(contrast, set_label, significant_axes, consensus_direction, best_global_FDR)],
  by = c("contrast", "set_label"), all.x = TRUE
)
b_plot[is.na(significant_axes), `:=`(significant_axes = 0L, consensus_direction = "none")]
b_plot[, contrast_label := factor(
  unname(gittings_labels[contrast]), levels = unname(gittings_labels[gittings_contrasts])
)]
b_plot[, programme_label := factor(set_label, levels = rev(programme_levels))]
b_plot[, number_colour := fifelse(consensus_direction == "none", pal[["muted"]], "white")]

p_b <- ggplot(b_plot, aes(contrast_label, programme_label)) +
  geom_point(
    aes(size = significant_axes, fill = consensus_direction),
    shape = 21, stroke = 0.45, colour = "white"
  ) +
  geom_text(
    aes(label = fifelse(significant_axes > 0, as.character(significant_axes), "0"),
        colour = number_colour),
    size = 1.85, fontface = "bold"
  ) +
  scale_colour_identity() +
  scale_size_continuous(range = c(1.4, 5.8), limits = c(0, 6), guide = "none") +
  scale_fill_manual(
    values = c(
      positive = pal[["positive"]], negative = pal[["negative"]],
      mixed = pal[["muted"]], none = "white"
    ),
    breaks = c("positive", "negative", "mixed", "none"),
    labels = c("positive", "negative", "mixed", "no significant axis"),
    name = "Consensus direction"
  ) +
  labs(
    title = "External C9 programme landscape",
    subtitle = "Global-FDR cell-type counts;\ncolour gives consensus direction.",
    x = NULL, y = NULL, tag = "b"
  ) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 5.5, lineheight = 0.92, margin = margin(t = 2)),
    axis.text.y = element_text(size = 5.55, lineheight = 0.95, margin = margin(r = 2)),
    panel.grid.major = element_line(colour = pal[["grid"]], linewidth = 0.30),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.key.width = grid::unit(3.1, "mm"),
    legend.spacing.x = grid::unit(1.3, "mm"),
    plot.tag.position = c(0.35, 1),
    plot.margin = margin(6, 6, 5, 8, unit = "pt")
  )

# Panel c — compact, exact donor-omission robustness summary.
robust <- data.table(
  metric = factor(c("≥90% direction retained", "100% direction retained"),
                  levels = rev(c("≥90% direction retained", "100% direction retained"))),
  fraction = c(68 / 68, 66 / 68),
  label = c("68/68", "66/68"),
  fill = c(pal[["teal"]], pal[["blue"]])
)
p_c <- ggplot(robust, aes(fraction, metric)) +
  geom_segment(aes(x = 0.90, xend = 1.00, yend = metric),
               linewidth = 3.8, colour = pal[["neutral_soft"]]) +
  geom_segment(aes(x = 0.90, xend = fraction, yend = metric, colour = fill),
               linewidth = 3.8) +
  geom_point(aes(colour = fill), size = 2.2) +
  geom_text(aes(label = label), hjust = -0.12, size = 2.25, fontface = "bold") +
  scale_colour_identity() +
  scale_x_continuous(
    limits = c(0.90, 1.035), breaks = c(0.90, 0.95, 1.00),
    labels = percent_format(accuracy = 1), expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Donor-LOO robustness",
    subtitle = "68 global-FDR axes;\n65 had 37 valid refits, 3 had 36.",
    x = "Donor-LOO direction agreement", y = NULL, tag = "c"
  ) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 5.8),
    axis.text.x = element_text(size = 5.8),
    panel.grid.major.x = element_line(colour = pal[["grid"]], linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.tag.position = c(0.08, 1),
    plot.margin = margin(6, 8, 6, 8, unit = "pt")
  )

# Panel d — cross-genetic boundary support without invented intervals.
panel_d[, axis_label := fifelse(
  population == "NEUNpos", "NEUN+", "OLIG2+"
)]
panel_d[, axis_label := paste(axis_label, set_label, sep = " · ")]
panel_d[, axis_label := factor(axis_label, levels = axis_label[order(NES)])]
panel_d[, direction := fifelse(NES > 0, "positive", "negative")]

p_d <- ggplot(panel_d, aes(NES, axis_label, colour = direction)) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.36, colour = pal[["muted"]]) +
  geom_segment(aes(x = 0, xend = NES, yend = axis_label), linewidth = 0.70) +
  geom_point(size = 2.65) +
  scale_colour_manual(values = c(positive = pal[["mauve"]], negative = pal[["blue"]]), guide = "none") +
  scale_x_continuous(limits = c(-2.05, 2.15), breaks = -2:2, expand = c(0, 0)) +
  labs(
    title = "GRN–FTD boundary",
    subtitle = "Three axes; each 9/9 LOO-stable. NES shown without CI.",
    x = "NES", y = NULL, tag = "d"
  ) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 5.35, lineheight = 0.92),
    panel.grid.major.x = element_line(colour = pal[["grid"]], linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.tag.position = c(0.33, 1),
    plot.margin = margin(6, 7, 6, 8, unit = "pt")
  )

# Panel e — retain the gene-level negative boundary without claiming replicated biomarkers.
gene_summary <- panel_e[, .(
  axes = .N,
  same_direction_n = sum(same_direction, na.rm = TRUE),
  same_direction_fraction = mean(same_direction, na.rm = TRUE),
  strict_axes = sum(
    same_direction & discovery_adjusted_p < 0.05 & replication_adjusted_p < 0.05,
    na.rm = TRUE
  )
), by = gene]
setorder(gene_summary, -same_direction_fraction, gene)
gene_summary[, rank_id := .I]
gene_summary[, block := fifelse(rank_id <= 9, "Top-ranked genes", "Remaining genes")]
gene_summary[, block := factor(block, levels = c("Top-ranked genes", "Remaining genes"))]
gene_summary[, gene := factor(gene, levels = rev(as.character(gene)))]
gene_summary[, label := sprintf("%d/9", same_direction_n)]

p_e <- ggplot(gene_summary, aes(same_direction_fraction, gene)) +
  geom_segment(aes(x = 0, xend = same_direction_fraction, yend = gene),
               linewidth = 1.55, colour = pal[["blue_soft"]]) +
  geom_point(size = 2.0, colour = pal[["blue"]]) +
  geom_text(aes(x = 1.025, label = label), hjust = 0, size = 1.9, colour = pal[["ink"]]) +
  scale_x_continuous(
    limits = c(0, 1.14), breaks = c(0, 0.5, 1),
    labels = percent_format(accuracy = 1), expand = c(0, 0)
  ) +
  labs(
    title = "Candidate-gene boundary",
    subtitle = "Direction agreement across 9 settings per gene.\nNo axis passes the strict gene gate (0/153).",
    x = "Same-direction matched axes", y = NULL, tag = "e"
  ) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 5.25),
    axis.text.x = element_text(size = 5.8),
    panel.grid.major.x = element_line(colour = pal[["grid"]], linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.tag.position = c(0.03, 1),
    plot.margin = margin(6, 6, 6, 24, unit = "pt")
  )

middle_row <- p_b | p_c
middle_row <- middle_row + plot_layout(widths = c(1, 1))
bottom_row <- p_d | p_e
bottom_row <- bottom_row + plot_layout(widths = c(1, 1))
figure3 <- p_a / middle_row / bottom_row +
  plot_layout(heights = c(1.08, 0.78, 1.18))

width_mm <- 190
height_mm <- 210
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
prefix <- file.path(draft_root, "Figure3_redesign_draft")

svglite::svglite(paste0(prefix, ".svg"), width = width_in, height = height_in, bg = "white")
print(figure3)
dev.off()

grDevices::cairo_pdf(
  paste0(prefix, ".pdf"), width = width_in, height = height_in,
  family = font_body, bg = "white"
)
print(figure3)
dev.off()

ragg::agg_tiff(
  paste0(prefix, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, background = "white", compression = "lzw"
)
print(figure3)
dev.off()

ragg::agg_png(
  paste0(prefix, "_preview.png"), width = width_in, height = height_in,
  units = "in", res = 300, background = "white"
)
print(figure3)
dev.off()

write_clean_tsv(exact, file.path(draft_source, "Figure3_panel_a_all_63_exact_axes.tsv"))
write_clean_tsv(b_plot, file.path(draft_source, "Figure3_panel_b_complete_landscape.tsv"))
write_clean_tsv(robust, file.path(draft_source, "Figure3_panel_c_robustness_summary.tsv"))
write_clean_tsv(gene_summary, file.path(draft_source, "Figure3_panel_e_gene_summary.tsv"))

qa_lines <- c(
  "Figure 3 redesign draft — statistical and interpretive boundaries",
  "Core conclusion: exact programme axes replicate more consistently than individual genes, with explicit stable heterogeneity and genetic-background boundaries.",
  "Panel a: all 63 exactly matched axes are shown. 39 are directionally concordant; 13 meet the prespecified dual-cohort FDR plus >=90% donor-LOO gate; 3 are significant and LOO-stable in both cohorts but discordant.",
  "Panel a effects are cohort-specific FGSEA NES values and are not pooled.",
  "Panel b: globally significant Gittings cell-type axes are counted by programme and contrast. The ALS-versus-FTD contrast is a partition-bridged difference-in-differences and is secondary support.",
  "Panel c: all 68 globally significant Gittings axes retain at least 90% LOO direction agreement; 66/68 retain 100%. LOO is internal robustness, not another cohort.",
  "Panel d: three 9-donor GRN-FTD axes are shown as NES points without fabricated uncertainty intervals; this is cross-genetic support, not direct C9 replication.",
  "Panel e: none of 153 candidate-gene axes passes the strict dual-cohort gene gate; the 17 genes remain discovery-derived contributors/hypotheses rather than replicated biomarkers.",
  "Public identifiers show no detected overlap between discovery and Gittings, but de-identification precludes proof of complete donor independence."
)
writeLines(qa_lines, file.path(draft_root, "Figure3_redesign_draft_QA.txt"))

cat("Rendered Figure 3 draft to:", draft_root, "\n")
cat("Exact axes:", nrow(exact), "\n")
cat("Same-direction axes:", sum(exact$same_direction), "\n")
cat("Strict replicated axes:", sum(exact$strict_replicated), "\n")
cat("Stable discordant axes:", sum(exact$stable_discordant), "\n")
cat("Globally significant Gittings axes in LOO panel:", nrow(panel_c), "\n")
cat("Strict dual-cohort candidate-gene axes:", sum(gene_summary$strict_axes), "\n")
