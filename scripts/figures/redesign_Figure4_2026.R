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
draft_root <- file.path(project_root, "outputs", "figures", "Figure4")
draft_source <- file.path(draft_root, "source_data")
dir.create(draft_source, recursive = TRUE, showWarnings = FALSE)

source_paths <- file.path(
  source_root, sprintf("Figure4_panel_%s_source_data.tsv", letters[1:5])
)
if (!all(file.exists(source_paths))) {
  stop("One or more frozen Figure 4 source files are missing.")
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
  neutral_soft = "#EEF2F6"
)

font_body <- "Arial"
font_display <- "Avenir Next"

theme_fig4 <- function(base_size = 7.2) {
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
        size = 6.35, colour = pal[["muted"]], lineheight = 1.05,
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

theme_set(theme_fig4())

write_clean_tsv <- function(x, path) {
  y <- copy(as.data.table(x))
  text_columns <- names(y)[vapply(y, function(z) is.character(z) || is.factor(z), logical(1))]
  for (column in text_columns) {
    value <- gsub("[\r\n\t]", " ", as.character(y[[column]]), perl = TRUE)
    set(y, j = column, value = value)
  }
  fwrite(y, path, sep = "\t", quote = FALSE, na = "")
}

wrap_genes <- function(x, width = 25) {
  paste(strwrap(paste(sort(unique(x)), collapse = "  ·  "), width = width), collapse = "\n")
}

# Figure contract:
# Core conclusion: a finite 17-gene discovery-derived candidate set links programme
# evidence to external published TDP-43 event-level context, within explicit assay limits.
# Archetype: asymmetric quantitative composite with a full-width hero matrix.
# Backend: R only. No pooled evidence statistic or causal regulatory edge is shown.

# Panel a — full finite evidence matrix; counts remain layer-specific.
layer_levels <- c(
  "Original 7", "Same-cell programmes", "External FTD", "Regional snRNA/spatial",
  "Same-donor ATAC", "Proteomics", "External/cross-modal"
)
layer_labels <- c(
  "Discovery\ncomparisons", "Same-cell\nprogrammes", "External FTD\ncontext",
  "Regional snRNA /\nspatial context", "Same-donor\nATAC", "Proteomics",
  "External /\ncross-modal"
)
panel_a[, layer_label := factor(layer_label, levels = layer_levels, labels = layer_labels)]
gene_order <- unique(panel_a[order(-evidence_score, gene), gene])
panel_a[, gene := factor(gene, levels = rev(gene_order))]
panel_a[, display_value := ifelse(supported, as.character(value), "")]
panel_a[, candidate_tier := fifelse(grepl("Tier 1", evidence_tier), "Tier 1 candidate", "Tier 2 candidate")]

p_a <- ggplot(panel_a, aes(layer_label, gene)) +
  geom_tile(aes(fill = supported), colour = "white", linewidth = 0.55) +
  geom_text(
    aes(label = display_value, colour = supported),
    size = 2.25, fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(`FALSE` = pal[["neutral_soft"]], `TRUE` = pal[["teal"]]),
    labels = c(`FALSE` = "not supported", `TRUE` = "supported"),
    name = NULL
  ) +
  scale_colour_manual(values = c(`FALSE` = pal[["muted"]], `TRUE` = "white"), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL, y = NULL,
    title = "Finite candidate evidence matrix",
    subtitle = "17 discovery-derived candidates; numbers are supporting units within each evidence layer",
    tag = "a"
  ) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 6.4, face = "bold", lineheight = 0.95, margin = margin(b = 3)),
    axis.text.y = element_text(size = 6.5, face = "bold", margin = margin(r = 3)),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.margin = margin(t = -2), legend.box.margin = margin(0, 0, 0, 0),
    plot.tag.position = c(0.12, 1),
    plot.margin = margin(6, 8, 2, 20, unit = "pt")
  ) +
  guides(fill = guide_legend(override.aes = list(colour = NA), keywidth = grid::unit(5, "mm")))

# Panel b — collapse source programme labels into the four manuscript-level families.
panel_b[, programme_group := fcase(
  set_label == "Mitochondrial energy", "Mitochondrial",
  set_label %chin% c("Protein folding", "Protein stabilization", "Heat response"), "Proteostasis",
  set_label == "Small GTPase", "Small-GTPase",
  set_label %chin% c("Axon development", "Synapse / junction assembly"), "Synaptic / axonal",
  default = "Other"
)]
programme_cards <- panel_b[programme_group != "Other", .(
  genes = wrap_genes(gene),
  n_genes = uniqueN(gene)
), by = programme_group]
card_positions <- data.table(
  programme_group = c("Mitochondrial", "Proteostasis", "Small-GTPase", "Synaptic / axonal"),
  xmin = c(0.03, 1.03, 0.03, 1.03), xmax = c(0.97, 1.97, 0.97, 1.97),
  ymin = c(1.03, 1.03, 0.03, 0.03), ymax = c(1.97, 1.97, 0.97, 0.97),
  fill = c(pal[["teal_soft"]], pal[["amber_soft"]], pal[["mauve_soft"]], pal[["blue_soft"]]),
  accent = c(pal[["teal"]], pal[["amber"]], pal[["mauve"]], pal[["blue"]])
)
programme_cards <- merge(programme_cards, card_positions, by = "programme_group", all.x = TRUE)
programme_cards[, x := (xmin + xmax) / 2]
programme_cards[, title_y := ymax - 0.19]
programme_cards[, genes_y := (ymin + ymax) / 2 - 0.17]

p_b <- ggplot(programme_cards) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = I(fill)),
    colour = "white", linewidth = 0.8
  ) +
  geom_segment(
    aes(x = xmin + 0.06, xend = xmax - 0.06, y = ymax - 0.34, yend = ymax - 0.34, colour = I(accent)),
    linewidth = 0.6
  ) +
  geom_text(
    aes(x = x, y = title_y, label = sprintf("%s  ·  n=%d", programme_group, n_genes), colour = I(accent)),
    family = font_display, fontface = "bold", size = 2.45
  ) +
  geom_text(
    aes(x = x, y = genes_y, label = genes),
    colour = pal[["ink"]], size = 1.82, lineheight = 1.06
  ) +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 2), expand = FALSE, clip = "off") +
  labs(
    title = "Programme contributor map",
    subtitle = "Shared membership is undirected and does not imply regulation",
    tag = "b"
  ) +
  theme_void(base_family = font_body) +
  theme(
    plot.title = element_text(family = font_display, face = "bold", size = 8.5, colour = pal[["ink"]], margin = margin(b = 2.5)),
    plot.subtitle = element_text(size = 6.35, colour = pal[["muted"]], margin = margin(b = 4.5)),
    plot.tag = element_text(family = font_display, face = "bold", size = 9.5, colour = pal[["ink"]]),
    plot.tag.position = c(0.24, 1),
    plot.margin = margin(6, 10, 6, 18, unit = "pt")
  )

# Panel c — external, published event-level anchors only.
panel_c[, event_status := ifelse(
  gittings_snrna_event_detected,
  "Detected in published C9 snRNA",
  "Established literature target; not detected in C9 snRNA"
)]
panel_c[, gene := factor(gene, levels = c("STMN2", "KALRN", "UNC13A"))]
panel_c[, fdr_label := sprintf("ΔPSI %.3f  ·  FDR %.0e", event_strength, min_LeafCutter_FDR)]
panel_c[, fdr_label := gsub("e-0", "e−", fdr_label, fixed = TRUE)]
panel_c[, fdr_label := gsub("e-", "e−", fdr_label, fixed = TRUE)]

p_c <- ggplot(panel_c, aes(event_strength, gene, colour = event_status)) +
  geom_segment(aes(x = 0, xend = event_strength, yend = gene), linewidth = 4.5, alpha = 0.18) +
  geom_point(size = 3.2) +
  geom_text(
    aes(label = fdr_label), hjust = -0.08, size = 2.1,
    colour = pal[["ink"]], fontface = "bold"
  ) +
  scale_colour_manual(values = c(
    "Detected in published C9 snRNA" = pal[["blue"]],
    "Established literature target; not detected in C9 snRNA" = pal[["violet"]]
  )) +
  scale_x_continuous(
    limits = c(0, 1.42), breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Published event-level splice effect (maximum |ΔPSI|)", y = NULL,
    colour = NULL, title = "External TDP-43 event anchors",
    subtitle = "Published values; no local junction or PSI inference",
    tag = "c"
  ) +
  theme(
    axis.text.y = element_text(face = "bold", size = 6.7),
    legend.position = "bottom", legend.direction = "vertical",
    legend.key.height = grid::unit(2.8, "mm"),
    plot.tag.position = c(0.02, 1),
    plot.margin = margin(6, 16, 6, 10, unit = "pt")
  )

# Panel d — capability audit prevents gene counts from being read as splicing inference.
dataset_labels <- c(
  "GSE219280 primary snRNA" = "Primary GSE219280\ngene-count input",
  "Gittings 2023 independent C9 spectrum" = "Published C9-spectrum\nsnRNA study",
  "Ma 2022 TDP-43-negative human neuronal nuclei" = "Published TDP-43-loss\nneuronal nuclei"
)
capability_labels <- c(
  "Junction input", "De novo\ninference", "Reproducible\npublished extraction"
)
panel_d[, dataset_short := factor(dataset_labels[dataset], levels = rev(unname(dataset_labels)))]
panel_d[, capability_short := factor(
  capability_label,
  levels = c("Junction input", "De novo inference", "Reproducible published extraction"),
  labels = capability_labels
)]
panel_d[, status := ifelse(available, "yes", "no")]

p_d <- ggplot(panel_d, aes(capability_short, dataset_short)) +
  geom_tile(aes(fill = available), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = status, colour = available), size = 2.25, fontface = "bold") +
  scale_fill_manual(values = c(`FALSE` = pal[["neutral_soft"]], `TRUE` = pal[["teal"]]), guide = "none") +
  scale_colour_manual(values = c(`FALSE` = pal[["muted"]], `TRUE` = "white"), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL, y = NULL, title = "Input-capability boundary",
    subtitle = "Primary input supports expression—not local junction, isoform or PSI testing",
    tag = "d"
  ) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 6.2, face = "bold", lineheight = 0.95, margin = margin(b = 3)),
    axis.text.y = element_text(size = 6.2, lineheight = 0.95, margin = margin(r = 3)),
    plot.tag.position = c(0.24, 1),
    plot.margin = margin(6, 14, 6, 18, unit = "pt")
  )

# Panel e — descriptive counts copied from the external publication.
panel_e[, disease_group := factor(
  disease_group,
  levels = c("Control", "C9-ALS", "C9-ALS-FTD", "C9-FTD")
)]
panel_e[, gene := factor(gene, levels = c("KALRN", "STMN2"))]

p_e <- ggplot(panel_e, aes(disease_group, junctions, fill = gene)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(
    aes(label = junctions), position = position_dodge(width = 0.72),
    vjust = -0.35, size = 2.15, fontface = "bold", colour = pal[["ink"]]
  ) +
  scale_fill_manual(values = c(KALRN = pal[["blue"]], STMN2 = pal[["amber"]])) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0, 0))) +
  labs(
    x = NULL, y = "Author-reported cryptic-exon junctions", fill = NULL,
    title = "Published C9 event summary",
    subtitle = "Descriptive author-reported totals; not independently recomputed",
    tag = "e"
  ) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 6.2),
    legend.position = "bottom", legend.direction = "horizontal",
    plot.tag.position = c(0.02, 1),
    plot.margin = margin(6, 10, 6, 12, unit = "pt")
  )

middle_row <- p_b | p_c
middle_row <- middle_row + plot_layout(widths = c(1, 1))
bottom_row <- p_d | p_e
bottom_row <- bottom_row + plot_layout(widths = c(1, 1))
figure4 <- p_a / middle_row / bottom_row +
  plot_layout(heights = c(1.30, 0.84, 0.94))

width_mm <- 190
height_mm <- 220
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
prefix <- file.path(draft_root, "Figure4_redesign_draft")

svglite::svglite(paste0(prefix, ".svg"), width = width_in, height = height_in, bg = "white")
print(figure4)
dev.off()

grDevices::cairo_pdf(
  paste0(prefix, ".pdf"), width = width_in, height = height_in,
  family = font_body, bg = "white"
)
print(figure4)
dev.off()

ragg::agg_tiff(
  paste0(prefix, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, background = "white", compression = "lzw"
)
print(figure4)
dev.off()

ragg::agg_png(
  paste0(prefix, "_preview.png"), width = width_in, height = height_in,
  units = "in", res = 300, background = "white"
)
print(figure4)
dev.off()

write_clean_tsv(programme_cards, file.path(draft_source, "Figure4_panel_b_four_programme_groups.tsv"))

qa_lines <- c(
  "Figure 4 redesign draft — statistical and interpretive boundaries",
  "Core conclusion: a finite 17-gene discovery-derived candidate set links programme evidence to external published TDP-43 event-level context within explicit assay limits.",
  "Panel a: 17 candidates by seven evidence layers; displayed numbers are layer-specific supporting units and are not a pooled evidence score or common test statistic.",
  "Panel b: membership is grouped into four manuscript-level programme families. Membership is undirected and does not imply regulation or causation.",
  "Panel c: STMN2, KALRN and UNC13A values are published event-level splice effects. STMN2 and KALRN were detected in published C9 snRNA; UNC13A is literature-established but was not detected in that dataset.",
  "Panel d: GSE219280 local input is gene-count expression data and cannot support local junction, isoform or PSI inference.",
  "Panel e: cryptic-exon junction counts are author-reported descriptive totals and were not independently recomputed.",
  "No gene is labelled as a validated driver, biomarker or causal target."
)
writeLines(qa_lines, file.path(draft_root, "Figure4_redesign_draft_QA.txt"))

cat("Rendered Figure 4 draft to:", draft_root, "\n")
cat("Candidate genes:", uniqueN(panel_a$gene), "\n")
cat("Evidence matrix cells:", nrow(panel_a), "\n")
cat("Programme groups:", nrow(programme_cards), "\n")
cat("Published splice-event genes:", nrow(panel_c), "\n")
cat("Capability audit cells:", nrow(panel_d), "\n")
cat("Published event-count bars:", nrow(panel_e), "\n")
