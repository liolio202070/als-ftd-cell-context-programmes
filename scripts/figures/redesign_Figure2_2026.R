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
draft_root <- file.path(project_root, "outputs", "figures", "Figure2")
draft_source_root <- file.path(draft_root, "source_data")
dir.create(draft_source_root, recursive = TRUE, showWarnings = FALSE)

source_paths <- file.path(
  source_root,
  sprintf("Figure2_panel_%s_source_data.tsv", c("a", "b", "c", "d"))
)
if (!all(file.exists(source_paths))) {
  stop("One or more frozen Figure 2 source-data files are missing.")
}
file.copy(source_paths, draft_source_root, overwrite = TRUE)

panel_a <- fread(source_paths[1])
panel_b <- fread(source_paths[2])
panel_c <- fread(source_paths[3])
panel_d <- fread(source_paths[4])

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
  negative = "#3D6FB1",
  positive = "#B75E68"
)

font_body <- "Arial"
font_display <- "Avenir Next"

theme_fig2 <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = font_body) +
    theme(
      text = element_text(colour = pal[["ink"]]),
      axis.line = element_line(linewidth = 0.35, colour = pal[["ink"]]),
      axis.ticks = element_line(linewidth = 0.35, colour = pal[["ink"]]),
      axis.ticks.length = grid::unit(1.2, "mm"),
      axis.title = element_text(size = 7.2, colour = pal[["ink"]]),
      axis.text = element_text(size = 6.6, colour = pal[["ink"]]),
      strip.background = element_blank(),
      strip.text = element_text(size = 7.0, face = "bold"),
      legend.title = element_text(size = 6.7, face = "bold"),
      legend.text = element_text(size = 6.3),
      plot.title = element_text(
        family = font_display, face = "bold", size = 8.5,
        colour = pal[["ink"]], margin = margin(b = 2.5)
      ),
      plot.subtitle = element_text(
        size = 6.4, colour = pal[["muted"]], lineheight = 1.05,
        margin = margin(b = 5)
      ),
      plot.tag = element_text(
        family = font_display, face = "bold", size = 9.5,
        colour = pal[["ink"]]
      ),
      plot.tag.position = c(0, 1),
      plot.margin = margin(4, 5, 4, 5, unit = "pt")
    )
}

theme_set(theme_fig2())

cell_order <- c("Astrocyte", "OPC", "Oligodendrocyte")
cell_short <- c(Astrocyte = "Astro", OPC = "OPC", Oligodendrocyte = "Oligo")
contrast_order <- c("ALS vs control", "FTD vs control", "ALS vs FTD")
contrast_colours <- c(
  "ALS vs control" = pal[["blue"]],
  "FTD vs control" = pal[["mauve"]],
  "ALS vs FTD" = pal[["violet"]]
)
programme_order <- c(
  "Axon development",
  "Synapse / junction assembly",
  "Small GTPase",
  "Protein stabilization",
  "Protein folding",
  "Heat response",
  "Mitochondrial energy"
)

# Panel a: one continuous quantitative matrix rather than three stacked heatmaps.
panel_a[, contrast_label := factor(contrast_label, levels = contrast_order)]
panel_a[, cell := factor(cell, levels = cell_order)]
panel_a[, set_label := factor(set_label, levels = rev(programme_order))]
panel_a[, x := (as.integer(contrast_label) - 1L) * 4L + as.integer(cell)]
panel_a[, y := as.integer(set_label)]

header <- data.table(
  contrast_label = contrast_order,
  xmin = c(0.5, 4.5, 8.5),
  xmax = c(3.5, 7.5, 11.5),
  x = c(2, 6, 10),
  fill = unname(contrast_colours[contrast_order])
)

max_nes <- ceiling(max(abs(panel_a$NES), na.rm = TRUE) * 2) / 2
p_a <- ggplot(panel_a, aes(x, y)) +
  geom_tile(aes(fill = NES), width = 0.94, height = 0.90,
            colour = "white", linewidth = 0.55) +
  geom_point(
    data = panel_a[global_primary_FDR < 0.05],
    shape = 21, size = 2.05, stroke = 0.55,
    fill = "white", colour = pal[["ink"]]
  ) +
  geom_rect(
    data = header,
    aes(xmin = xmin, xmax = xmax, ymin = 7.67, ymax = 7.92),
    inherit.aes = FALSE, fill = header$fill, colour = NA
  ) +
  geom_text(
    data = header,
    aes(x = x, y = 8.20, label = contrast_label),
    inherit.aes = FALSE, family = font_display, fontface = "bold",
    size = 2.75, colour = pal[["ink"]]
  ) +
  scale_x_continuous(
    breaks = c(1:3, 5:7, 9:11),
    labels = rep(unname(cell_short[cell_order]), 3),
    expand = expansion(add = 0.1)
  ) +
  scale_y_continuous(
    breaks = seq_along(programme_order),
    labels = rev(programme_order),
    expand = expansion(add = c(0.05, 0.08))
  ) +
  scale_fill_gradient2(
    low = pal[["negative"]], mid = "#FAFAF7", high = pal[["positive"]],
    midpoint = 0, limits = c(-max_nes, max_nes), oob = squish,
    breaks = c(-2, -1, 0, 1, 2),
    name = "Normalized enrichment score (NES)",
    guide = guide_colourbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(36, "mm"), barheight = grid::unit(2.7, "mm"),
      ticks.colour = pal[["ink"]], frame.colour = pal[["grid"]]
    )
  ) +
  coord_cartesian(xlim = c(0.45, 11.55), ylim = c(0.45, 8.42), clip = "off") +
  labs(
    title = "Programme effects separate by cell type and phenotype contrast",
    subtitle = "Seven prespecified mechanism families × three cell types × three contrasts; white rings mark global primary FDR < 0.05 (35/63 axes).",
    x = NULL, y = NULL, tag = "a"
  ) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 6.3, face = "bold", margin = margin(t = 2)),
    axis.text.y = element_text(size = 6.8, margin = margin(r = 3)),
    legend.position = "bottom",
    legend.margin = margin(t = -2),
    panel.background = element_rect(fill = pal[["soft"]], colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(7, 6, 0, 8, unit = "pt")
  )

# Panel b: descriptive classes across the 21 programme-cell combinations.
class_summary <- panel_b[, .N, by = disease_control_pattern]
class_summary[, class_label := fifelse(
  disease_control_pattern == "shared_direction",
  "Shared direction",
  "Opposite direction"
)]
class_summary[, class_label := factor(
  class_label,
  levels = c("Opposite direction", "Shared direction")
)]
class_summary[, proportion := N / sum(N)]
class_summary[, fill := fifelse(
  disease_control_pattern == "shared_direction", pal[["teal"]], pal[["amber"]]
)]
class_summary[, count_label := sprintf("%d of 21  ·  %d%%", N, round(100 * proportion))]

p_b <- ggplot(class_summary, aes(N, class_label)) +
  geom_col(aes(fill = fill), width = 0.54, colour = NA) +
  geom_text(
    aes(label = count_label), hjust = -0.08,
    family = font_display, fontface = "bold", size = 2.55,
    colour = pal[["ink"]]
  ) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, 24.5), breaks = c(0, 7, 14, 21), expand = c(0, 0)) +
  labs(
    title = "Phenotype-direction classes",
    subtitle = paste0(
      "19 shared, 2 opposite; all direct ALS-versus-FTD contrasts were FTD-higher.\n",
      "Descriptive classification; no separate test."
    ),
    x = "Programme–cell combinations", y = NULL, tag = "b"
  ) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6.1, lineheight = 0.95),
    panel.grid.major.x = element_line(colour = pal[["grid"]], linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(6, 8, 6, 8, unit = "pt")
  )

# Panel c: threshold retention plus donor-omission stability.
threshold_summary <- panel_c[, .(
  lower = mean(direction_retained, na.rm = TRUE),
  middle = mean(direction_retained, na.rm = TRUE),
  upper = mean(direction_retained, na.rm = TRUE),
  label = sprintf("%d/%d", sum(direction_retained, na.rm = TRUE), .N)
), by = .(threshold, minimum_nuclei)]
threshold_summary[, row_label := sprintf("Threshold ≥%d nuclei", minimum_nuclei)]
threshold_summary[, group := "Threshold sensitivity"]

loo_summary <- panel_a[, .(
  lower = min(loo_direction_retained_fraction, na.rm = TRUE),
  middle = median(loo_direction_retained_fraction, na.rm = TRUE),
  upper = max(loo_direction_retained_fraction, na.rm = TRUE),
  label = sprintf("median %d%%", round(100 * median(loo_direction_retained_fraction, na.rm = TRUE)))
), by = cell]
loo_summary[, row_label := paste0("Donor LOO · ", unname(cell_short[as.character(cell)]))]
loo_summary[, group := "Donor leave-one-out"]

robustness <- rbindlist(list(
  threshold_summary[, .(row_label, group, lower, middle, upper, label)],
  loo_summary[, .(row_label, group, lower, middle, upper, label)]
), use.names = TRUE)
robustness[, row_label := factor(
  row_label,
  levels = rev(c(
    "Threshold ≥20 nuclei", "Threshold ≥50 nuclei",
    "Donor LOO · Astro", "Donor LOO · OPC", "Donor LOO · Oligo"
  ))
)]
robustness[, colour := fifelse(group == "Threshold sensitivity", pal[["blue"]], pal[["teal"]])]

p_c <- ggplot(robustness, aes(middle, row_label)) +
  geom_vline(xintercept = 0.9, linetype = "22", linewidth = 0.35, colour = pal[["muted"]]) +
  geom_segment(aes(x = lower, xend = upper, yend = row_label, colour = colour), linewidth = 1.0) +
  geom_point(aes(colour = colour), size = 2.3) +
  geom_text(aes(x = 1.015, label = label), hjust = 0, size = 2.12, colour = pal[["ink"]]) +
  scale_colour_identity() +
  scale_x_continuous(
    limits = c(0.65, 1.13), breaks = c(0.7, 0.8, 0.9, 1.0),
    labels = percent_format(accuracy = 1), expand = c(0, 0)
  ) +
  labs(
    title = "Direction robustness",
    subtitle = paste0(
      "Threshold sensitivity and donor leave-one-out.\n",
      "Points: mean or median; bars: range across axes."
    ),
    x = "Direction retained", y = NULL, tag = "c"
  ) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6.0),
    panel.grid.major.x = element_line(colour = pal[["grid"]], linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(6, 7, 6, 8, unit = "pt")
  )

# Panel d: focus the forest plot on the eight primary FDR-significant CLR effects.
cell_pretty <- c(
  Endo = "Endothelial",
  OPC = "OPC",
  Inh_PVALB = "PVALB inhibitory",
  Exc_intermediate = "Intermediate excitatory",
  Exc_upper = "Upper-layer excitatory"
)
composition_sig <- panel_d[
  analysis_role == "primary" & !is.na(global_primary_FDR) & global_primary_FDR < 0.05
]
composition_sig[, cell_pretty := unname(cell_pretty[cell_type])]
composition_sig[is.na(cell_pretty), cell_pretty := cell]
composition_sig[, contrast_label := factor(
  contrast_label, levels = c("FTD vs control", "ALS vs FTD")
)]
composition_sig <- composition_sig[order(contrast_label, clr_effect)]

p_d <- ggplot(
  composition_sig,
  aes(clr_effect, reorder(cell_pretty, clr_effect), colour = contrast_label)
) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.40, colour = pal[["muted"]]) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    orientation = "y", width = 0.16, linewidth = 0.62
  ) +
  geom_point(size = 2.45) +
  scale_colour_manual(values = contrast_colours[c("FTD vs control", "ALS vs FTD")], guide = "none") +
  facet_wrap(~contrast_label, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = -2:2, limits = c(-2.35, 2.8), expand = c(0, 0)) +
  labs(
    title = "Captured-nucleus composition context",
    subtitle = "Eight of 42 primary CLR contrasts pass global FDR < 0.05; model-based 95% CI.",
    x = "Donor-level CLR effect", y = NULL, tag = "d"
  ) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6.2),
    strip.text = element_text(
      family = font_display, face = "bold", size = 7.2,
      colour = pal[["ink"]], margin = margin(b = 2)
    ),
    panel.grid.major.x = element_line(colour = pal[["grid"]], linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(6, 6, 6, 8, unit = "pt")
  )

figure2 <- p_a / (p_b | p_c) / p_d +
  plot_layout(heights = c(1.48, 0.72, 0.88))

width_mm <- 183
height_mm <- 165
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
prefix <- file.path(draft_root, "Figure2_redesign_draft")

svglite::svglite(paste0(prefix, ".svg"), width = width_in, height = height_in, bg = "white")
print(figure2)
dev.off()

grDevices::cairo_pdf(
  paste0(prefix, ".pdf"), width = width_in, height = height_in,
  family = font_body, bg = "white"
)
print(figure2)
dev.off()

ragg::agg_tiff(
  paste0(prefix, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, background = "white", compression = "lzw"
)
print(figure2)
dev.off()

ragg::agg_png(
  paste0(prefix, "_preview.png"), width = width_in, height = height_in,
  units = "in", res = 300, background = "white"
)
print(figure2)
dev.off()

qa_lines <- c(
  "Figure 2 redesign draft — statistical and interpretive boundaries",
  "Core conclusion: programme directions vary by cell type and phenotype contrast rather than forming one uniform abnormality.",
  "Panel a: 63 frozen axes; 35 pass BH global primary FDR < 0.05. Biological n = 16 astrocyte donors and 17 OPC/oligodendrocyte donors.",
  "Panel b: descriptive classes only; 19/21 shared disease-control direction and 2/21 opposite direction, all with FTD-higher direct phenotype direction.",
  "Panel c: 63/63 directions retained at both 20- and 50-nucleus sensitivity thresholds. LOO is internal robustness, not external replication.",
  "Panel d: only the eight of 42 primary CLR composition contrasts passing global FDR < 0.05 are displayed; all 42 primary and 28 secondary interaction rows remain in Source Data.",
  "CLR effects describe relative captured-nucleus composition, not absolute histological abundance.",
  "No pooled effect is calculated across NES and CLR scales."
)
writeLines(qa_lines, file.path(draft_root, "Figure2_redesign_draft_QA.txt"))

cat("Rendered Figure 2 draft to:", draft_root, "\n")
cat("Primary programme axes:", nrow(panel_a), "\n")
cat("Programme axes passing global FDR < 0.05:", sum(panel_a$global_primary_FDR < 0.05), "\n")
cat("Primary CLR contrasts shown:", nrow(composition_sig), "of", panel_d[analysis_role == "primary", .N], "\n")
