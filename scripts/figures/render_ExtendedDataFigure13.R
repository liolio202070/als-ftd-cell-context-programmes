#!/usr/bin/env Rscript

# Rebuild Extended Data Figure 13 from the frozen panel-level source data.

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

source_root <- file.path(project_root, "data", "source_data", "extended_figures")
output_root <- file.path(project_root, "outputs", "figures", "ExtendedDataFigure13")
source_output <- file.path(output_root, "source_data")
dir.create(source_output, recursive = TRUE, showWarnings = FALSE)

source_paths <- file.path(
  source_root,
  sprintf("ExtendedDataFigure13_panel_%s_source_data.tsv", letters[1:5])
)
if (!all(file.exists(source_paths))) stop("One or more frozen Extended Data Figure 13 inputs are missing.")
file.copy(source_paths, source_output, overwrite = TRUE)

groups <- fread(source_paths[1])
counts <- fread(source_paths[2])
hm <- fread(source_paths[3])
gd <- fread(source_paths[4], colClasses = list(character = "fc_lab"))
gate <- fread(source_paths[5])

palette_contract <- c(
  neutral_dark = "#2B2B2B", neutral_mid = "#787878", neutral_light = "#D9D9D9",
  control = "#6F6F6F", ALS = "#3F72AF", FTD = "#A85B73",
  discovery = "#3F72AF", replication = "#2E9281", support = "#D18A3B",
  negative = "#3F72AF", zero = "#F2F2F2", positive = "#B55B5B"
)

theme_nature_submission <- function(base_size = 6.5, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "black"),
      axis.ticks = element_line(linewidth = 0.32, colour = "black"),
      axis.ticks.length = grid::unit(1.2, "mm"),
      axis.title = element_text(size = base_size, colour = "black"),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.7),
      legend.key.height = grid::unit(3.0, "mm"),
      legend.key.width = grid::unit(3.3, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size - 0.2, face = "bold"),
      plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.3, colour = palette_contract[["neutral_mid"]]),
      plot.tag = element_text(size = 8, face = "bold"),
      plot.tag.position = c(0, 1),
      panel.grid = element_blank(),
      plot.margin = margin(3, 4, 3, 4, unit = "pt")
    )
}
theme_set(theme_nature_submission())

groups[, group := factor(group, levels = c("Control", "ALS", "sporadic ALS-FTD", "C9 ALS-FTD"))]
pal13 <- c("Control" = "#777777", "ALS" = "#3F72AF", "sporadic ALS-FTD" = "#A85B73", "C9 ALS-FTD" = "#C7899B")
e13a <- ggplot(groups, aes(group, donors, fill = group)) +
  geom_col(width = 0.68) + geom_text(aes(label = share), vjust = -0.32, size = 2.0) +
  scale_fill_manual(values = pal13, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = NULL, y = "Donors", title = "External diagnostic cohort", subtitle = "M1 motor cortex; same-nucleus RNA + ATAC") +
  theme(axis.text.x = element_text(angle = 24, hjust = 1))

contrast_labels <- c(ALS_vs_Control = "ALS vs control", ALSFTD_vs_Control = "ALS-FTD vs control", ALS_vs_ALSFTD = "ALS vs ALS-FTD")
counts[, contrast_label := factor(as.character(contrast_label), levels = unname(contrast_labels))]
counts[, gate := factor(gate, levels = c("Base within-cohort", "Final incl. SVA12", "Three-cohort full gate"))]
e13b <- ggplot(counts, aes(contrast_label, count, fill = gate)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.72) +
  geom_text(aes(label = count), position = position_dodge(width = 0.82), vjust = -0.28, size = 2.0) +
  scale_fill_manual(values = c("#9AC3B8", "#2E9281", "#D18A3B"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Programme axes", title = "Robustness and transportability",
       subtitle = "Zero at three-cohort gate is a retained result") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom",
        legend.key.width = grid::unit(3, "mm"), legend.key.height = grid::unit(2.2, "mm"))

fam_order <- c("Mitochondrial energy", "Protein folding", "Protein stabilization", "Heat response", "Small GTPase", "Synapse/junction", "Axon development")
cell_order <- c("Astrocyte", "Excitatory neuron", "Inhibitory neuron", "Microglia", "OPC", "Oligodendrocyte")
hm[, contrast_label := factor(as.character(contrast_label), levels = unname(contrast_labels))]
hm[, family_label := factor(family_label, levels = rev(fam_order))]
hm[, cell_type_label := factor(cell_type_label, levels = cell_order)]
lim13 <- max(abs(hm$NES), na.rm = TRUE)
e13c <- ggplot(hm, aes(cell_type_label, family_label, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_point(data = hm[final_strict_Ruf_axis == TRUE], shape = 21, fill = "white", colour = "black", size = 1.1, stroke = 0.32) +
  facet_wrap(~contrast_label, nrow = 1) +
  scale_fill_gradient2(low = "#3F72AF", mid = "white", high = "#B55B5B", midpoint = 0,
                       limits = c(-lim13, lim13), oob = squish, name = "NES") +
  labs(x = NULL, y = NULL, title = sprintf("Frozen seven-family test across 126 axes (%d final)", sum(hm$final_strict_Ruf_axis == TRUE)),
       subtitle = "White circles: final FDR + donor-LOO + technical + age + SVA12 gates") +
  theme(axis.text.x = element_text(angle = 48, hjust = 1), legend.position = "bottom",
        legend.key.width = grid::unit(3.5, "mm"), legend.key.height = grid::unit(2.2, "mm"))

gd[, contrast_label := as.character(contrast_label)]
gd[, axis := factor(axis, levels = rev(unique(axis[order(logFC)])))]
e13d <- ggplot(gd, aes(logFC, axis, colour = contrast_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = palette_contract[["neutral_mid"]]) +
  geom_point(size = 2) +
  geom_text(aes(label = fc_lab), hjust = -0.45, size = 1.7, show.legend = FALSE) +
  scale_colour_manual(values = c("ALS-FTD vs control" = "#A85B73", "ALS vs ALS-FTD" = "#3F72AF"), name = NULL) +
  coord_cartesian(xlim = c(min(gd$logFC) - 0.12, max(gd$logFC) + 0.42)) +
  labs(x = "RNA log2 fold change", y = NULL, title = "Fourteen final frozen-gene axes",
       subtitle = "All involve ALS-FTD; none is an ALS-vs-control marker") +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 5.1))

gate[, evidence_gate := factor(evidence_gate, levels = rev(evidence_gate))]
e13e <- ggplot(gate, aes(count, evidence_gate)) +
  geom_col(fill = "#545454", width = 0.62) +
  geom_text(aes(label = count), hjust = -0.25, size = 2.0) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Passing axes", y = NULL, title = "Joint-evidence gate",
       subtitle = "STARD13 promoter: FDR 0.047\ntechnical-adjusted FDR 0.052") +
  theme(axis.text.y = element_text(size = 5.8), plot.subtitle = element_text(size = 5.8, lineheight = 1.02))

figure <- (e13a | e13b) / e13c / (e13d | e13e) +
  plot_layout(heights = c(0.85, 1.45, 1.10)) +
  plot_annotation(tag_levels = "a")

width_in <- 190 / 25.4
height_in <- 178 / 25.4
stem <- file.path(output_root, "ExtendedDataFigure13")
svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in, bg = "white"); print(figure); dev.off()
cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial", bg = "white"); print(figure); dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in, units = "in", res = 600,
               background = "white", compression = "lzw"); print(figure); dev.off()
ragg::agg_png(paste0(stem, "_preview.png"), width = width_in, height = height_in, units = "in", res = 300,
              background = "white"); print(figure); dev.off()

writeLines(c(
  "Ruf 2026: 79-donor external diagnostic same-nucleus cohort; no C9-ALS arm.",
  "40 of 126 frozen programme axes pass the full five-gate chain; the three-cohort ALS gate retains zero axes.",
  "14 frozen genes pass all RNA gates; all involve ALS-FTD comparisons.",
  "ATAC: one primary promoter peak (STARD13, FDR 0.047) not surviving technical adjustment (0.052); zero complete RNA-ATAC chains.",
  "RNA and ATAC derive from the same nuclei and are not independent observations."
), file.path(output_root, "ExtendedDataFigure13_QA.txt"))

message("Rendered Extended Data Figure 13 to: ", output_root)
