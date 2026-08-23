#!/usr/bin/env Rscript
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
project_root <- analysis_project_root()
configure_project_library(project_root)
library(data.table)
library(ggplot2)
library(patchwork)
library(svglite)
library(ragg)

ruf <- file.path(project_root, "results", "Ruf_2026_frozen_validation")
pkg <- file.path(project_root, "results", "submission_readiness_v1", "final_submission_package")
sup <- file.path(pkg, "supplement")
fig <- file.path(pkg, "figures")
prev <- file.path(pkg, "previews")
src <- file.path(pkg, "source_data")
dir.create(sup, recursive=TRUE, showWarnings=FALSE)
dir.create(fig, recursive=TRUE, showWarnings=FALSE)
dir.create(prev, recursive=TRUE, showWarnings=FALSE)
dir.create(src, recursive=TRUE, showWarnings=FALSE)

prog_all <- fread(file.path(ruf, "final_strict_Ruf_programme_axes.tsv"))
gene_all <- fread(file.path(ruf, "final_strict_Ruf_single_gene_axes.tsv"))
tri <- fread(file.path(ruf, "final_three_cohort_ALS_programme_replication.tsv"))
atac <- fread(file.path(ruf, "ATAC_targeted_results.tsv"))
joint <- fread(file.path(ruf, "RNA_ATAC_joint_results.tsv"))
c9 <- fread(file.path(ruf, "RNA_C9_ALSFTD_sensitivity_7_family.tsv"))

prog <- prog_all[final_strict_Ruf_axis == TRUE]
genes <- gene_all[final_strict_single_gene == TRUE]

ed <- rbindlist(list(
  cbind(data.table(section="programme_final_strict"), prog),
  cbind(data.table(section="single_gene_final_strict"), genes),
  cbind(data.table(section="three_cohort_ALS_gate"), tri),
  cbind(data.table(section="ATAC_frozen_promoter_audit"), atac),
  cbind(data.table(section="RNA_ATAC_joint_gate"), joint),
  cbind(data.table(section="C9_ALSFTD_sensitivity"), c9)
), fill=TRUE, use.names=TRUE)
fwrite(ed, file.path(sup, "ExtendedData13_Ruf_2026_frozen_same_nucleus_validation.tsv"), sep="\t", quote=FALSE, na="")

groups <- data.table(group=c("Control","ALS","sporadic ALS-FTD","C9 ALS-FTD"), donors=c(32,30,10,7))
groups[, group := factor(group, levels=c("Control","ALS","sporadic ALS-FTD","C9 ALS-FTD"))]
fwrite(groups, file.path(src, "ExtendedDataFigure13_panel_a_source_data.tsv"), sep="\t")

fam_order <- c("Mitochondrial energy","Protein folding","Protein stabilization","Heat response","Small GTPase","Synapse/junction","Axon development")
cell_order <- c("Astrocyte","Excitatory neuron","Inhibitory neuron","Microglia","OPC","Oligodendrocyte")
contrast_labels <- c(ALS_vs_Control="ALS vs control", ALSFTD_vs_Control="ALS-FTD vs control", ALS_vs_ALSFTD="ALS vs ALS-FTD")
hm <- copy(prog_all)
hm[, contrast_label := factor(contrast_labels[contrast], levels=unname(contrast_labels))]
hm[, family_label := factor(family_label, levels=rev(fam_order))]
hm[, cell_type_label := factor(cell_type_label, levels=cell_order)]
fwrite(hm, file.path(src, "ExtendedDataFigure13_panel_b_source_data.tsv"), sep="\t", quote=FALSE, na="")

counts <- rbind(
  prog_all[, .(count=sum(strict_external_programme_replication == TRUE)), by=contrast][, gate := "Base within-cohort"],
  prog_all[, .(count=sum(final_strict_Ruf_axis == TRUE)), by=contrast][, gate := "Final incl. SVA12"],
  data.table(contrast="ALS_vs_Control", count=sum(tri$final_strict_three_cohort_same_direction == TRUE), gate="Three-cohort full gate")
)
counts[, contrast_label := contrast_labels[contrast]]
fwrite(counts, file.path(src, "ExtendedDataFigure13_panel_c_source_data.tsv"), sep="\t", quote=FALSE, na="")

gd <- copy(genes)
gd[, contrast_label := contrast_labels[contrast]]
gd[, axis := paste(cell_type_label, gene, sep=" · ")]
fwrite(gd, file.path(src, "ExtendedDataFigure13_panel_d_source_data.tsv"), sep="\t", quote=FALSE, na="")

gate <- data.table(
  evidence_gate=c("Final RNA gene axes","Primary peak-wide ATAC","Technical-adjusted peak-wide ATAC","Complete RNA-ATAC chain"),
  count=c(nrow(genes), sum(atac$FDR < 0.05, na.rm=TRUE), sum(atac$technical_peak_wide_FDR < 0.05, na.rm=TRUE), sum(joint$complete_joint_chain == TRUE, na.rm=TRUE))
)
fwrite(gate, file.path(src, "ExtendedDataFigure13_panel_e_source_data.tsv"), sep="\t")

pal <- c("Control"="#777777","ALS"="#3F72AF","sporadic ALS-FTD"="#A85B73","C9 ALS-FTD"="#C7899B")
theme_pub <- theme_classic(base_family="Arial", base_size=8) + theme(plot.tag=element_text(face="bold", size=11), strip.background=element_blank(), strip.text=element_text(face="bold"), legend.position="bottom")

p1 <- ggplot(groups, aes(group, donors, fill=group)) + geom_col(width=.7) + geom_text(aes(label=donors), vjust=-.3, size=2.7) + scale_fill_manual(values=pal, guide="none") + scale_y_continuous(expand=expansion(mult=c(0,.12))) + labs(x=NULL,y="Donors",title="External diagnostic cohort") + theme_pub + theme(axis.text.x=element_text(angle=25,hjust=1))

p2 <- ggplot(hm, aes(cell_type_label, family_label, fill=NES)) + geom_tile(colour="white", linewidth=.2) + geom_point(data=hm[final_strict_Ruf_axis==TRUE], shape=21, fill="white", colour="black", size=1.2, stroke=.35) + facet_wrap(~contrast_label, nrow=1) + scale_fill_gradient2(low="#3F72AF",mid="white",high="#B55B5B",midpoint=0,name="NES") + labs(x=NULL,y=NULL,title="Frozen seven-family test",subtitle="White circle: final FDR + LOO + technical + age + SVA12 gate") + theme_pub + theme(axis.text.x=element_text(angle=50,hjust=1),legend.position="right")

counts[, gate := factor(gate, levels=c("Base within-cohort","Final incl. SVA12","Three-cohort full gate"))]
p3 <- ggplot(counts, aes(contrast_label, count, fill=gate)) + geom_col(position=position_dodge(width=.8),width=.72) + geom_text(aes(label=count),position=position_dodge(width=.8),vjust=-.25,size=2.5) + scale_fill_manual(values=c("#9AC3B8","#2E9281","#D18A3B")) + scale_y_continuous(expand=expansion(mult=c(0,.15))) + labs(x=NULL,y="Programme axes",fill=NULL,title="Robustness and transportability") + theme_pub + theme(axis.text.x=element_text(angle=22,hjust=1))

gd[, axis := factor(axis, levels=rev(unique(axis[order(logFC)])))]
p4 <- ggplot(gd, aes(logFC, axis, colour=contrast_label)) + geom_vline(xintercept=0,colour="#BDBDBD",linewidth=.3) + geom_point(size=2) + scale_colour_manual(values=c("ALS-FTD vs control"="#A85B73","ALS vs ALS-FTD"="#3F72AF")) + labs(x="RNA log2 fold change",y=NULL,colour=NULL,title="Fourteen final frozen-gene axes") + theme_pub + theme(legend.position="bottom")

gate[, evidence_gate := factor(evidence_gate, levels=rev(evidence_gate))]
p5 <- ggplot(gate, aes(count,evidence_gate)) + geom_col(fill="#545454",width=.65) + geom_text(aes(label=count),hjust=-.25,size=2.7) + scale_x_continuous(expand=expansion(mult=c(0,.18))) + labs(x="Passing axes",y=NULL,title="Joint-evidence gate",subtitle="No complete chromatin-to-expression chain") + theme_pub

layout <- (p1 | p3) / p2 / (p4 | p5) + plot_annotation(tag_levels="a") + plot_layout(heights=c(.75,1.5,1.25))
base <- file.path(fig,"ExtendedDataFigure13")
svglite(paste0(base,".svg"),width=7.2,height=9.0); print(layout); dev.off()
cairo_pdf(paste0(base,".pdf"),width=7.2,height=9.0,family="Arial"); print(layout); dev.off()
agg_tiff(paste0(base,".tiff"),width=7.2,height=9.0,units="in",res=600,compression="lzw",background="white"); print(layout); dev.off()
agg_png(file.path(prev,"ExtendedDataFigure13.png"),width=7.2,height=9.0,units="in",res=220,background="white"); print(layout); dev.off()

cat("programme_rows\t",nrow(prog),"\n",sep="")
cat("gene_rows\t",nrow(genes),"\n",sep="")
cat("three_cohort_full\t",sum(tri$final_strict_three_cohort_same_direction == TRUE),"\n",sep="")
cat("complete_joint\t",sum(joint$complete_joint_chain == TRUE),"\n",sep="")

qa <- data.table(
  check=c("programme_final_40","gene_final_14","three_cohort_zero","primary_ATAC_one","technical_ATAC_zero","complete_joint_zero","ED13_rows_618","ED13_unique_section_present"),
  observed=c(nrow(prog),nrow(genes),sum(tri$final_strict_three_cohort_same_direction == TRUE),sum(atac$FDR < 0.05,na.rm=TRUE),sum(atac$technical_peak_wide_FDR < 0.05,na.rm=TRUE),sum(joint$complete_joint_chain == TRUE),nrow(ed),uniqueN(ed$section)),
  expected=c(40,14,0,1,0,0,618,6)
)
qa[, pass := observed == expected]
fwrite(qa, file.path(pkg,"qa","Ruf_2026_submission_integration_QA.tsv"), sep="\t")
if (!all(qa$pass)) stop("Ruf submission integration QA failed")
