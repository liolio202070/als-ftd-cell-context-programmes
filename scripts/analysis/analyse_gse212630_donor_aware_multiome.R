#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ArchR); library(data.table); library(edgeR); library(limma)
  library(Matrix); library(SummarizedExperiment); library(GenomicRanges)
})

set.seed(212630)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[1])
source(file.path(dirname(script_path), "..", "lib", "release_paths.R"))
root <- analysis_project_root()
configure_project_library(root)
raw <- file.path(root, "data/raw/GSE212630_multiome")
out <- file.path(root, "results/GSE212630_multiome")
arrow_dir <- file.path(raw, "ArchR_arrows")
proj_dir <- file.path(out, "ArchR_project")
status_file <- file.path(out, "status.tsv")
dir.create(arrow_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(proj_dir, recursive=TRUE, showWarnings=FALSE)
setwd(out)

status <- function(stage,state,detail,progress) {
  x <- data.table(timestamp=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),
                  stage=stage,state=state,detail=detail,progress=progress)
  tmp <- paste0(status_file,".tmp.",Sys.getpid()); fwrite(x,tmp,sep="\t")
  file.rename(tmp,status_file)
}

manifest <- fread(file.path(out,"GSE212630_paired_donor_manifest.tsv"))
queue <- fread(file.path(out,"GSE212630_full_file_queue.tsv"))
identity <- fread(file.path(out,"GSE212630_full_barcode_identity_audit.tsv"))
stopifnot(nrow(manifest)==26L, identity[identity_qc_pass==TRUE,.N]==26L)

# Frozen before inspecting GSE212630 outcomes: the 17 strict leading-edge genes
# and their seven previously defined mechanism families.
source_shortlist <- file.path(root,"results/leading_edge_evidence_v2/leading_edge_driver_shortlist.tsv")
targets <- fread(source_shortlist)[strict_candidate==TRUE,
  .(gene, original7_family_labels, evidence_tier, spectrum_class)]
stopifnot(nrow(targets)==17L)
target_file <- file.path(out,"GSE212630_predefined_17_genes_7_families.tsv")
fwrite(targets,target_file,sep="\t")
target_genes <- targets$gene
family_long <- targets[, .(family=unlist(strsplit(original7_family_labels,";",fixed=TRUE))), by=gene]
family_long <- family_long[nzchar(family)]
expected_families <- c("Axon development","Heat response","Mitochondrial energy",
 "Protein folding","Protein stabilization","Small GTPase","Synapse / junction assembly")
stopifnot(setequal(unique(family_long$family),expected_families))

addArchRThreads(threads=min(6L,parallel::detectCores()))
addArchRGenome("hg38")

frag_q <- queue[file_type=="atac_fragments"]
h5_q <- queue[file_type=="filtered_feature_bc_matrix_h5"]
frag_q <- frag_q[match(manifest$pair_id,pair_id)]
h5_q <- h5_q[match(manifest$pair_id,pair_id)]
stopifnot(identical(frag_q$pair_id,h5_q$pair_id))

valid_barcodes <- setNames(lapply(h5_q$local_path,function(p) {
  rhdf5::h5read(p,"/matrix/barcodes")
}), h5_q$pair_id)
arrow_files <- file.path(arrow_dir,paste0(frag_q$pair_id,".arrow"))
arrow_ok <- vapply(arrow_files,function(p) {
  file.exists(p) && isTRUE(tryCatch(nrow(rhdf5::h5ls(p))>0,error=function(e) FALSE))
},logical(1))
if (any(file.exists(arrow_files) & !arrow_ok)) {
  bad <- arrow_files[file.exists(arrow_files) & !arrow_ok]
  for (p in bad) file.rename(p,paste0(p,".invalid_",format(Sys.time(),"%Y%m%d%H%M%S")))
}
if (!all(arrow_ok)) {
  status("archr_arrow","running","creating 26 same-nucleus Arrow files; TSS>=4; fragments>=1000",75)
  missing <- which(!arrow_ok)
  old <- setwd(arrow_dir); on.exit(setwd(old),add=TRUE)
  made <- createArrowFiles(
    inputFiles=setNames(frag_q$local_path[missing],frag_q$pair_id[missing]),
    sampleNames=frag_q$pair_id[missing], outputNames=frag_q$pair_id[missing],
    validBarcodes=valid_barcodes[missing], minTSS=4, minFrags=1000,
    addTileMat=TRUE, addGeneScoreMat=TRUE, force=FALSE,
    threads=getArchRThreads()
  )
  setwd(old)
}
stopifnot(length(arrow_files)==26L,all(file.exists(arrow_files)))

saved_project <- file.path(proj_dir,"Save-ArchR-Project.rds")
if (file.exists(saved_project)) {
  proj <- loadArchRProject(proj_dir,showLogo=FALSE)
} else {
  status("archr_project","running","joining RNA H5 to ATAC by shared cell barcode",78)
  proj <- ArchRProject(ArrowFiles=arrow_files,outputDirectory=proj_dir,copyArrows=FALSE)
  seRNA <- import10xFeatureMatrix(input=h5_q$local_path,names=h5_q$pair_id,
                                  strictMatch=FALSE,featureType="Gene Expression")
  proj <- addGeneExpressionMatrix(input=proj,seRNA=seRNA,strictMatch=FALSE,force=TRUE)
  proj <- proj[proj$cellNames[!is.na(proj$Gex_nUMI)],]
}

pair_map <- setNames(manifest$ptdp_group,manifest$pair_id)
order_map <- setNames(manifest$ptdp_order,manifest$pair_id)
published_pairs <- manifest[published_celltype_metadata_available==TRUE,pair_id]
sample_vec <- as.character(proj$Sample)
proj$ptdp_group <- unname(pair_map[sample_vec])
proj$ptdp_order <- unname(order_map[sample_vec])
proj$cohort <- ifelse(sample_vec %in% published_pairs,"Emory","Mayo")

published_meta <- fread(file.path(raw,
 "reference_code_c9alsftd_multiome/meta/Emory_cell_metadata.txt"))
label_map <- setNames(published_meta$major_types,published_meta$perSample_barcodes)
labels <- unname(label_map[proj$cellNames]); labels[is.na(labels)] <- "Unassigned"
proj$major_types <- labels
proj$donor_celltype <- paste(proj$Sample,proj$major_types,sep="__")

cell_qc <- as.data.table(as.data.frame(getCellColData(
  proj,select=c("Sample","ptdp_group","ptdp_order","cohort","major_types",
                "donor_celltype","TSSEnrichment","nFrags","Gex_nUMI"))))
cell_qc[,cell_name:=rownames(as.data.frame(getCellColData(proj,select="Sample")))]
donor_qc <- cell_qc[,.(n_cells=.N,median_TSS=median(TSSEnrichment),
 median_fragments=median(nFrags),median_RNA_UMI=median(Gex_nUMI)),
 by=.(pair_id=Sample,ptdp_group,ptdp_order,cohort)]
fwrite(donor_qc,file.path(out,"GSE212630_donor_multiome_QC.tsv"),sep="\t")
celltype_qc <- cell_qc[major_types!="Unassigned",.(n_cells=.N,
 median_TSS=median(TSSEnrichment),median_fragments=median(nFrags),
 median_RNA_UMI=median(Gex_nUMI)),
 by=.(pair_id=Sample,ptdp_group,ptdp_order,cohort,major_types,donor_celltype)]
fwrite(celltype_qc,file.path(out,"GSE212630_donor_celltype_multiome_QC.tsv"),sep="\t")

primary_cells <- proj$cellNames[proj$major_types %in% setdiff(unique(published_meta$major_types),"Unassigned")]
primary <- proj[primary_cells,]
status("peak_calling","running","Emory annotated same-nucleus cells; donor-aware groups; unified MACS3 peak set",81)
macs <- file.path(root,"cache/macs3_venv/bin/macs3")
stopifnot(file.exists(macs))
if (!"PeakMatrix" %in% getAvailableMatrices(primary)) {
  primary <- addGroupCoverages(primary,groupBy="major_types",sampleLabels="Sample",
    minCells=40,minReplicates=2,maxReplicates=5,force=TRUE)
  primary <- addReproduciblePeakSet(primary,groupBy="major_types",pathToMacs2=macs,
    cutOff=0.01,reproducibility="2",force=TRUE)
  primary <- addPeakMatrix(primary,force=TRUE)
}
# A prior interrupted run can leave a valid PeakMatrix in all Arrow files but
# not the in-memory ArchR project peakSet. Recover the exact unified peak rows
# from the stored matrix rather than rerunning MACS3 or treating peakSet=NULL as
# a biological failure.
if (is.null(getPeakSet(primary))) {
  status("peak_recovery","running",
    "recovering unified peak coordinates from completed Arrow PeakMatrix",84)
  recovery_se <- getMatrixFromArrow(ArrowFile=getArrowFiles(primary)[1],
                                    useMatrix="PeakMatrix")
  recovered_peaks <- rowRanges(recovery_se)
  stopifnot(inherits(recovered_peaks,"GRanges"), length(recovered_peaks)>0L)
  primary <- addPeakSet(primary,peakSet=recovered_peaks,force=TRUE)
}

count_assay <- function(se) {
  nm <- assayNames(se); if ("sum" %in% nm) assay(se,"sum") else assay(se,nm[1])
}
row_ids <- function(se,type) {
  if (type=="RNA") as.character(rowData(se)$name) else {
    rd <- as.data.frame(rowData(se))
    if (all(c("seqnames","start","end") %in% colnames(rd))) {
      paste0(rd$seqnames,":",rd$start,"-",rd$end)
    } else if (!is.null(rownames(rd))) {
      rownames(rd)
    } else stop("Peak group matrix has no genomic row identifiers")
  }
}
make_design <- function(meta,include_cohort=FALSE) {
  m <- copy(meta)
  for (v in c("median_TSS","median_fragments","median_RNA_UMI")) {
    z <- paste0("z_",v); m[[z]] <- as.numeric(scale(log1p(m[[v]])))
  }
  rhs <- c("ptdp_order","z_median_TSS","z_median_fragments","z_median_RNA_UMI")
  if(include_cohort) rhs <- c(rhs,"cohort")
  d <- model.matrix(as.formula(paste("~",paste(rhs,collapse="+"))),m)
  if(qr(d)$rank<ncol(d)) stop("Rank-deficient predeclared technical model")
  list(meta=m,design=d,trend_col=which(colnames(d)=="ptdp_order"),
       formula=paste(rhs,collapse="+"))
}
fit_edgeR <- function(counts,features,meta,include_cohort=FALSE,layer,cell_type) {
  design_obj <- make_design(meta,include_cohort); d <- design_obj$design
  y <- DGEList(counts=counts); keep <- filterByExpr(y,design=d,min.count=5)
  y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE])
  y <- estimateDisp(y,d,robust=TRUE); fit <- glmQLFit(y,d,robust=TRUE)
  q <- glmQLFTest(fit,coef=design_obj$trend_col); z <- topTags(q,n=Inf,sort.by="none")$table
  ans <- data.table(feature=features[keep],log2FC_per_pTDP_level=z$logFC,
    logCPM=z$logCPM,F=z$F,p_value=z$PValue,FDR=z$FDR,
    layer=layer,cell_type=cell_type,n_donors=ncol(counts),
    model=design_obj$formula)
  list(result=ans,logCPM=cpm(y,log=TRUE,prior.count=2),meta=design_obj$meta)
}
loo_direction <- function(logcpm,meta,features,include_cohort=FALSE,layer,cell_type) {
  wanted <- intersect(features,rownames(logcpm)); if(!length(wanted)) return(data.table())
  full <- make_design(meta,include_cohort); X <- full$design; k <- full$trend_col
  beta_full <- apply(logcpm[wanted,,drop=FALSE],1,function(y) coef(lm.fit(X,y))[k])
  signs <- sapply(seq_len(nrow(meta)),function(j) {
    o <- make_design(meta[-j],include_cohort); Xj <- o$design; kj <- o$trend_col
    apply(logcpm[wanted,-j,drop=FALSE],1,function(y) sign(coef(lm.fit(Xj,y))[kj]))
  })
  data.table(feature=wanted,full_beta=beta_full,
    loo_same_direction_fraction=rowMeans(signs==sign(beta_full),na.rm=TRUE),
    loo_min_sign=apply(signs,1,min,na.rm=TRUE),loo_max_sign=apply(signs,1,max,na.rm=TRUE),
    layer=layer,cell_type=cell_type,n_leaveouts=nrow(meta))
}

rna_full_se <- getGroupSE(proj,useMatrix="GeneExpressionMatrix",groupBy="Sample",
                           divideN=FALSE,threads=getArchRThreads())
rna_full_counts <- count_assay(rna_full_se); rna_features <- row_ids(rna_full_se,"RNA")
rownames(rna_full_counts) <- rna_features
colnames(rna_full_counts) <- sub("^Sample_","",colnames(rna_full_counts))
meta_full <- merge(data.table(pair_id=colnames(rna_full_counts)),donor_qc,by="pair_id",sort=FALSE)
rna_full_counts <- rna_full_counts[,meta_full$pair_id,drop=FALSE]
rna_full_fit <- fit_edgeR(rna_full_counts,rna_features,meta_full,TRUE,"RNA","whole_tissue")
fwrite(rna_full_fit$result,file.path(out,"GSE212630_RNA_whole_tissue_donor_trend.tsv"),sep="\t")
rna_loo <- list(loo_direction(rna_full_fit$logCPM,rna_full_fit$meta,target_genes,TRUE,"RNA","whole_tissue"))

rna_ct_se <- getGroupSE(primary,useMatrix="GeneExpressionMatrix",groupBy="donor_celltype",
                         divideN=FALSE,threads=getArchRThreads())
rna_ct_counts <- count_assay(rna_ct_se); rownames(rna_ct_counts) <- row_ids(rna_ct_se,"RNA")
colnames(rna_ct_counts) <- sub("^donor_celltype_","",colnames(rna_ct_counts))
rna_ct_results <- list(); rna_ct_log <- list(); rna_ct_meta <- list()
for(ct in sort(unique(celltype_qc$major_types))) {
  m <- celltype_qc[major_types==ct & n_cells>=50]
  groups <- intersect(paste(m$pair_id,ct,sep="__"),colnames(rna_ct_counts))
  m <- m[match(sub(paste0("__",ct,"$"),"",groups),pair_id)]
  if(nrow(m)<8 || uniqueN(m$ptdp_group)<4) next
  fit <- fit_edgeR(rna_ct_counts[,groups,drop=FALSE],rownames(rna_ct_counts),m,FALSE,"RNA",ct)
  rna_ct_results[[ct]] <- fit$result; rna_ct_log[[ct]] <- fit$logCPM; rna_ct_meta[[ct]] <- fit$meta
  rna_loo[[length(rna_loo)+1L]] <- loo_direction(fit$logCPM,fit$meta,target_genes,FALSE,"RNA",ct)
}
rna_ct_all <- rbindlist(rna_ct_results,fill=TRUE)
fwrite(rna_ct_all,file.path(out,"GSE212630_RNA_celltype_donor_trend.tsv"),sep="\t")

# Apply the Emory-derived unified peak set to the complete same-nucleus cohort.
peak_set_for_all <- getPeakSet(primary)
stopifnot(inherits(peak_set_for_all,"GRanges"),length(peak_set_for_all)>0L)
proj <- addPeakSet(proj,peakSet=peak_set_for_all,force=TRUE)
if (!"PeakMatrix" %in% getAvailableMatrices(proj)) proj <- addPeakMatrix(proj,force=TRUE)
peak_ct_se <- getGroupSE(primary,useMatrix="PeakMatrix",groupBy="donor_celltype",
                          divideN=FALSE,threads=getArchRThreads())
peak_ct_counts <- count_assay(peak_ct_se); peak_features <- row_ids(peak_ct_se,"ATAC")
rownames(peak_ct_counts) <- peak_features
colnames(peak_ct_counts) <- sub("^donor_celltype_","",colnames(peak_ct_counts))
peak_results <- list(); peak_log <- list(); peak_meta <- list()
for(ct in names(rna_ct_meta)) {
  m <- celltype_qc[major_types==ct & n_cells>=50]
  groups <- intersect(paste(m$pair_id,ct,sep="__"),colnames(peak_ct_counts))
  m <- m[match(sub(paste0("__",ct,"$"),"",groups),pair_id)]
  fit <- fit_edgeR(peak_ct_counts[,groups,drop=FALSE],peak_features,m,FALSE,"ATAC_peak",ct)
  peak_results[[ct]] <- fit$result; peak_log[[ct]] <- fit$logCPM; peak_meta[[ct]] <- fit$meta
}
peak_all <- rbindlist(peak_results,fill=TRUE)
fwrite(peak_all,file.path(out,"GSE212630_ATAC_celltype_donor_DAR_trend.tsv"),sep="\t")
status("joint_peak_gene","running",
 "RNA and ATAC donor models complete; mapping frozen 17 genes to unified PeakMatrix rows",91)

# Same-nucleus donor-metacell peak-gene links for the frozen 17 genes only.
genes <- getGeneAnnotation(primary)$genes
sym_col <- intersect(c("symbol","name","gene_name"),colnames(mcols(genes)))[1]
stopifnot(!is.na(sym_col)); gene_symbols <- as.character(mcols(genes)[[sym_col]])
names(genes) <- gene_symbols
peaks <- getPeakSet(primary)
# PeakMatrix row identifiers can differ from freshly formatted GRanges strings
# (for example, coordinate convention or internal naming). The PeakMatrix and
# getPeakSet rows share an invariant order, so use the actual modeled matrix row
# identifiers for all indexing and retain the GRanges only for distances.
stopifnot(length(peaks) == length(peak_features))
peak_ids <- peak_features
links <- list()
for(ct in intersect(names(rna_ct_log),names(peak_log))) {
  common_groups <- intersect(colnames(rna_ct_log[[ct]]),colnames(peak_log[[ct]]))
  rmta <- rna_ct_meta[[ct]][match(sub(paste0("__",ct,"$"),"",common_groups),pair_id)]
  tech <- model.matrix(~scale(log1p(median_TSS))+scale(log1p(median_fragments))+
                         scale(log1p(median_RNA_UMI)),rmta)
  for(g in intersect(target_genes,intersect(rownames(rna_ct_log[[ct]]),names(genes)))) {
    window <- resize(genes[g],width=1,fix="start")
    window <- resize(window,width=500001,fix="center")
    idx <- queryHits(findOverlaps(peaks,window,ignore.strand=TRUE))
    if(!length(idx)) next
    d <- distance(peaks[idx],resize(genes[g],width=1,fix="start"),ignore.strand=TRUE)
    idx <- idx[order(d)][seq_len(min(50L,length(idx)))]
    gy <- as.numeric(rna_ct_log[[ct]][g,common_groups]); gr <- residuals(lm.fit(tech,gy))
    for(pi in idx) {
      if (pi > nrow(peak_log[[ct]])) next
      py <- as.numeric(peak_log[[ct]][pi,common_groups]); pr <- residuals(lm.fit(tech,py))
      tst <- suppressWarnings(cor.test(gr,pr,method="spearman",exact=FALSE))
      links[[length(links)+1L]] <- data.table(cell_type=ct,gene=g,peak=peak_ids[pi],
        distance_to_gene=distance(peaks[pi],resize(genes[g],width=1,fix="start"),
                                  ignore.strand=TRUE),
        rho=unname(tst$estimate),p_value=tst$p.value,
        n_donors=length(common_groups),aggregation="same_nucleus_donor_celltype_metacell",
        adjustment="residualized_TSS_fragments_RNA_UMI")
    }
  }
}
link_dt <- rbindlist(links,fill=TRUE)
if(nrow(link_dt)) link_dt[,FDR:=p.adjust(p_value,"BH"),by=cell_type]
fwrite(link_dt,file.path(out,"GSE212630_same_nucleus_peak_gene_links_17_targets.tsv"),sep="\t")

# LOO for the nearest unified peak to each frozen target gene.
nearest_targets <- list()
for(g in intersect(target_genes,names(genes))) {
  hit <- distanceToNearest(resize(genes[g],1,"start"),peaks,ignore.strand=TRUE)
  if(length(hit)) nearest_targets[[g]] <- peak_ids[subjectHits(hit)[1]]
}
nearest_dt <- data.table(gene=names(nearest_targets),peak=unlist(nearest_targets))
fwrite(nearest_dt,file.path(out,"GSE212630_predefined_nearest_peak_map.tsv"),sep="\t")
atac_loo <- list()
for(ct in names(peak_log)) atac_loo[[ct]] <- loo_direction(
  peak_log[[ct]],peak_meta[[ct]],nearest_dt$peak,FALSE,"ATAC_peak",ct)
fwrite(rbindlist(c(rna_loo,atac_loo),fill=TRUE),
       file.path(out,"GSE212630_target_leave_one_donor_out_stability.tsv"),sep="\t")

# Seven family scores use only the frozen 17-gene membership; no post-hoc genes.
family_results <- list()
for(ct in names(rna_ct_log)) {
  L <- rna_ct_log[[ct]]; m <- rna_ct_meta[[ct]]; D <- make_design(m,FALSE)
  for(fam in expected_families) {
    gs <- intersect(family_long[family==fam,gene],rownames(L)); if(!length(gs)) next
    score <- colMeans(L[gs,,drop=FALSE]); fit <- lm.fit(D$design,score)
    se <- sqrt(diag(chol2inv(qr.R(fit$qr))) * sum(fit$residuals^2)/fit$df.residual)
    k <- D$trend_col; t <- fit$coefficients[k]/se[k]
    family_results[[length(family_results)+1L]] <- data.table(
      cell_type=ct,family=fam,n_genes=length(gs),genes=paste(gs,collapse=";"),
      beta_per_pTDP_level=fit$coefficients[k],t=t,p_value=2*pt(abs(t),fit$df.residual,lower.tail=FALSE),
      n_donors=nrow(m),model=D$formula)
  }
}
family_dt <- rbindlist(family_results); family_dt[,global_FDR:=p.adjust(p_value,"BH")]
fwrite(family_dt,file.path(out,"GSE212630_RNA_7_family_donor_trend.tsv"),sep="\t")

# Motif interpretation is deliberately withheld unless a donor-level chromVAR
# model can be fit with the same QC covariates and stable LOO direction.
motif_qc <- data.table(stage="motif_chromVAR",run=FALSE,
 reason="chromVAR donor-level model not pre-specified/validated in this pipeline; no motif or TF-activity claim",
 required_before_interpretation="donor-level chromVAR trend adjusted for TSS/fragments/RNA UMI plus LOO direction stability")
fwrite(motif_qc,file.path(out,"GSE212630_motif_interpretation_gate.tsv"),sep="\t")

target_rna <- rbindlist(list(rna_full_fit$result,rna_ct_all))[feature %in% target_genes]
target_rna[,target_global_FDR:=p.adjust(p_value,"BH")]
target_atac <- peak_all[feature %in% nearest_dt$peak]
target_atac[,target_global_FDR:=p.adjust(p_value,"BH")]
fwrite(target_rna,file.path(out,"GSE212630_predefined_17_gene_RNA_results.tsv"),sep="\t")
fwrite(target_atac,file.path(out,"GSE212630_predefined_17_gene_nearest_peak_results.tsv"),sep="\t")

qc_summary <- data.table(
  metric=c("donors","pTDP_groups","Emory_annotated_donors","same_nucleus_cells",
           "unified_peaks","RNA_target_genes","mechanism_families","motif_interpreted"),
  value=c(nrow(donor_qc),uniqueN(donor_qc$ptdp_group),length(published_pairs),nrow(cell_qc),
          length(peaks),length(target_genes),length(expected_families),0),
  interpretation=c("full donor cohort","Control/TDPneg/TDPmed/TDPhigh",
    "primary cell-type analysis","RNA and ATAC shared-barcode cells","MACS3 common peak set",
    "frozen before outcome inspection","frozen before outcome inspection","gate not passed; no TF claim")
)
fwrite(qc_summary,file.path(out,"GSE212630_analysis_QC_summary.tsv"),sep="\t")
saveArchRProject(proj,outputDirectory=proj_dir,overwrite=TRUE,load=FALSE,dropCells=FALSE)
status("multiome_analysis","complete",
 "donor-aware RNA pseudobulk+ATAC DAR+same-nucleus metacell peak-gene+technical adjustment+target LOO complete; motif withheld",99)
