# Statistical reproducibility index

This archive indexes the exact session-information files that already existed when the analysis modules were completed. Missing historical environments are marked `unavailable_historical`; current installed versions are never back-attributed as an exact missing session.

## GSE212630 pTDP design

The primary annotated Emory cell-type models used 14 biological donors: Control 5, TDPneg 3, TDPmed 3 and TDPhigh 3. Ordinal pathology was coded 0, 1, 2 and 3, respectively, and the reported coefficient is a one-step linear trend. The small, unbalanced groups limit formal assessment of non-linearity. The complete public cohort contains 26 donors (7/6/5/8 across the four groups); the 12 Mayo donors without the frozen published major-cell-type annotation were omitted from the primary cell-type models but included in the cohort-adjusted whole-tissue model.

## Multiplicity

The exact BH scopes are frozen in `supplement/multiplicity_and_FDR_dictionary.tsv`. Candidate RNA `target_global_FDR` spans 119 tests, including 17 whole-tissue and 102 cell-type rows. Peak–gene correlation FDR is BH separately within each cell type. Genome-wide DAR FDR is returned within each cell-type edgeR model. These fields are not interchangeable with transcriptome-wide RNA FDR or programme global-primary FDR.

## Source Data

Individual Source Data files do not all repeat statistical unit, estimand, donor n and correction-family columns. `source_data/panel_statistical_dictionary.tsv` supplies the panel-keyed definitions and identifies when donor n is explicit versus legend/contract-defined. `source_data/source_data_field_dictionary.tsv` defines recurring fields. This dictionary supersedes any broader statement that every TSV is self-contained.

## Software versions

The retained runtime verifies R 4.6.1, ArchR 1.0.3, edgeR 4.10.1, limma 3.68.4, data.table 1.18.4, fgsea 1.38.0, ggplot2 4.0.3 and MACS3 3.0.4. Exact historical sessionInfo is unavailable for GSE212630 multiome, GSE219280 upstream pseudobulk DE and GSE153960 bulk/deconvolution modules. The retained GSE212630 script imports ArchR/edgeR/limma/data.table and points to the retained MACS3 executable; this supports reconstruction but is not represented as a contemporaneously captured full session.
