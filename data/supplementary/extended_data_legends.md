# Extended Data and Supplementary legends

All supplementary entities are tabular Extended Data outputs generated deterministically from frozen project results. Biological donors, not cells or nuclei, define inferential sample size. Blank fields mean not applicable or not present in the originating frozen output; no statistic was imputed.

## ED01 | pseudobulk qc and donor eligibility

Donor-stratum coverage and 20/30/50-nucleus comparison eligibility.
Source mapping: `results/GSE219280_pseudobulk_DE/sample_inclusion_by_group.tsv`; `results/GSE219280_same_celltype_fixed_programs/comparison_QC.tsv`.

## ED02 | full discovery family matrix and region interactions

Full primary seven-family matrix and shortlist gene region interactions; pathway and gene estimands remain labelled separately.
Source mapping: `results/GSE219280_same_celltype_fixed_programs/fixed_family_fgsea_results.tsv`; `results/GSE219280_pseudobulk_DE/all_gene_contrast_results.tsv`.

## ED03 | discovery threshold and donor loo sensitivity

All family-axis threshold comparisons and donor leave-one-out stability summaries.
Source mapping: `results/GSE219280_same_celltype_fixed_programs/threshold_sensitivity_by_program.tsv`; `results/GSE219280_same_celltype_fixed_programs_LOO/leave_one_donor_out_program_stability.tsv`.

## ED04 | full gittings family results and partition audit

Complete Gittings family output, asset/partition donor mapping, and skipped models.
Source mapping: `results/Gittings_2023_C9_snRNA/fixed_7_family_fgsea_results.tsv`; `results/Gittings_2023_C9_snRNA/CELLxGENE_API_asset_donor_manifest.tsv`; `results/Gittings_2023_C9_snRNA/program_model_skipped_celltypes.tsv`.

## ED05 | stable discordant cross cohort axes

Both-cohort global-FDR significant and >=90% LOO-stable axes with opposite directions.
Source mapping: `results/submission_cross_cohort_evidence_v1/exact_estimand_program_replication.tsv`.

## ED06 | grn ftd cross genetic boundary analysis

Three primary global-FDR GRN-FTD axes with matched LOO and RIN complete-case sensitivity.
Source mapping: `results/GSE163122_FTD_GRN_validation/fixed_family_fgsea_results.tsv`; `results/GSE163122_FTD_GRN_validation/leave_one_donor_out/LOO_program_stability.tsv`; `results/GSE163122_FTD_GRN_validation/RIN_complete_case_sensitivity/main_vs_complete_case_program_comparison.tsv`.

## ED07 | full 17 gene cohort specific effects

All frozen 17-gene cohort-specific effects without pooled meta-analysis.
Source mapping: `results/submission_cross_cohort_evidence_v1/unified_17_gene_evidence_forest_data.tsv`.

## ED08 | tdp43 splicing provenance and capability audit

Published event-level evidence, local assay capability, and author-reported donor tests.
Source mapping: `results/TDP43_splicing_v2/predefined_key_gene_splicing_evidence.tsv`; `results/TDP43_splicing_v2/data_availability_and_claim_audit.tsv`; `results/TDP43_splicing_v2/Gittings_2023_author_reported_donor_tests.tsv`.

## ED09 | gse219279 technical adjustment and motif gate

Technical sensitivity of same-donor ATAC and explicit motif/TF-activity capability boundary.
Source mapping: `results/snATAC_GSE219279/regulatory_v2/technical_adjustment_robustness_summary.tsv`; `results/snATAC_GSE219279/regulatory_v2/RNA_ATAC_targeted_concordance_summary.tsv`; `results/snATAC_GSE219279/regulatory_v2/project_capability_audit.tsv`.

## ED10 | gse212630 multiome qc dar and peak gene gate

Multiome QC, full donor-level DAR output, same-nucleus candidate correlations, and motif gate.
Source mapping: `results/GSE212630_multiome/GSE212630_analysis_QC_summary.tsv`; `results/GSE212630_multiome/GSE212630_ATAC_celltype_donor_DAR_trend.tsv`; `results/GSE212630_multiome/GSE212630_same_nucleus_peak_gene_links_17_targets.tsv`; `results/GSE212630_multiome/GSE212630_motif_interpretation_gate.tsv`.

## ED11 | full donor level composition and loo

All donor-level CLR differential-abundance models and leave-one-donor-out summaries.
Source mapping: `results/pre_submission_composition_abundance/GSE219280_CLR_composition_results.tsv`; `results/pre_submission_composition_abundance/GSE219280_CLR_composition_LOO.tsv`; `results/pre_submission_composition_abundance/Gittings2023_CLR_composition_results.tsv`; `results/pre_submission_composition_abundance/Gittings2023_CLR_composition_LOO.tsv`; `results/pre_submission_composition_abundance/GSE212630_CLR_composition_pTDP_results.tsv`; `results/pre_submission_composition_abundance/GSE212630_CLR_composition_pTDP_LOO.tsv`.

## ED12 | gse153960 bulk composition adjusted sensitivity

Full bulk primary-versus-composition-adjusted gene-effect sensitivity.
Source mapping: `results/GSE153960_composition_sensitivity/primary_vs_composition_adjusted_summary.tsv`; `results/GSE153960_composition_sensitivity/primary_vs_composition_adjusted_gene_effects.tsv`.

## ED13 | Ruf 2026 frozen same-nucleus validation

Frozen external diagnostic validation in 79 motor-cortex donors and 180,016 same-nucleus RNA–ATAC profiles. The table combines 40 final strict programme axes, 14 final strict frozen-gene axes, the complete three-cohort ALS gate, frozen promoter-peak tests, RNA–ATAC joint evidence and the seven-donor C9 ALS–FTD sensitivity. Zero complete three-cohort ALS axes and zero complete RNA–ATAC chains are retained explicitly. Source mapping: `results/Ruf_2026_frozen_validation/final_strict_Ruf_programme_axes.tsv`; `final_strict_Ruf_single_gene_axes.tsv`; `final_three_cohort_ALS_programme_replication.tsv`; `ATAC_targeted_results.tsv`; `RNA_ATAC_joint_results.tsv`; `RNA_C9_ALSFTD_sensitivity_7_family.tsv`.
