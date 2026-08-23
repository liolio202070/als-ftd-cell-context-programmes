# Cell-context programmes across the ALS-FTD spectrum

This repository is the frozen code, figure-source-data and compact supplementary-data release for the manuscript **“Cell-context programmes across the ALS-FTD spectrum under external same-nucleus multiome testing.”**

The study reanalyses public human-brain transcriptomic, chromatin-accessibility, paired multiome, perturbation and ALS GWAS resources. The biological donor is retained as the inferential unit. The archive does not redistribute large third-party raw matrices; those must be obtained from the original repositories under their own terms.

## Repository contents

- `figures/`: locked submission figures (TIFF plus editable PDF/SVG where available).
- `data/source_data/`: panel-level Source Data for Figures 1–6 and Extended Data Figures 13–14.
- `data/supplementary/`: compact supplementary tables and statistical dictionaries.
- `data/figure_inputs/`: additional frozen inputs needed by selected figure scripts.
- `scripts/figures/`: portable scripts for Figures 2–6 and Extended Data Figures 13–14.
- `scripts/analysis/`: donor-aware analysis workflows; full execution requires reconstruction of the public raw-data tree.
- `scripts/qa/`: figure-bundle and manuscript numeric-freeze checks.
- `manuscript/`: synchronized manuscript, Word file and bibliography.
- `metadata/`: dataset and release metadata templates.

## Rebuild the quantitative figures

R 4.6.1 and the package versions in `environment/R-packages.tsv` were used for the frozen release. If the packages are installed in a project-specific library, expose it as `ALS_FTD_R_LIBRARY`.

```bash
export ALS_FTD_R_LIBRARY=/path/to/R_library
Rscript scripts/figures/render_all_figures.R
Rscript scripts/qa/qa_figure_bundle.R
Rscript scripts/qa/qa_manuscript_numeric_freeze.R
```

Outputs are written below `outputs/figures/`; locked submission exports remain unchanged in `figures/`.

Figures 2–5, the final fused Figure 6 and Extended Data Figure 13 reproduce the locked TIFF files byte-for-byte in the frozen environment. Extended Data Figure 14 reproduces the same dimensions, displayed values and visual layout; its exported file hash differs. Figure 1 is an author-final overview schematic supplied as a locked TIFF with its dataset-role Source Data, but no exact editable master or deterministic rendering script is asserted. See `docs/REPRODUCIBILITY_LEVELS.md`.

## Re-run the complete analyses

Set `ALS_FTD_PROJECT_ROOT` to a reconstructed project directory containing the public inputs described in `metadata/DATASET_MANIFEST.tsv`. The scripts default to this repository root when the variable is absent.

```bash
export ALS_FTD_PROJECT_ROOT=/path/to/reconstructed_project
export ALS_FTD_R_LIBRARY=/path/to/R_library
Rscript scripts/analysis/build_GSE219280_pseudobulk.R
```

The analysis scripts are the frozen computational record, not a promise that third-party data downloads or high-memory workflows will run without the source repositories, platform-specific tools and compute resources documented by each module. Derived Source Data and compact supplementary results needed to audit the paper are included directly.

## Data and code availability

Ready-to-paste statements are in `DATA_AVAILABILITY.md` and `CODE_AVAILABILITY.md`. DOI, repository URL, creator list and license remain placeholders until the authors approve an external GitHub/Zenodo release. No external upload has been performed from this draft.

## License status

No reuse license has yet been selected. Until the authors add one, normal copyright restrictions apply. See `LICENSE_SELECTION_REQUIRED.md` for the recommended split between code and author-generated data/figures.

