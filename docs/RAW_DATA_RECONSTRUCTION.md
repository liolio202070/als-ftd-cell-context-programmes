# Raw-data reconstruction

Large third-party inputs are intentionally excluded. To execute the complete analysis workflows, download the public files from the repositories listed in `metadata/DATASET_MANIFEST.tsv` and reconstruct the directories expected by the individual scripts below `data/raw/`, `data/metadata/` and `external/`.

Use the release root itself or set:

```bash
export ALS_FTD_PROJECT_ROOT=/path/to/reconstructed_project
export ALS_FTD_R_LIBRARY=/path/to/R_library
```

The major expected locations include:

- `data/raw/` for GSE219279/280, GSE212630, GSE163122, Ruf, Gittings and published splicing inputs.
- `data/metadata/GSE219280_nuclei_metadata.tsv` for discovery-cell annotations.
- `external/magma/` for MAGMA, the GRCh37 reference panel, GENCODE v19 and GCST90027164 summary statistics.
- `external/` for the two TDP-43 perturbation count files.

The full public raw-data tree can exceed 200 GB and some multiome steps require high memory, C++17, ArchR and command-line genomic tools. Figure reproduction does not require these raw files; it uses the frozen panel tables distributed in this repository.

Do not place participant-identifiable or controlled-access data in this release. This study package is designed around public, de-identified resources only.

