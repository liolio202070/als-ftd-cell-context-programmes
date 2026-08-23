# Reproducibility levels

The release distinguishes exact file reproduction from scientific and visual reproduction.

| Item | Frozen input included | Portable render script | Verified release result |
|---|---:|---:|---|
| Figure 1 | Yes | No | Author-final TIFF and dataset/evidence-role tables supplied; editable provenance remains to be archived |
| Figure 2 | Yes | Yes | TIFF SHA-256 exactly matches locked submission file |
| Figure 3 | Yes | Yes | TIFF SHA-256 exactly matches locked submission file |
| Figure 4 | Yes | Yes | TIFF SHA-256 exactly matches locked submission file |
| Figure 5 | Yes | Yes | TIFF SHA-256 exactly matches locked submission file |
| Figure 6 | Yes | Yes | Final fused TIFF SHA-256 exactly matches locked submission file |
| Extended Data Figure 13 | Yes | Yes | TIFF SHA-256 exactly matches locked submission file |
| Extended Data Figure 14 | Yes | Yes | Same dimensions, values and visual layout; export hash differs |

“Exact” refers to byte-for-byte equality of the 600-dpi TIFF under the frozen local R environment. It does not imply that every operating system, font stack or later package version will emit identical PDF/SVG metadata. Scientific reproducibility is supported by panel-level Source Data, explicit analysis boundaries and numeric-freeze QA.

Figure 1 should not be advertised as code-generated until its editable source (for example SVG, draw.io or presentation master) is deposited and linked to the release.

