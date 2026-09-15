# dTMP-transfer

Analysis code, figure data, and input references for
*Contact-dependent intercellular dTMP transfer sustains tumor cell proliferation*.
Each figure folder describes its data source, analysis method, required inputs,
and outputs.

## Figure guide

| Panel | Status | Files and inputs |
| --- | --- | --- |
| [Fig. 2A; Extended Fig. 3A/B/D/E](fig2a_ext3_scrnaseq/README.md) | Provenance record; external input needed | GSE148071 scRNA-seq scripts. Processed expression matrices are not included. |
| [Fig. 2B/D; Extended Fig. 3F](fig2b_d_ext3f_visium/README.md) | Provenance record; external input needed | Visium analysis script and panel data. The 10x feature-barcode H5 is not included. |
| [Extended Fig. 3C](ext3c_tcga/README.md) | Re-runnable | TCGA-LUAD summary counts, R plotting script, and SVG. |
| [Extended Fig. 4A-C; Fig. 5A](ext4_fig5_rnaseq/README.md) | Source files and provenance; external input needed | Bulk RNA-seq commands, count and DESeq2 tables. FASTQ files are not included. |
| [Extended Fig. 7C](ext7c_depmap/README.md) | Re-runnable | DepMap 24Q4 source IDs/checksums, generation script, TPM matrix, and figure exports. |
| [Fig. 4J](fig4j_cellstates/README.md) | Provenance record; external input needed | TYMS/TK1 cell-state classification, cell2location, and NMF code. Input data and trained models are not included. |
| [Extended Fig. 10B/D/E/F](ext10_amplicon/README.md) | Re-runnable figures; external read-level input needed | Summary tables, plotting scripts, reference metadata, and read-level input checksums. |

“Re-runnable” means the listed output can be generated from included data or
specified public downloads. Analyses that need external inputs list those
requirements in their figure README.

## Validation

Run from the repository root with Python 3:

```bash
python3 --version
python3 scripts/validate_repository.py
python3 ext7c_depmap/scripts/validate_figure.py
```

The checks verify listed source-file presence, table dimensions, figure counts,
and selected panel metrics. A failure identifies a missing file or inconsistent
result to review in the corresponding figure folder before using it.
In a Git checkout, record the code revision with `git rev-parse HEAD`.

## Code availability

Repository URL: [https://github.com/paoloceppilab/dTMP-transfer](https://github.com/paoloceppilab/dTMP-transfer).
No repository-wide license is asserted; follow the attribution and usage terms
of the underlying data and third-party code.
