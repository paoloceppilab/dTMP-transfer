# dTMP-transfer

Analysis scripts and figure source data accompanying
*Contact-dependent intercellular dTMP transfer sustains tumor cell proliferation*.

## Analysis guide

| Analysis | Included material |
| --- | --- |
| [Patient-level scRNA-seq](fig2a_ext3_scrnaseq/README.md) | GSE148071 cell-state, MKI67, UMAP and marker scripts |
| [Visium TYMS/TK1](fig2b_d_ext3f_visium/README.md) | Spatial detection-state and neighborhood code and panel data |
| [TCGA-LUAD alterations](ext3c_tcga/README.md) | TYMS/TK1 alteration summary, R script and SVG |
| [Human TYMS bulk RNA-seq](ext4_fig5_rnaseq/README.md) | Count and DESeq2 tables, sample metadata and command transcript |
| [DepMap connexins](ext7c_depmap/README.md) | DepMap 24Q4 data references, 8 × 17 matrix, plotting code and exports |
| [Spatial cell-state abundance](fig4j_cellstates/README.md) | Cell classification, cell2location and NMF source scripts |
| [Mouse KP/KPTT bulk RNA-seq](mouse_kp_kptt_rnaseq/README.md) | Sample sheets, count and DESeq2 tables, workflow transcript and FastQC |
| [Kras/Trp53 amplicons](ext10_amplicon/README.md) | Read-processing and plotting scripts, references, summary tables and QC |

Folder and output names retain the archived source naming. Use the analysis
description and each folder's provenance record to identify the relevant code.
Each analysis README specifies its inputs and environment. Raw sequencing
datasets and trained models are obtained separately as described there.

## Validate the archive

From the repository root, using Python 3.9 or later:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r env/requirements-validation.lock.txt
.venv/bin/python scripts/validate_repository.py
.venv/bin/python ext7c_depmap/scripts/validate_figure.py
.venv/bin/python scripts/validate_release.py
```

Checks cover file integrity, workbook/sample dimensions, figure readability,
panel counts, metadata and relative links. They work in a Git checkout and in
an extracted release ZIP. A failed check reports the affected file or invariant.
These are archive/QC checks; they do not perform full read-level reanalysis or
retrain the spatial models.

To regenerate the TCGA summary SVG with base R:

```bash
Rscript ext3c_tcga/source/plot_extended_3c.R --output /tmp/dtmp_tcga_pie.svg
```

Record `python3 -V`, `python3 -m pip freeze`, `Rscript -e 'sessionInfo()'`,
`uname -a` and, in a Git checkout, `git rev-parse HEAD` for a rerun.
Historical analysis versions and seeds are listed in
[env/software_versions.tsv](env/software_versions.tsv), separately from the
archive-validation environment.

## Code availability and citation

The public repository is https://github.com/paoloceppilab/dTMP-transfer.
Publication release: **v1.0.0**. Citation metadata is in
[CITATION.cff](CITATION.cff); Zenodo metadata is in `.zenodo.json`.

## License and attribution

Laboratory-authored code and associated code documentation are available under
the [MIT license](LICENSE). Data, source workbooks and third-party-derived code
retain their own terms; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
