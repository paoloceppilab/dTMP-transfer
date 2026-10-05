# dTMP-transfer

Analysis scripts and figure source data accompanying
*Contact-dependent intercellular dTMP transfer sustains tumor cell proliferation*.

## Analysis guide

The archive follows the five supplied bioinformatics Methods sections.
The source document checksum and section-to-directory mapping are recorded in
[analysis_scope.json](analysis_scope.json).

| Methods section | Included material |
| --- | --- |
| [RNA-sequencing](ext4_fig5_rnaseq/README.md) | Human A549 count and DESeq2 tables, sample metadata and command transcript |
| [Single-cell RNA-sequencing](fig2a_ext3_scrnaseq/README.md) | GSE148071 cell-state, MKI67, UMAP and marker scripts |
| [Visium neighborhood analysis](fig2b_d_ext3f_visium/README.md) | Spatial detection-state and neighborhood code and panel data |
| [Spatial transcriptomics deconvolution](fig4j_cellstates/README.md) | Reference classification, cell2location and NMF source scripts |
| [Amplicon sequencing](ext10_amplicon/README.md) | Kras/Trp53 read-processing and plotting scripts, references, tables and QC |

Folder and output names retain the archived source naming. Use the Methods
section and each folder's provenance record to identify the relevant code.
Each analysis README specifies inputs and execution limits. Raw sequencing
datasets and trained models are obtained separately as described there.

## Validate the archive

From the repository root, using Python 3.9 or later:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r env/requirements-validation.lock.txt
.venv/bin/python scripts/validate_repository.py
.venv/bin/python scripts/validate_release.py
```

Checks cover all file hashes, workbook dimensions, Visium category counts and
radius coverage, amplicon denominators, Methods scope, citation metadata and
relative links. Expected values are recorded in the validators and provenance
manifests. A failed check identifies the affected file or invariant; restore
the matching source or confirm its provenance before using it.

These checks run in a Git checkout and an extracted ZIP. They validate archive
integrity and retained scientific outputs; they do not repeat read alignment,
differential expression, or spatial-model training.

Record `python3 -V`, `python3 -m pip freeze`, `Rscript -e 'sessionInfo()'`,
`uname -a` and, in a Git checkout, `git rev-parse HEAD` during a rerun.
Historical analysis versions and seeds are in
[env/software_versions.tsv](env/software_versions.tsv), separate from the
archive-validation lockfile.

## Code availability and citation

Public repository: https://github.com/paoloceppilab/dTMP-transfer.
Publication package: **v1.0.1**. Citation metadata is in
[CITATION.cff](CITATION.cff); Zenodo metadata is in `.zenodo.json`.

## License and attribution

Laboratory-authored code and its documentation are available under the
[MIT license](LICENSE). Data, workbooks and third-party-derived code retain
their own terms; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
