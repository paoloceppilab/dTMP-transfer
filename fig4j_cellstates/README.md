# TYMS/TK1 spatial cell-state abundance analysis

**Author:** Mert Demirdizen ([mert@bmb.sdu.dk](mailto:mert@bmb.sdu.dk))

The analysis uses four TYMS/TK1 detection states in tumor epithelial cells,
cell2location abundance estimates in spatial sections, and a seven-factor
nonnegative matrix factorization (NMF) of those estimates.

## Data and code

- scRNA-seq reference: `E-MTAB-13526`.
- Spatial transcriptomics: `E-MTAB-13530`.
- Upstream method: [De Zuani et al., Nature Communications (2024)](https://doi.org/10.1038/s41467-024-48700-8) and [source code at commit 51a8997](https://gitlab.com/cvejic-group/lung/-/tree/51a8997ebffc8a5b529736f45b71f06d16af31ac).
- [Source files](source/) and [analysis provenance](provenance_manifest.json).

The [classification notebook](source/MalignantCellClassification_TYMS_TK1_cycling.ipynb)
reads the external AnnData file from `DTMP_FIG4J_DATA_DIR` (default `data/`)
and writes outputs to `DTMP_FIG4J_OUTPUT_DIR` (default `results/`). It retains
the recorded cell-count output but omits installation and system logs.

The notebook classifies cells from the three specified tumor epithelial
populations into TYMS/TK1 detection states using count values greater than
zero. `estimate_expression_signatures.py` fits the reference model;
`run_cell2location.py` maps cell-state abundances to spatial spots;
`save_all_samples_merged.py` attaches images and coordinates; and
`colocation.py` filters spots and fits NMF. `config.py` defines the reference
labels, exclusions, and tumor sections. The Visium helper is
`visium_qc_and_visualisation.py`.

The source code specifies a minimum 800-count spot filter and evaluates 7–10
NMF factors; the displayed analysis uses seven factors. Five NMF estimator
runs are represented in the model settings. Source code includes
`torch.manual_seed(42)`; the NMF random initialization is not fully controlled
by an explicit seed in the supplied source. Record package versions and random
states for any rerun.

## Interpretation and execution limits

The figure describes spatial co-occurrence of **modeled cell-state
abundances**. It does not measure physical contact, dTMP transfer, or
cell-to-cell exchange directly. A zero TYMS/TK1 transcript count does not by
itself demonstrate loss of enzyme function.

Full execution requires the AnnData reference, section-level Visium inputs,
and model outputs, which are not included here. The exact training environment
is not specified as an installable lockfile. `python3 scripts/validate_repository.py`
checks that the listed source files are present; it does not rerun cell2location or NMF.
