# GSE148071 scRNA-seq analysis code

**Author:** Vignesh Ramesh ([vramesh@bmb.sdu.dk](mailto:vramesh@bmb.sdu.dk))

Code for MKI67 and Extended Data Fig. 3A/B/D/E using patient-level
single-cell RNA-seq data from [GSE148071](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148071).
The processed `GSM*_P*_exp.txt` expression matrices are not included.

## Analyses and files

- [Pooled script](scripts/24_05_16_scRNA_Mod_script_Pooled_new.txt): Seurat
  preprocessing, tumor epithelial TYMS/TK1/MKI67 states, UMAP, proportions,
  and heatmap outputs.
- `scripts/24_05_16_scRNA_Mod_script_P*_new.txt`: patient-level cell-state
  classification and marker analyses for 36 patients with identified
  cancer-cell clusters.
- `scripts/MKI67_P*.txt`: MKI67-positive subsets for 21 patients
  with more than 50 MKI67-positive cancer cells.
- [Source-file list](provenance_manifest.json) and the software-version
  listing in `scripts/Readme.txt`.

The software listing includes R 4.3.1, Seurat 5.2.1, SeuratObject 5.0.2,
and Matrix 1.7.0. Record `sessionInfo()` for an independent run. The code
uses the source matrix filenames and expects them in its working directory;
configure paths explicitly before running on another machine.

P19 and P27 were excluded because the processed matrices contained duplicate
gene-name rows. P7, P11, P37, and P42 did not have identified cancer-cell
clusters in the patient-level analysis. These exclusions affect the patient
unit of analysis and should be retained in comparisons.

Run `python3 scripts/validate_repository.py` from the repository root to
check that the 59 listed source files are present. Missing files must be
restored before reviewing the analysis code.
