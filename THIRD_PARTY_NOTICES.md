# License scope and third-party attribution

The root MIT license applies to laboratory-authored code and its code
documentation. It does not relicense raw/public datasets, source workbooks,
experimental figure material or pre-existing third-party code.

| Source | Attribution and retained terms |
| --- | --- |
| `fig4j_cellstates/source/` | Cell2location mapping/helpers derive from De Zuani et al., Nature Communications (2024), https://doi.org/10.1038/s41467-024-48700-8, and https://gitlab.com/cvejic-group/lung at commit `51a8997ebffc8a5b529736f45b71f06d16af31ac`. Original third-party rights remain with their authors; the root MIT license covers laboratory-authored changes only. |
| GSE148071 | Wu et al., https://doi.org/10.1038/s41467-021-22801-0; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148071. Dataset access and reuse follow the source record. |
| Visium input | 10x Genomics CytAssist FFPE human lung squamous-cell carcinoma dataset; source URL, filename and SHA-256 are in `fig2b_d_ext3f_visium/provenance_manifest.json`. Source dataset terms apply. |
| TCGA/cBioPortal | TCGA-LUAD summary through cBioPortal; source and acquisition information are in `ext3c_tcga/provenance_manifest.json`. Source data terms apply. |
| DepMap | DepMap 24Q4 public release; original download URLs and checksums are in `ext7c_depmap/provenance_manifest.json`. Retain the source dataset attribution and reuse terms. |
| Human/mouse RNA-seq and amplicon sources | Laboratory experimental tables and presentation/source files are supplied as scientific source material. The software MIT license does not grant separate rights over the underlying experimental data or third-party resource content. |

Dependencies (including Seurat, DESeq2, Rsubread, cell2location and plotting
libraries) are distributed separately under their respective licenses.
