# Visium TYMS/TK1 spatial panels

**Author:** Mert Demirdizen ([mert@bmb.sdu.dk](mailto:mert@bmb.sdu.dk))

Analysis code and figure data for Fig. 2B/D and Extended Data Fig. 3F.
The input is the 10x Genomics
`CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma` Visium dataset.
The filtered feature-barcode H5 is not included; its SHA-256 value and size
are in [provenance_manifest.json](provenance_manifest.json).

## Method

[analysis.R](source/analysis.R) reads raw count-layer values for `TYMS` and
`TK1`. A spot is positive for a gene when its count is greater than zero.
The script counts four detection states, uses `dbscan::frNN` on low-resolution
tissue coordinates for radii 1–10, and keeps neighborhoods with the modal
number of neighbors at each radius. It sets `set.seed(123)`.

The R environment requires Seurat, tidyverse, dbscan, magrittr, and
multcompView. The manuscript method reports R 4.4.0. The script's
`Load10X_Spatial` call expects the H5 one directory above its working
directory; record any path change during a rerun.

## Outputs and checks

- [Spot-state counts](results/tyms_tk1_counts.csv): 1,707 TYMS+/TK1+,
  549 TYMS+/TK1−, 1,273 TYMS−/TK1−, and 329 TYMS−/TK1+ spots; total 3,858.
- [Neighborhood data](results/neighborhood_data_points.csv): 160 rows over
  ten radii.
- [Spatial maps](results/extended_3f/) and combined plot exports in
  `results/`.

Run `python3 scripts/validate_repository.py` from the repository root to
check file hashes, state counts, neighborhood rows, and radii. A mismatch
indicates a different input or filtering rule and should be resolved before
using the panel data.
