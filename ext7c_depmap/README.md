# Extended Fig. 7C: DepMap connexin expression

**Author:** Mert Demirdizen ([mert@bmb.sdu.dk](mailto:mert@bmb.sdu.dk))

The figure reports RNA expression for 17 connexin genes in eight lung-cancer
cell lines using the [DepMap 24Q4 Public release](https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248).
The unit of analysis is a cell-line RNA profile. The heatmap is descriptive;
there is no replicate-level hypothesis test.

## Data and units

- [All-gene RNA expression file](https://ndownloader.figshare.com/files/51065360):
  `OmicsExpressionAllGenesTPMLogp1Profile.csv`.
- [Default model-profile file](https://ndownloader.figshare.com/files/51065339):
  `OmicsDefaultModelProfiles.csv`.
- Release file IDs and checksums for the downloaded inputs are in the
  [input-source record](publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap/docs/source_urls_and_checksums.md).

DepMap values are `log2(TPM + 1)`. The displayed TPM matrix applies
`TPM = 2 ** input_value - 1`. The script requires the all-gene expression
file because a protein-coding-only table can omit connexin genes. `GJD3`
must match `GJD3 (ENSG00000183153)`.

The [TPM matrix](publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap/tables/depmap_24q4_tpm_matrix.tsv)
has eight rows and 17 genes. [PNG](publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap/figures/extended_fig7c_depmap_24q4_tpm.png),
[PDF](publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap/figures/extended_fig7c_depmap_24q4_tpm.pdf),
and [SVG](publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap/figures/extended_fig7c_depmap_24q4_tpm.svg)
exports are supplied. The [package README](publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap/README.md)
links the generation script, additional exports, methods, and QC tables.

## Reproduce and validate

From the repository root, install the dependencies specified in the package
README and run in a disposable copy so generated run dates do not replace
the committed outputs:

```bash
python3 -V
depmap_rerun_root=$(mktemp -d /tmp/extended_7c_depmap_rerun.XXXXXX)
python3 -m pip freeze > "$depmap_rerun_root/python_packages.txt"
git rev-parse HEAD > "$depmap_rerun_root/source_commit.txt"
cp -R ext7c_depmap "$depmap_rerun_root/ext7c_depmap"
cd "$depmap_rerun_root"
python3 ext7c_depmap/publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap/scripts/generate_extended_fig7c_depmap_24q4_tpm_connexin_heatmap.py --force
```

The script downloads both public files or accepts `--raw-cache-dir` with
validated local CSVs. Confirm an 8 × 17 matrix, no missing or negative TPM
values, and matching source SHA-256 values. From the repository root, run
`python3 ext7c_depmap/scripts/validate_figure.py` to
check the committed tables, figure formats, source files, and local links.

The plotted genes are matched to pinned symbol/Ensembl pairs in the generation
script. The displayed values come solely from DepMap 24Q4. The primary figure
and title-free export are the two supplied layout versions.
