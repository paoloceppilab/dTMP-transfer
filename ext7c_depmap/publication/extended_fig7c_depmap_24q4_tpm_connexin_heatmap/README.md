# Extended Data Fig. 7C | DepMap 24Q4 connexin expression

**Author:** Mert Demirdizen ([mert@bmb.sdu.dk](mailto:mert@bmb.sdu.dk))

The primary connexin mRNA heatmap for Extended Data Fig. 7C was generated from the
official DepMap 24Q4 Public all-gene RNA expression release.

## Primary outputs

- [Final PNG](figures/extended_fig7c_depmap_24q4_tpm.png), [PDF](figures/extended_fig7c_depmap_24q4_tpm.pdf), and [SVG](figures/extended_fig7c_depmap_24q4_tpm.svg).
- [Title-free PNG](figures/extended_fig7c_depmap_24q4_tpm_no_title.png), with PDF/SVG companions.
- [DepMap TPM matrix](tables/depmap_24q4_tpm_matrix.tsv), [source log2(TPM + 1) matrix](tables/depmap_24q4_log2_tpm_plus1_matrix.tsv), [gene provenance](tables/gene_availability_provenance.tsv), and [profile mapping](tables/sample_profile_mapping.tsv).
- [Methods](methods/depmap_24q4_tpm_methods.md), [caption](methods/extended_fig7c_depmap_24q4_tpm_caption.md), [QC report](qc/qc_report.md), and [run manifest](run_manifest.json).

## Source and units

The [DepMap 24Q4 Public release](https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248)
provides [`OmicsExpressionAllGenesTPMLogp1Profile.csv`](https://ndownloader.figshare.com/files/51065360)
and [`OmicsDefaultModelProfiles.csv`](https://ndownloader.figshare.com/files/51065339).
The source values are `log2(TPM + 1)`; tables and figures report TPM after
`TPM = 2 ** source_value - 1`. Both input SHA-256 hashes and sizes are pinned in
[source_urls_and_checksums.md](docs/source_urls_and_checksums.md).

The final matrix has eight cell lines × 17 genes, no missing values, and an
explicit match for `GJD3 (ENSG00000183153)`. The primary TPM-native color
scale runs from 0 to the observed maximum, 91.71 TPM, without clipping.

## Re-run

From the repository root:

```bash
python3 ext7c_depmap/publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap/scripts/generate_extended_fig7c_depmap_24q4_tpm_connexin_heatmap.py --force
```

The optional `--raw-cache-dir` points to a directory holding validated source
CSVs; otherwise the script downloads them. `--keep-raw` retains local raw copies
under `raw/`. Raw release files are not tracked by Git. The recorded run's
software versions are in `run_manifest.json`; a new run records its own
versions, source validation, QC, and output paths. Matrix and figure generation
has no random component. For verification, run the command from a temporary
copy of the repository as described in the figure-level README; `--force` in the
committed checkout rewrites run-specific records and requires the provenance
hashes to be reviewed and updated.
