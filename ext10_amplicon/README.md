# Extended Data Fig. 10B/D/E/F: Kras and Trp53 amplicon analysis

**Author:** Mert Demirdizen ([mert@bmb.sdu.dk](mailto:mert@bmb.sdu.dk))

Code, summary tables, and figure exports for mouse `Kras`/`Trp53` amplicon
sequencing. The analysis includes 17 specimens: 2 normal lungs, 10 KP tumors,
and 5 KPTT tumors. Reference and primer/guide definitions are in
[references/amplicons.tsv](references/amplicons.tsv); the reference build is
GRCm38/mm10.

## Panel methods

- **Extended Data Fig. 10B:** [editing summary](tables/extended_10b_editing_summary.tsv)
  and [plotting code](source/plot_figure7b_like.py). `kras_hdr_exact_pct`
  measures exact matching to the eight-base HDR barcode in processed
  forward-oriented reads. `trp53_indel_pct` excludes pure right-edge
  truncation artifacts. `KP_811` Trp53 is missing after QC.
- **Extended Data Fig. 10D:** [cut-site analysis](source/analyze_kras_cutsite_spectrum.py)
  classifies `Kras` read pairs in a ±20 bp window around the cut site after
  amplicon base 91. Displayed non-HDR denominators exclude HDR-positive
  pairs. Small summary tables and plots are supplied; the full pair-call
  table is not included.
- **Extended Data Fig. 10E/F:** [codon-window analysis](source/analyze_kras_codon_window.py)
  and [plotting code](source/plot_kras_g12_g13_hotspot_frequency_nonhdr.py)
  summarize `Kras` codon 12/13 mutations after excluding strict HDR-barcode
  positive alleles.

Production software records include CRISPResso 2.3.3, FastQC 0.12.1,
MultiQC 1.33, and Python 3.10.20. Plot run-info files in `metadata/` record
local Python and plotting-library versions. Capture `python -V`,
`python -m pip freeze`, and the Git commit for any rerun.

## Validation and input limits

Run `python3 scripts/validate_repository.py` from the repository root. It
checks that all listed source files are present, the 17-row editing summary, the 15-tumor
cut-site subset, reported denominators, corrected non-HDR percentages, and
category totals. A mismatch indicates changed reads, filtering, or category
classification and must be resolved before replacing a figure.

FASTQ files and the large read-pair classification table are not included.
Their checksums are recorded in the corresponding run-info JSON files, so
read-level reprocessing requires separately obtaining those exact inputs.
The supplied summary tables support downstream figure reproduction.
