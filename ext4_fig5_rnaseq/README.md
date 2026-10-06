# Human TYMS bulk RNA-seq source files

**Author:** Vignesh Ramesh ([vramesh@bmb.sdu.dk](mailto:vramesh@bmb.sdu.dk))

Source files and analysis commands for Extended Data Fig. 4A-C and Fig. 5A.

## Design and methods

Study-generated A549 bulk RNA-seq: [GSE271721](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271721), 12 libraries across four conditions.

The [analysis command transcript](source/22_08_23_Script_TS.txt) specifies
FastQC 0.11.9, STAR with GRCh38/Ensembl 105, paired-end Rsubread
featureCounts, and DESeq2 design `~ Type`. Twelve libraries represent NTC,
NTC coculture, TS-KO, and TS-KO coculture, with three replicates per group.
The differential-testing unit is the library/sample. Contrasts are TS versus
NTC, TS_CC versus NTC_CC, and TS_CC versus TS.

`metadata/` contains sample and GEO submission records. `tables/` contains
count, DESeq2, and source-data workbooks.
[Prism numeric source data](tables/gap_junctions_prism_source.tsv) preserves
8 genes and 16 decimal values from the original GraphPad table; decimal commas
are converted to points without rounding. The
[provenance manifest](provenance_manifest.json) lists file roles and paths.
Raw FASTQ/BAM files and STAR indices are not included, so complete
read-level processing requires the external sequencing files.

## Validation

Run `python3 scripts/validate_repository.py` from the repository root. The
check covers 14 source records, required workbook sheets and rows, software
and contrast terms, and excluded raw-file types. A missing file or failing table
check requires confirming the input release and contrast definition before
using the result. Capture `python3 -V`, `git rev-parse HEAD`, and R
`sessionInfo()` during any rerun.
