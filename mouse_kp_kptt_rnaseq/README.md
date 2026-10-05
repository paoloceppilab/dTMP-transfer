# Mouse KP/KPTT bulk RNA-seq

**Source analysis:** Vignesh Ramesh, University of Southern Denmark.

The archived workflow compares KP and KPTT mouse lung tumors using FastQC,
Trimmomatic 0.39, STAR, Ensembl GRCm39 release 108, Rsubread featureCounts,
and DESeq2 design `~ Type`, contrast `KPTT` versus `KP`.

## Inputs and artifacts

- `source/24_01_16_mice_KPTT_Script.txt`: historical workflow commands.
- `metadata/deseq2_sample_info/`: all-sample and three-per-group sample sheets.
- `tables/`: featureCounts, DESeq2 input/results and the source workbook.
- `metadata/fastqc_html/`: per-read FastQC reports.
- `metadata/raw_fastq_md5/`: recorded raw FASTQ MD5 checksums.
- `provenance_manifest.json`: file roles, checksums and QC expectations.

The all-sample analysis contains 3 KP and 5 KPTT libraries (8 total). The
three-per-group analysis contains 3 KP and 3 KPTT libraries (6 total). The unit
of differential-expression analysis is the tumor/library. Raw FASTQ, BAM files
and STAR indices are supplied separately for read-level reanalysis. Configure
input paths in a copy of the historical transcript for a new machine.

Historical Methods versions: STAR 2.7.9a, R 4.2.1, Rsubread 2.10.5, and
DESeq2 1.36.0. No random seed is required for the archived DESeq2 contrast.
Capture `STAR --version`, `Rscript -e 'sessionInfo()'` and the reference FASTA/GTF
checksums during a rerun. Use CPU batch jobs for alignment and differential
expression; set thread counts and memory to the input size.

## QC

Run `python3 scripts/validate_release.py` from the repository root. It verifies
file hashes, the expected workflow commands, eight/six distinct sample IDs,
the sample-sheet groups and the source workbook dimensions: 56,981 rows for
each count/result sheet and nine rows for `SI`, including the header.
Any mismatch requires checking the source file and analysis branch before use.
