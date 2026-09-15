# Extended Fig. 3C: TCGA-LUAD TYMS/TK1 alterations

**Author:** Mert Demirdizen ([mert@bmb.sdu.dk](mailto:mert@bmb.sdu.dk))

This panel summarizes TYMS/TK1 amplification and missense-mutation categories
for 496 TCGA-LUAD cases with the required genomic data. The summary source is
NCI GDC TCGA via cBioPortal, accessed 2025-05-16. The unit of analysis is a
case; the six plotted categories are mutually exclusive.

## Input, method, and output

- [Category counts](data/alteration_summary.csv): labels and integer counts.
- [Plotting script](source/plot_extended_3c.R): base R, no external R packages.
- [Figure SVG](results/extended_3c_pie.svg).

From the repository root:

```bash
Rscript ext3c_tcga/source/plot_extended_3c.R \
  --data ext3c_tcga/data/alteration_summary.csv \
  --output ext3c_tcga/results/extended_3c_pie.svg
Rscript --version
```

The script checks required columns, category order, integer counts, and a
total of 496 before rendering. Percentages are `count / 496 × 100`, rounded
to two decimal places. No random seed is needed. Run
`python3 scripts/validate_repository.py` for source-file and category-count checks.
If the total or category labels differ, confirm the TCGA-LUAD cohort and the
alteration-classification rules before plotting.
