# Kras and Trp53 amplicon references

**Author:** Mert Demirdizen ([mert@bmb.sdu.dk](mailto:mert@bmb.sdu.dk))

The assay uses GRCm38/mm10 (C57BL/6J) coordinates. The wild-type amplicon,
primer, guide, cut-site, and expected `Kras` HDR sequences are specified in
[amplicons.tsv](amplicons.tsv). These sequences are required to orient the
reads and classify edit outcomes.

Wild-type amplicons were checked by primer mapping and alignment of
WT-dominant normal-lung reads. The expected `Kras` HDR sequence uses the
dominant donor-positive 67-base prefix observed in libraries `L713` and
`L718`; bases
68–134 come from the wild-type amplicon because the available read lengths
were 67 and 51 bases. The [derivation script](../source/derive_kras_expected_hdr_from_reads.py)
documents this rule. The full HDR amplicon beyond the observed prefix is
therefore an explicit construction, not a read-level observation.

Check that each primer pair maps to one plausible amplicon, the guide lies
within it, and the HDR sequence differs from wild type only at designed donor
positions. If these checks fail, do not reuse the reference in CRISPResso.
`KP_811` Trp53 is excluded from quantified read-level results after QC; the
editing summary records its missing value.
