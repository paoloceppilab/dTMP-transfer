# Reference Derivation for WT Amplicons

**Author:** Mert Demirdizen ([mert@bmb.sdu.dk](mailto:mert@bmb.sdu.dk))

## Goal
- Derive the exact wild-type amplicon sequence for `Kras` and `Trp53` in the assay orientation required by CRISPResso2.

## Recommended source of truth
- Primary source: mouse reference assembly from the C57BL/6J strain.
- Validation source: normal-lung FASTQs and other WT-dominant libraries from this run.

## Why C57BL/6J reference
- The mouse reference assembly primary chromosomes are derived from the C57BL/6J strain.
- For current work, use `GRCm38/mm10` as the locked reference build for primer mapping, amplicon extraction, and downstream CRISPResso2 inputs.

## Procedure
1. Map each primer pair to the mouse reference genome.
2. Confirm each primer pair maps uniquely and yields a single plausible amplicon on `GRCm38/mm10`.
3. Extract the genomic sequence from the forward-primer start to the reverse-primer reverse-complement end.
4. Orient the final sequence so it matches the sequenced amplicon/guide orientation used for CRISPResso2.
5. Confirm the guide lies inside the amplicon and infer the SpCas9 cut site as 3 bp upstream of the PAM.
6. Cross-check the candidate WT amplicon against normal-lung reads:
   - the start of many reads should match the amplicon orientation
   - the dominant allele in normal samples should align as WT
7. For `Kras`, create `expected_hdr_amplicon_seq` by introducing the designed donor edits into the WT amplicon sequence.

## QC checks
- Primer pair is unique or overwhelmingly dominant at the intended locus.
- Guide sequence is present in the WT amplicon.
- WT amplicon length is compatible with the sequenced read length and paired-end overlap.
- Normal samples show predominantly WT sequence.
- `Kras` HDR amplicon differs from WT only at the intended donor-edited positions.

## Current status
- WT amplicons for `Kras` and `Trp53` have been filled into `amplicons.tsv`.
- `Kras expected_hdr_amplicon_seq` has been derived empirically from donor-positive tumor reads and added to `amplicons.tsv`.

## HDR derivation note
- The `Kras` HDR donor motif was observed directly in tumor FASTQs on UCloud.
- A reproducible derivation script is stored at `/work/alignment/scripts/derive_kras_expected_hdr_from_reads.py`.
- Saved derivation summaries:
  - `/work/alignment/results/summary/kras_hdr_derivation_L713.json`
  - `/work/alignment/results/summary/kras_hdr_derivation_L718.json`
- In both `L713` and `L718`, the dominant donor-positive 67 bp prefix was identical.
- Because the available read lengths are 67 bp and 51 bp, the full 134 bp HDR amplicon is completed by copying bases 68-134 from the WT amplicon.

## Remaining blocker
- `L742_S32` FASTQs are still missing remotely.
