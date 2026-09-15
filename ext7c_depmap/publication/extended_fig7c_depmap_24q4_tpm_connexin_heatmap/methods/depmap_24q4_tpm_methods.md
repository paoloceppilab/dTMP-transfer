# DepMap 24Q4 TPM Methods

Connexin mRNA abundance was re-derived from the official DepMap 24Q4 Public Figshare+ release (https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248). The all-gene RNA expression source file was `OmicsExpressionAllGenesTPMLogp1Profile.csv` (Figshare file `51065360`), and the default model/profile map was `OmicsDefaultModelProfiles.csv` (Figshare file `51065339`).

DepMap expression values were treated as `log2(TPM + 1)`. Values were converted to TPM before tabulation or plotting using `TPM = 2 ** log2_tpm_plus_1 - 1`. No RPKM matrix and no protein-coding-only TPM matrix was used.

The eight lung cancer cell lines were plotted in the predefined order: A549, Calu-1, NCI-H23, SK-MES-1, NCI-H520, NCI-H1299, BEN, NCI-H838. Cell lines were mapped to DepMap model IDs using the supplied identifiers, then to default RNA profile IDs using rows with `ProfileType == "RNA"` in `OmicsDefaultModelProfiles.csv`.

The 17 connexin genes were plotted in the predefined order: GJA1, GJA10, GJA3, GJA5, GJA8, GJA9, GJB2, GJB4, GJB5, GJB6, GJB7, GJC1, GJC2, GJC3, GJD2, GJD3, GJD4. Gene columns were matched by HGNC symbol and Ensembl gene ID. `GJD3` was required to match `GJD3 (ENSG00000183153)`; no `GJC1` substitution was permitted.

The primary heatmap uses a TPM-native scale with `vmin=0` and `vmax=91.71`, the observed maximum in the DepMap 24Q4 TPM matrix. Its DepMap values are not clipped.

Quality control required all 17 genes to be present, all 8 cell lines to map to default RNA profiles, no missing values in the final 8 x 17 TPM matrix, successful PNG/PDF/SVG parsing, and checksum validation of both DepMap source files.
