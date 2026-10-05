# DepMap 24Q4 TPM Methods

Connexin mRNA abundance was derived from the official [DepMap 24Q4 Public](https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248) release. The all-gene RNA expression file was `OmicsExpressionAllGenesTPMLogp1Profile.csv` (Figshare file ID `51065360`), and the default model/profile map was `OmicsDefaultModelProfiles.csv` (Figshare file ID `51065339`). Both raw files were validated against pinned sizes and SHA-256 hashes before processing.

The source expression values are `log2(TPM + 1)`. Values were inverse-transformed using `TPM = 2 ** source_value - 1` before saving the plotted table or drawing the heatmap. No RPKM or protein-coding-only matrix was used for the primary figure.

The predefined lung cancer cell-line order was A549, Calu-1, NCI-H23, SK-MES-1, NCI-H520, NCI-H1299, BEN, and NCI-H838. Supplied DepMap model IDs were matched to default RNA profile IDs using `ProfileType == "RNA"` in `OmicsDefaultModelProfiles.csv`. The gene order was GJA1, GJA10, GJA3, GJA5, GJA8, GJA9, GJB2, GJB4, GJB5, GJB6, GJB7, GJC1, GJC2, GJC3, GJD2, GJD3, and GJD4. Gene columns were checked by HGNC symbol and Ensembl gene ID; GJD3 was required to match `GJD3 (ENSG00000183153)` with no gene substitution.

The primary heatmap uses a TPM-native scale from 0 to 91.71 TPM, the maximum of the plotted DepMap matrix, so no plotted DepMap value is clipped. QC requires all eight default RNA profiles and 17 genes to map uniquely, an 8 × 17 finite nonnegative TPM matrix without missing values, valid source checksums, and parseable PNG/PDF/SVG exports. This is a descriptive cell-line heatmap; no replicate-level statistical test or multiple-testing inference is reported.
