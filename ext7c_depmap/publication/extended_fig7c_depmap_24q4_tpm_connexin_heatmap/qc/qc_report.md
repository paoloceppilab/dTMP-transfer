# QC Report

## Required Checks

| check | passed |
| --- | --- |
| all_17_genes_present | True |
| gjd3_exact_column_present | True |
| all_8_cell_lines_map_to_default_rna_profiles | True |
| final_tpm_matrix_shape_is_8_by_17 | True |
| final_tpm_matrix_has_no_missing_values | True |
| tpm_native_color_scale_vmin_is_0 | True |
| tpm_native_color_scale_vmax_is_observed_depmap_maximum | True |
| depmap_tpm_values_not_clipped_in_primary_figure | True |
| png_pdf_svg_files_non_empty_and_renderable | True |

## Gene Availability

| gene | present | depmap_column_id | depmap_symbol | depmap_ensembl_id | expected_ensembl_id | match_type |
| --- | --- | --- | --- | --- | --- | --- |
| GJA1 | True | GJA1 (ENSG00000152661) | GJA1 | ENSG00000152661 | ENSG00000152661 | symbol_and_ensembl |
| GJA10 | True | GJA10 (ENSG00000135355) | GJA10 | ENSG00000135355 | ENSG00000135355 | symbol_and_ensembl |
| GJA3 | True | GJA3 (ENSG00000121743) | GJA3 | ENSG00000121743 | ENSG00000121743 | symbol_and_ensembl |
| GJA5 | True | GJA5 (ENSG00000265107) | GJA5 | ENSG00000265107 | ENSG00000265107 | symbol_and_ensembl |
| GJA8 | True | GJA8 (ENSG00000121634) | GJA8 | ENSG00000121634 | ENSG00000121634 | symbol_and_ensembl |
| GJA9 | True | GJA9 (ENSG00000131233) | GJA9 | ENSG00000131233 | ENSG00000131233 | symbol_and_ensembl |
| GJB2 | True | GJB2 (ENSG00000165474) | GJB2 | ENSG00000165474 | ENSG00000165474 | symbol_and_ensembl |
| GJB4 | True | GJB4 (ENSG00000189433) | GJB4 | ENSG00000189433 | ENSG00000189433 | symbol_and_ensembl |
| GJB5 | True | GJB5 (ENSG00000189280) | GJB5 | ENSG00000189280 | ENSG00000189280 | symbol_and_ensembl |
| GJB6 | True | GJB6 (ENSG00000121742) | GJB6 | ENSG00000121742 | ENSG00000121742 | symbol_and_ensembl |
| GJB7 | True | GJB7 (ENSG00000164411) | GJB7 | ENSG00000164411 | ENSG00000164411 | symbol_and_ensembl |
| GJC1 | True | GJC1 (ENSG00000182963) | GJC1 | ENSG00000182963 | ENSG00000182963 | symbol_and_ensembl |
| GJC2 | True | GJC2 (ENSG00000198835) | GJC2 | ENSG00000198835 | ENSG00000198835 | symbol_and_ensembl |
| GJC3 | True | GJC3 (ENSG00000176402) | GJC3 | ENSG00000176402 | ENSG00000176402 | symbol_and_ensembl |
| GJD2 | True | GJD2 (ENSG00000159248) | GJD2 | ENSG00000159248 | ENSG00000159248 | symbol_and_ensembl |
| GJD3 | True | GJD3 (ENSG00000183153) | GJD3 | ENSG00000183153 | ENSG00000183153 | symbol_and_ensembl |
| GJD4 | True | GJD4 (ENSG00000177291) | GJD4 | ENSG00000177291 | ENSG00000177291 | symbol_and_ensembl |

## Sample/Profile Mapping

| cell_line | depmap_model_id | depmap_rna_profile_id |
| --- | --- | --- |
| A549 | ACH-000681 | PR-Em6pVT |
| Calu-1 | ACH-000511 | PR-8led3R |
| NCI-H23 | ACH-000900 | PR-idtmVT |
| SK-MES-1 | ACH-000665 | PR-sVtE0j |
| NCI-H520 | ACH-000395 | PR-4FMaFX |
| NCI-H1299 | ACH-000510 | PR-Er9ReY |
| BEN | ACH-000603 | PR-m1IHlG |
| NCI-H838 | ACH-000416 | PR-AexEWR |

## Figure File Validation

| file | suffix | size_bytes | non_empty | renderable | details |
| --- | --- | --- | --- | --- | --- |
| ext7c_depmap/publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap/figures/extended_fig7c_depmap_24q4_tpm.png | png | 162976 | True | True | {"format": "PNG", "height_px": 2010, "mode": "RGBA", "width_px": 4080} |
| ext7c_depmap/publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap/figures/extended_fig7c_depmap_24q4_tpm.pdf | pdf | 22183 | True | True | {"pages": 1} |
| ext7c_depmap/publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap/figures/extended_fig7c_depmap_24q4_tpm.svg | svg | 19084 | True | True | {"root_tag": "{http://www.w3.org/2000/svg}svg"} |
| ext7c_depmap/publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap/figures/extended_fig7c_depmap_24q4_tpm_no_title.png | png | 137337 | True | True | {"format": "PNG", "height_px": 2010, "mode": "RGBA", "width_px": 4080} |
| ext7c_depmap/publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap/figures/extended_fig7c_depmap_24q4_tpm_no_title.pdf | pdf | 19989 | True | True | {"pages": 1} |
| ext7c_depmap/publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap/figures/extended_fig7c_depmap_24q4_tpm_no_title.svg | svg | 18739 | True | True | {"root_tag": "{http://www.w3.org/2000/svg}svg"} |
