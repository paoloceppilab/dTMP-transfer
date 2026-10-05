#!/usr/bin/env python3
"""Validate tracked figure files and key panel invariants.

Run with Python >=3.9 from any directory; standard library only.
"""

import csv
import json
import math
import re
import subprocess
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

root = Path(__file__).resolve().parents[1]

def _tracked_paths(figure_dir_name):
    figure_dir = root / figure_dir_name
    if (root / ".git").exists():
        try:
            raw_paths = subprocess.check_output(
                ["git", "ls-files", "--cached", "-z", "--", figure_dir_name],
                cwd=root,
            )
        except (OSError, subprocess.CalledProcessError) as exc:
            raise SystemExit("cannot list tracked figure files with Git") from exc
        paths = set()
        for raw in raw_paths.split(b"\0"):
            if not raw:
                continue
            path = Path(raw.decode("utf-8", "surrogateescape"))
            if path == Path(figure_dir_name) / "provenance_manifest.json":
                continue
            paths.add(path.relative_to(figure_dir_name).as_posix())
        return paths
    return {
        path.relative_to(figure_dir).as_posix()
        for path in figure_dir.rglob("*")
        if path.is_file() and path.name != "provenance_manifest.json"
    }

def _xlsx_col_to_int(col):
    value = 0
    for char in col:
        value = value * 26 + ord(char) - ord("A") + 1
    return value

def _xlsx_dimensions(path):
    main_ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    rel_ns = "http://schemas.openxmlformats.org/package/2006/relationships"
    office_rel_ns = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    ns = {"main": main_ns, "rel": rel_ns}
    dimensions = {}
    with zipfile.ZipFile(path) as workbook_zip:
        workbook = ET.fromstring(workbook_zip.read("xl/workbook.xml"))
        rels = ET.fromstring(workbook_zip.read("xl/_rels/workbook.xml.rels"))
        rel_targets = {
            rel.attrib["Id"]: rel.attrib["Target"]
            for rel in rels.findall("rel:Relationship", ns)
        }
        for sheet in workbook.findall("main:sheets/main:sheet", ns):
            sheet_name = sheet.attrib["name"]
            rel_id = sheet.attrib[f"{{{office_rel_ns}}}id"]
            target = rel_targets[rel_id]
            sheet_path = target.lstrip("/") if target.startswith("/") else f"xl/{target}"
            sheet_xml = ET.fromstring(workbook_zip.read(sheet_path))
            dimension = sheet_xml.find("main:dimension", ns)
            ref = dimension.attrib.get("ref", "") if dimension is not None else ""
            match = re.search(r"([A-Z]+)([0-9]+)$", ref.split(":")[-1])
            if not match:
                max_row = 0
                max_col = 0
                for cell in sheet_xml.findall(".//main:c", ns):
                    cell_ref = cell.attrib.get("r", "")
                    cell_match = re.match(r"([A-Z]+)([0-9]+)", cell_ref)
                    if cell_match:
                        max_col = max(max_col, _xlsx_col_to_int(cell_match.group(1)))
                        max_row = max(max_row, int(cell_match.group(2)))
            else:
                max_col = _xlsx_col_to_int(match.group(1))
                max_row = int(match.group(2))
            dimensions[sheet_name] = (max_row, max_col)
    return dimensions

def _assert_xlsx(path, required_sheets, expected_rows=None):
    dimensions = _xlsx_dimensions(path)
    missing = set(required_sheets) - set(dimensions)
    if missing:
        raise SystemExit(f"missing workbook sheets in {path}: {sorted(missing)}")
    for sheet_name, expected_row_count in (expected_rows or {}).items():
        observed_row_count = dimensions[sheet_name][0]
        if observed_row_count != expected_row_count:
            raise SystemExit(
                f"unexpected row count in {path}:{sheet_name}: "
                f"{observed_row_count}"
            )

expected_manifest_record_counts = {
    "fig4j_cellstates": 8,
    "fig2b_d_ext3f_visium": 10,
    "fig2a_ext3_scrnaseq": 60,
    "ext4_fig5_rnaseq": 15,
    "ext10_amplicon": 52,
}

for figure_dir_name in (
    "fig4j_cellstates",
    "fig2b_d_ext3f_visium",
    "fig2a_ext3_scrnaseq",
    "ext4_fig5_rnaseq",
    "ext10_amplicon",
):
    figure_dir = root / figure_dir_name
    manifest = json.loads((figure_dir / "provenance_manifest.json").read_text())
    observed_count = len(manifest["source_files"])
    expected_count = expected_manifest_record_counts[figure_dir_name]
    if manifest.get("source_file_count") != observed_count:
        raise SystemExit(f"{figure_dir_name} source_file_count does not match its source list")
    if observed_count != expected_count:
        raise SystemExit(
            f"{figure_dir_name} source-record count changed: "
            f"expected {expected_count}, observed {observed_count}"
        )
    manifest_paths = {rec["relative_local_path"] for rec in manifest["source_files"]}
    if len(manifest_paths) != observed_count:
        raise SystemExit(f"{figure_dir_name} has duplicate source paths")
    if figure_dir_name in {
        "ext4_fig5_rnaseq",
        "ext10_amplicon",
    }:
        tracked_paths = _tracked_paths(figure_dir_name)
        if manifest_paths != tracked_paths:
            raise SystemExit(
                f"{figure_dir_name} manifest does not cover every "
                "non-self source file"
            )
    for rec in manifest["source_files"]:
        relative = Path(rec["relative_local_path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise SystemExit(f"unsafe source path: {figure_dir_name}/{relative}")
        path = figure_dir / relative
        if not path.is_file():
            raise SystemExit(f"missing source file: {figure_dir_name}/{relative}")
    print(f"OK {figure_dir_name}: {observed_count} source files present")

spatial_dir = root / "fig2b_d_ext3f_visium"
with (spatial_dir / "results/tyms_tk1_counts.csv").open(newline="") as handle:
    counts = {row["Category"]: int(row["Counts"]) for row in csv.DictReader(handle)}
expected_counts = {
    "TYMS_pos_TK1_pos": 1707,
    "TYMS_pos_TK1_neg": 549,
    "TYMS_neg_TK1_neg": 1273,
    "TYMS_neg_TK1_pos": 329,
}
if counts != expected_counts:
    raise SystemExit(f"unexpected Visium TYMS/TK1 counts: {counts}")

with (spatial_dir / "results/neighborhood_data_points.csv").open(newline="") as handle:
    rows = list(csv.DictReader(handle))
if len(rows) != 160:
    raise SystemExit(f"unexpected Visium neighborhood row count: {len(rows)}")
radii = sorted({int(float(row["radius"])) for row in rows})
if radii != list(range(1, 11)):
    raise SystemExit(f"unexpected Visium radii: {radii}")
for row in rows:
    radius = int(float(row["radius"]))
    expected_neighbors = {
        1: 6,
        2: 18,
        3: 36,
        4: 60,
        5: 90,
        6: 138,
        7: 186,
        8: 240,
        9: 300,
        10: 366,
    }[radius]
    observed_neighbors = int(float(row["numneigh"]))
    if observed_neighbors != expected_neighbors:
        raise SystemExit(
            f"unexpected Visium numneigh at radius {radius}: {observed_neighbors}"
        )

ext10_dir = root / "ext10_amplicon"
with (ext10_dir / "tables/extended_10b_editing_summary.tsv").open(newline="") as handle:
    ext10b_rows = list(csv.DictReader(handle, delimiter="\t"))
if len(ext10b_rows) != 17:
    raise SystemExit(f"unexpected Extended Data Fig. 10B specimen count: {len(ext10b_rows)}")
required_10b_fields = {
    "specimen_id",
    "cohort",
    "kras_hdr_exact_pct",
    "trp53_indel_pct",
    "trp53_qc_flag",
}
missing_10b_fields = required_10b_fields - set(ext10b_rows[0])
if missing_10b_fields:
    raise SystemExit(f"missing Extended Data Fig. 10B fields: {sorted(missing_10b_fields)}")
trp53_values = {row["specimen_id"]: row["trp53_indel_pct"] for row in ext10b_rows}
trp53_flags = {row["specimen_id"]: row["trp53_qc_flag"] for row in ext10b_rows}
if trp53_values.get("KP_811") != "NA":
    raise SystemExit(f"KP_811 Trp53 should be NA: {trp53_values.get('KP_811')}")
if trp53_flags.get("KP_811") != "missing_after_qc":
    raise SystemExit(f"unexpected KP_811 Trp53 flag: {trp53_flags.get('KP_811')}")
quantified_trp53 = [row for row in ext10b_rows if row["trp53_indel_pct"] != "NA"]
if len(quantified_trp53) != 16:
    raise SystemExit(f"unexpected quantified Trp53 sample count: {len(quantified_trp53)}")

threebin = json.loads(
    (
        ext10_dir
        / "metadata/kras_cutsite_pm20_threebin_background83_tumor_only_plot_run_info.json"
    ).read_text()
)
if threebin["aggregate_qc"]["plot_sample_count"] != 15:
    raise SystemExit("Extended Data Fig. 10D tumor-only sample count is not 15")
if threebin["aggregate_qc"]["displayed_nonhdr_classifiable_pairs"] != 8879020:
    raise SystemExit("Extended Data Fig. 10D non-HDR denominator is not 8,879,020")

frame = json.loads(
    (
        ext10_dir
        / "metadata/kras_cutsite_pm20_indel_frame_background83_tumor_only_plot_run_info.json"
    ).read_text()
)
frame_pct = frame["aggregate_percentages_of_displayed_nonhdr_classifiable"]
expected_frame_pct = {
    "frameshift_indel": 19.326986536802487,
    "in_frame_indel": 6.009582138569346,
    "other_indel_frame_status": 0.0,
}
for key, expected in expected_frame_pct.items():
    observed = frame_pct[key]
    if not math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-9):
        raise SystemExit(f"unexpected Extended Data Fig. 10D {key}: {observed}")

hotspot = json.loads(
    (
        ext10_dir
        / "metadata/kras_g12_g13_hotspot_frequency_tumor_only_nonhdr_run_info.json"
    ).read_text()
)
if hotspot["corrected_nonhdr"]["aggregate_denominator_allele_reads"] != 19487513:
    raise SystemExit("Extended Data Fig. 10E/F corrected denominator is not 19,487,513")
if hotspot["hdr_positive_removed"]["tumor_strict_hdr_barcode_yes_reads"] != 2594:
    raise SystemExit("Extended Data Fig. 10E/F HDR-positive removed read count is not 2,594")
if hotspot["corrected_nonhdr"]["aggregate_counts"]["G12D"] != 12680:
    raise SystemExit("Extended Data Fig. 10E/F corrected G12D read count is not 12,680")
g12d_pct = hotspot["corrected_nonhdr"]["aggregate_percentages"]["G12D"]
if not math.isclose(g12d_pct, 0.06506730745991035, rel_tol=0.0, abs_tol=1e-12):
    raise SystemExit(f"unexpected Extended Data Fig. 10E/F corrected G12D percent: {g12d_pct}")
with (
    ext10_dir
    / "tables/extended_10e_10f_kras_g12_g13_hotspot_frequency_tumor_only_nonhdr_figure_ready.tsv"
).open(newline="") as handle:
    hotspot_rows = list(csv.DictReader(handle, delimiter="\t"))
if len(hotspot_rows) != 15:
    raise SystemExit(f"unexpected Extended Data Fig. 10E/F sample count: {len(hotspot_rows)}")
for row in hotspot_rows:
    category_sum = float(row["category_pct_sum"])
    if not math.isclose(category_sum, 100.0, rel_tol=0.0, abs_tol=1e-9):
        raise SystemExit(
            f"Extended Data Fig. 10E/F category percentages do not sum to 100 for "
            f"{row['specimen_id']}: {category_sum}"
        )

rnaseq_dir = root / "ext4_fig5_rnaseq"
_assert_xlsx(
    rnaseq_dir / "tables/figure_5a_rnaseq_degs_source.xlsx",
    ["Raw_Counts", "SI", "TS vs NTC", "TScc vs NTCcc", "TScc vs TS"],
    {
        "Raw_Counts": 61542,
        "TS vs NTC": 61542,
        "TScc vs NTCcc": 61542,
        "TScc vs TS": 61542,
    },
)
_assert_xlsx(
    rnaseq_dir / "tables/upstream_deseq2_results.xlsx",
    ["Feature_Counts", "Raw_Counts", "SI", "TS_NTC", "TScc_NTCcc", "TScc_TS"],
    {
        "Feature_Counts": 61542,
        "Raw_Counts": 61542,
        "TS_NTC": 61542,
        "TScc_NTCcc": 61542,
        "TScc_TS": 61542,
    },
)
script_text = (rnaseq_dir / "source/22_08_23_Script_TS.txt").read_text()
required_script_terms = [
    "fastqc *.fastq.gz",
    "fastqc_v0.11.9",
    "Homo_sapiens.GRCh38.105.gtf",
    "featureCounts",
    "library(DESeq2)",
    "results(dds,contrast=c(\"Type\",\"TS\",\"NTC\"))",
    "results(dds,contrast=c(\"Type\",\"TS_CC\",\"NTC_CC\"))",
    "results(dds,contrast=c(\"Type\",\"TS_CC\",\"TS\"))",
]
missing_terms = [term for term in required_script_terms if term not in script_text]
if missing_terms:
    raise SystemExit(f"upstream RNA-seq script is missing terms: {missing_terms}")
for relative_path in _tracked_paths("ext4_fig5_rnaseq"):
    path = rnaseq_dir / relative_path
    lower_name = path.name.lower()
    lower_parts = {part.lower() for part in path.relative_to(rnaseq_dir).parts}
    if lower_name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        raise SystemExit(f"raw FASTQ should not be tracked: {path}")
    if lower_name.endswith((".bam", ".sam", ".cram", ".bai", ".crai")):
        raise SystemExit(f"alignment file should not be tracked: {path}")
    if lower_name in {".ds_store", "thumbs.db"} or lower_name.startswith("~$"):
        raise SystemExit(f"transient system/lock file should not be tracked: {path}")
    if lower_parts & {"genome", "star_index", "star_output", "bam", "fastq"}:
        raise SystemExit(f"STAR index/output or raw-data directory should not be tracked: {path}")
print("OK RNA-seq source files: workbook/script/exclusion checks passed")
