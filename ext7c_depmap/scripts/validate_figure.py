#!/usr/bin/env python3
"""Validate Extended Fig. 7C inputs and reported DepMap figure."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import xml.etree.ElementTree as ET
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote

from PIL import Image
from pypdf import PdfReader

FIGURE_DIR_NAME = "ext7c_depmap"
PACKAGE_REL = Path("publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap")


def fail(message: str) -> None:
    raise SystemExit(message)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def tracked_figure_paths(repo_root: Path) -> set[str]:
    figure_dir = Path(FIGURE_DIR_NAME)
    if not (repo_root / ".git").exists():
        source_root = repo_root / figure_dir
        return {
            path.relative_to(source_root).as_posix()
            for path in source_root.rglob("*")
            if path.is_file() and path.name != "provenance_manifest.json"
        }
    output = subprocess.check_output(
        ["git", "ls-files", "--cached", "--", FIGURE_DIR_NAME],
        cwd=repo_root,
        text=True,
    )
    return {
        Path(path).relative_to(figure_dir).as_posix()
        for path in output.splitlines()
        if path != (figure_dir / "provenance_manifest.json").as_posix()
    }


def validate_provenance(repo_root: Path, figure_dir: Path) -> dict:
    manifest = json.loads((figure_dir / "provenance_manifest.json").read_text())
    records = manifest["source_files"]
    relative_paths = [rec["relative_local_path"] for rec in records]
    if manifest.get("source_file_count") != len(records):
        fail("provenance source_file_count does not match its source list")
    if len(relative_paths) != len(set(relative_paths)):
        fail("duplicate provenance paths")
    if set(relative_paths) != tracked_figure_paths(repo_root):
        fail("provenance/tracked-file coverage mismatch")
    for record in records:
        relative = Path(record["relative_local_path"])
        if relative.is_absolute() or ".." in relative.parts:
            fail(f"unsafe provenance path: {relative}")
        path = figure_dir / relative
        if not path.is_file():
            fail(f"missing provenance file: {path}")
    return manifest


def validate_depmap_tables(package: Path, run: dict) -> None:
    genes = run["genes"]
    cell_lines = run["cell_lines"]
    if len(genes) != 17 or len(set(genes)) != 17:
        fail("expected 17 distinct connexin genes")
    if len(cell_lines) != 8 or len(set(cell_lines)) != 8:
        fail("expected eight distinct cell lines")
    if run["required_gjd3_column"] != "GJD3 (ENSG00000183153)":
        fail("required GJD3 Ensembl match changed")
    header, tpm_rows = read_tsv(package / "tables/depmap_24q4_tpm_matrix.tsv")
    log_header, log_rows = read_tsv(package / "tables/depmap_24q4_log2_tpm_plus1_matrix.tsv")
    expected_header = ["cell_line", *genes]
    if header != expected_header or log_header != expected_header:
        fail("DepMap matrix columns/order mismatch")
    if len(tpm_rows) != 8 or len(log_rows) != 8:
        fail("DepMap matrix must be 8 x 17")
    if [row["cell_line"] for row in tpm_rows] != cell_lines:
        fail("TPM cell-line order/mapping mismatch")
    if [row["cell_line"] for row in log_rows] != cell_lines:
        fail("source log2 matrix cell-line order mismatch")
    maximum = 0.0
    for tpm_row, log_row in zip(tpm_rows, log_rows):
        for gene in genes:
            try:
                tpm = float(tpm_row[gene])
                log_value = float(log_row[gene])
            except (KeyError, TypeError, ValueError):
                fail(f"missing/non-numeric value: {tpm_row['cell_line']} {gene}")
            if not all(math.isfinite(value) and value >= 0 for value in (tpm, log_value)):
                fail(f"non-finite or negative value: {tpm_row['cell_line']} {gene}")
            if abs(tpm - (2**log_value - 1)) > 0.0051:
                fail(f"inverse-transform mismatch: {tpm_row['cell_line']} {gene}")
            maximum = max(maximum, tpm)
    primary_scale = run["color_scales"]["primary_tpm_native"]
    if primary_scale["vmin"] != 0 or not math.isclose(primary_scale["vmax"], maximum, abs_tol=1e-9):
        fail("primary TPM color scale does not match matrix maximum")
    _, genes_provenance = read_tsv(package / "tables/gene_availability_provenance.tsv")
    gjd3 = [row for row in genes_provenance if row["gene"] == "GJD3"]
    if len(gjd3) != 1 or gjd3[0]["depmap_column_id"] != run["required_gjd3_column"]:
        fail("GJD3 gene provenance mismatch")
    _, mapping = read_tsv(package / "tables/sample_profile_mapping.tsv")
    if len(mapping) != 8 or [row["cell_line"] for row in mapping] != cell_lines:
        fail("default RNA profile mapping mismatch")
    if len({row["depmap_rna_profile_id"] for row in mapping}) != 8:
        fail("RNA profile IDs must be unique")


def validate_outputs(repo_root: Path, package: Path, run: dict) -> None:
    if run["project_root"] != ".":
        fail("run manifest project_root must be repository relative")
    script = repo_root / run["script_path"]
    if not script.is_file() or not script.is_relative_to(package):
        fail("run manifest script_path is invalid")
    if "ext7c_depmap/publication/" not in run["rerun_command_from_project_root"]:
        fail("rerun command points to a stale figure path")

    def check_paths(value: object) -> None:
        if isinstance(value, dict):
            for item in value.values():
                check_paths(item)
        elif isinstance(value, list):
            for item in value:
                check_paths(item)
        elif isinstance(value, str):
            path = Path(value)
            if path.is_absolute() or ".." in path.parts or not (repo_root / path).is_file():
                fail(f"invalid run manifest output path: {value}")

    check_paths(run["outputs"])
    figures = sorted((package / "figures").iterdir())
    if len(figures) != 6:
        fail(f"expected six figure exports, found {len(figures)}")
    formats = {".png": 0, ".pdf": 0, ".svg": 0}
    for figure in figures:
        if figure.stat().st_size == 0 or figure.suffix not in formats:
            fail(f"empty or unexpected figure export: {figure}")
        formats[figure.suffix] += 1
        if figure.suffix == ".png":
            with Image.open(figure) as image:
                image.verify()
        elif figure.suffix == ".pdf":
            if len(PdfReader(figure).pages) != 1:
                fail(f"PDF page count mismatch: {figure}")
        elif not ET.parse(figure).getroot().tag.endswith("svg"):
            fail(f"invalid SVG: {figure}")
    if set(formats.values()) != {2}:
        fail(f"expected two exports per format: {formats}")


class LocalLinks(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.targets: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        for name, value in attrs:
            if name in {"href", "src"} and value:
                self.targets.append(value)


def validate_links(repo_root: Path, figure_dir: Path, package: Path) -> None:
    markdown_files = [
        repo_root / "README.md",
        repo_root / "ext4_fig5_rnaseq/README.md",
        figure_dir / "README.md",
        package / "README.md",
    ]
    for markdown in markdown_files:
        text = markdown.read_text()
        for target in re.findall(r"\]\(([^)]+)\)", text):
            target = unquote(target.strip("<> ").split("#", 1)[0])
            if not target or "://" in target or target.startswith(("mailto:", "#")):
                continue
            if not (markdown.parent / target).exists():
                fail(f"broken local Markdown link: {markdown}: {target}")
    index = package / "index.html"
    links = LocalLinks()
    links.feed(index.read_text())
    for target in links.targets:
        if "://" not in target and not (package / unquote(target)).exists():
            fail(f"broken package index link: {target}")
    stale_prefix = "ext4_fig5_rnaseq/publication/extended_fig7c_depmap_24q4_tpm_connexin_heatmap"
    for path in package.rglob("*"):
        if path.is_file() and path.suffix in {".md", ".json", ".py", ".html", ".tsv"}:
            if stale_prefix in path.read_text(errors="replace"):
                fail(f"stale active package path: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root", type=Path, default=Path(__file__).resolve().parents[2],
        help="Git repository root (default: inferred from this script)",
    )
    args = parser.parse_args()
    repo_root = args.repo_root.resolve()
    figure_dir = repo_root / FIGURE_DIR_NAME
    package = figure_dir / PACKAGE_REL
    manifest = validate_provenance(repo_root, figure_dir)
    run = json.loads((package / "run_manifest.json").read_text())
    validate_depmap_tables(package, run)
    validate_outputs(repo_root, package, run)
    validate_links(repo_root, figure_dir, package)
    print(f"OK Extended Fig. 7C: {len(manifest['source_files'])} source files, 8x17 TPM, six figures, and links")


if __name__ == "__main__":
    main()
