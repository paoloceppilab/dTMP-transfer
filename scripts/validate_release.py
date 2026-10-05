#!/usr/bin/env python3
"""Validate publication metadata, mouse RNA-seq sources and release integrity.

Default operation is read-only. --write-manifest explicitly creates the
release file inventory before freezing the publication commit.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import re
import subprocess
from collections import Counter
from pathlib import Path
from urllib.parse import unquote, urlsplit

import openpyxl
import yaml

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = Path("release_manifest.json")
EXPECTED_CREATORS = [
    ("Paolo", "Ceppi"), ("Mohammad A.", "Siddiqui"),
    ("Anne Mette A.", "Rømer"), ("Vignesh", "Ramesh"), ("Mert", "Demirdizen"),
]
ANALYSIS_DIRS = {
    "fig2a_ext3_scrnaseq", "fig2b_d_ext3f_visium", "ext3c_tcga",
    "ext4_fig5_rnaseq", "ext7c_depmap", "fig4j_cellstates",
    "ext10_amplicon", "mouse_kp_kptt_rnaseq",
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def files() -> list[Path]:
    if (ROOT / ".git").exists():
        raw = subprocess.check_output(["git", "ls-files", "--cached", "-z"], cwd=ROOT)
        return sorted(Path(x.decode()) for x in raw.split(b"\0") if x and Path(x.decode()) != MANIFEST)
    excluded = {".git", ".venv", ".venv-validation", "__pycache__", ".pytest_cache"}
    return sorted(p.relative_to(ROOT) for p in ROOT.rglob("*")
                  if p.is_file() and p.relative_to(ROOT) != MANIFEST
                  and not excluded.intersection(p.relative_to(ROOT).parts)
                  and p.suffix not in {".pyc", ".pyo"} and p.name != ".DS_Store")


def fail(message: str) -> None:
    raise SystemExit("FAIL: " + message)


def safe_path(path: Path) -> None:
    if path.is_absolute() or ".." in path.parts:
        fail(f"unsafe inventory path: {path}")
    if path.parts[0] not in ANALYSIS_DIRS | {"scripts", "env"} and len(path.parts) != 1:
        fail(f"unexpected directory in release: {path}")
    if path.suffix.lower() in {".docx", ".bam", ".bai", ".cram", ".sam", ".h5ad", ".h5", ".rds", ".pt", ".pth", ".pem", ".key"}:
        fail(f"excluded file type in release: {path}")
    if path.name.lower().endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        fail(f"raw sequencing file in release: {path}")
    if path.name.startswith(".env"):
        fail(f"local credential/configuration file in release: {path}")


def check_metadata() -> None:
    zenodo = json.loads((ROOT / ".zenodo.json").read_text())
    cff = yaml.safe_load((ROOT / "CITATION.cff").read_text())
    expected = [f"{last}, {first}" for first, last in EXPECTED_CREATORS]
    if [c["name"] for c in zenodo["creators"]] != expected:
        fail("Zenodo creators/order do not match the agreed five authors")
    if [(c["given-names"], c["family-names"]) for c in cff["authors"]] != EXPECTED_CREATORS:
        fail("CITATION.cff authors/order mismatch")
    if zenodo["title"] != cff["title"] or str(zenodo["version"]) != "1.0.0" or str(cff["version"]) != "1.0.0":
        fail("citation title/version mismatch")
    if zenodo["upload_type"] != "software" or zenodo["license"] != "mit" or cff["license"] != "MIT":
        fail("software resource type or MIT metadata mismatch")
    if "MIT License" not in (ROOT / "LICENSE").read_text() or not (ROOT / "THIRD_PARTY_NOTICES.md").is_file():
        fail("license text or third-party attribution missing")
    logging.info("metadata: five ordered creators, title, version and license agree")


def check_mouse() -> None:
    directory = ROOT / "mouse_kp_kptt_rnaseq"
    manifest = json.loads((directory / "provenance_manifest.json").read_text())
    if len(manifest["source_files"]) != manifest["source_file_count"]:
        fail("mouse manifest count mismatch")
    recorded = {rec["relative_local_path"] for rec in manifest["source_files"]}
    actual = {p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file() and p.name != "provenance_manifest.json"}
    if recorded != actual:
        fail("mouse source inventory coverage mismatch")
    for record in manifest["source_files"]:
        p = directory / record["relative_local_path"]
        if digest(p) != record["sha256"]:
            fail(f"mouse source checksum mismatch: {p.relative_to(ROOT)}")
    for filename, expected_counts in [("all_samples_si.txt", {"KP": 3, "KPTT": 5}), ("three_samples_si.txt", {"KP": 3, "KPTT": 3})]:
        with (directory / "metadata/deseq2_sample_info" / filename).open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if len({r["Samples"] for r in rows}) != len(rows) or Counter(r["Type"] for r in rows) != Counter(expected_counts):
            fail(f"mouse sample IDs/group counts mismatch: {filename}")
    text = (directory / "source/24_01_16_mice_KPTT_Script.txt").read_text()
    for term in manifest["qc_expectations"]["workflow_script_required_terms"]:
        if term not in text:
            fail(f"mouse source workflow term missing: {term}")
    book = openpyxl.load_workbook(directory / manifest["qc_expectations"]["deseq2_workbook"]["relative_local_path"], read_only=True, data_only=True)
    try:
        for sheet, expected in manifest["qc_expectations"]["deseq2_workbook"]["expected_rows"].items():
            if sheet not in book or book[sheet].max_row != expected:
                fail(f"mouse workbook rows mismatch: {sheet}")
    finally:
        book.close()
    logging.info("mouse RNA-seq: %d sources, eight/six samples and workbook dimensions verified", len(recorded))


def check_links(paths: list[Path]) -> None:
    count = 0
    for rel in paths:
        if rel.suffix != ".md":
            continue
        for target in re.findall(r"\]\(([^)]+)\)", (ROOT / rel).read_text()):
            target = target.strip("<> ").split(" ", 1)[0]
            if not target or urlsplit(target).scheme or target.startswith("#"):
                continue
            target = unquote(target.split("#", 1)[0])
            if not (ROOT / rel.parent / target).exists():
                fail(f"broken relative link in {rel}: {target}")
            count += 1
    logging.info("relative Markdown links: %d targets verified", count)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    paths = files()
    for rel in paths:
        safe_path(rel)
        if not (ROOT / rel).is_file():
            fail(f"missing release file: {rel}")
    if args.write_manifest:
        records = [{"path": p.as_posix(), "size_bytes": (ROOT / p).stat().st_size, "sha256": digest(ROOT / p)} for p in paths]
        data = {"schema_version": "1.0", "release_version": "v1.0.0", "hash_algorithm": "SHA256", "self_exclusion": MANIFEST.as_posix(), "files": records}
        (ROOT / MANIFEST).write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
    manifest = json.loads((ROOT / MANIFEST).read_text())
    records = manifest["files"]
    if len({r["path"] for r in records}) != len(records) or set(r["path"] for r in records) != {p.as_posix() for p in paths}:
        fail("release inventory coverage or uniqueness mismatch")
    for rec in records:
        rel = Path(rec["path"])
        safe_path(rel)
        p = ROOT / rel
        if p.stat().st_size != rec["size_bytes"] or digest(p) != rec["sha256"]:
            fail(f"release file integrity mismatch: {rel}")
    check_metadata()
    check_mouse()
    check_links(paths)
    logging.info("release integrity: %d files, %d bytes; all SHA-256 checks pass", len(records), sum(r["size_bytes"] for r in records))


if __name__ == "__main__":
    main()
