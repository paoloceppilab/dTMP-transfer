#!/usr/bin/env python3
"""Validate publication metadata, Methods scope and release integrity.

Default operation is read-only. --write-manifest explicitly creates the
release file inventory before freezing the publication commit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import re
import subprocess
from pathlib import Path
from urllib.parse import unquote, urlsplit

import yaml

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = Path("release_manifest.json")
EXPECTED_CREATORS = [('Mert', 'Demirdizen'), ('Vignesh', 'Ramesh'), ('Mohammad A.', 'Siddiqui'), ('Anne Mette A.', 'Rømer'), ('Paolo', 'Ceppi')]
ANALYSIS_DIRS = {'ext10_amplicon', 'ext4_fig5_rnaseq', 'fig2a_ext3_scrnaseq', 'fig2b_d_ext3f_visium', 'fig4j_cellstates'}


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
    if path.suffix.lower() in {".pdf", ".svg", ".png", ".jpg", ".jpeg", ".tif", ".tiff", ".eps", ".ppt", ".pptx", ".pzfx", ".docx", ".bam", ".bai", ".cram", ".sam", ".h5ad", ".h5", ".rds", ".pt", ".pth", ".pem", ".key"}:
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
    if zenodo["title"] != cff["title"] or str(zenodo["version"]) != "1.0.4" or str(cff["version"]) != "1.0.4":
        fail("citation title/version mismatch")
    if zenodo["upload_type"] != "software" or zenodo["license"] != "mit" or cff["license"] != "MIT":
        fail("software resource type or MIT metadata mismatch")
    if "MIT License" not in (ROOT / "LICENSE").read_text() or not (ROOT / "THIRD_PARTY_NOTICES.md").is_file():
        fail("license text or third-party attribution missing")
    logging.info("metadata: five ordered creators, title, version and license agree")


def check_scope() -> None:
    scope = json.loads((ROOT / "analysis_scope.json").read_text())
    expected = {
        "RNA-sequencing": "ext4_fig5_rnaseq",
        "Analysis of Single-cell RNA-sequencing data": "fig2a_ext3_scrnaseq",
        "Neighborhood analysis of Visium spatial sequencing": "fig2b_d_ext3f_visium",
        "Spatial Transcriptomics Deconvolution": "fig4j_cellstates",
        "Amplicon Sequencing": "ext10_amplicon",
    }
    observed = {a["methods_heading"]: a["directory"] for a in scope["analyses"]}
    if observed != expected or len(scope["analyses"]) != 5:
        fail("Methods scope does not match the five supplied sections")
    if scope["release_version"] != "1.0.4" or not re.fullmatch(r"[0-9a-f]{64}", scope["scope_source"]["sha256"]):
        fail("Methods source hash or version is invalid")
    directories = {p.name for p in ROOT.iterdir() if p.is_dir()
                   and p.name not in {".git", "scripts", "env", ".venv", ".venv-validation", "__pycache__", ".pytest_cache"}}
    if directories != ANALYSIS_DIRS:
        fail(f"analysis directories differ from Methods scope: {sorted(directories)}")
    for name in ANALYSIS_DIRS:
        folder = ROOT / name
        manifest = json.loads((folder / "provenance_manifest.json").read_text())
        recorded = {rec["relative_local_path"] for rec in manifest["source_files"]}
        actual = {p.relative_to(folder).as_posix() for p in folder.rglob("*")
                  if p.is_file() and p.name != "provenance_manifest.json"
                  and "__pycache__" not in p.parts and p.suffix not in {".pyc", ".pyo"}}
        if recorded != actual or len(recorded) != manifest["source_file_count"]:
            fail(f"analysis inventory coverage mismatch: {name}")
        for record in manifest["source_files"]:
            if digest(folder / record["relative_local_path"]) != record["sha256"]:
                fail(f"analysis source checksum mismatch: {name}/{record['relative_local_path']}")
    logging.info("Methods scope: exactly five analysis directories; all source inventories and hashes verified")


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
        data = {"schema_version": "1.0", "release_version": "v1.0.4", "hash_algorithm": "SHA256", "self_exclusion": MANIFEST.as_posix(), "files": records}
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
    check_scope()
    check_links(paths)
    logging.info("release integrity: %d files, %d bytes; all SHA-256 checks pass", len(records), sum(r["size_bytes"] for r in records))


if __name__ == "__main__":
    main()
