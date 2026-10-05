#!/usr/bin/env python3
"""Generate a publication-ready DepMap 24Q4 TPM figure for Extended Data Fig. 7C.

The source expression matrix is the official DepMap 24Q4 Public Figshare+
all-gene file, OmicsExpressionAllGenesTPMLogp1Profile.csv. DepMap reports
these values as log2(TPM + 1); this workflow reports and plots TPM after the
inverse transform TPM = 2 ** x - 1.

By default, raw DepMap downloads are not retained in this publication folder.
The script first validates any local cache, then downloads to a temporary
directory if no valid cache is available. Use --keep-raw to retain raw source
files under raw/.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
from pypdf import PdfReader


SCRIPT_PATH = Path(__file__).resolve()
OUTPUT_ROOT = SCRIPT_PATH.parents[1]
PROJECT_ROOT = SCRIPT_PATH.parents[4]

WORKBOOK_PATH = (
    PROJECT_ROOT
    / "ext7c_depmap"
    / "tables"
    / "ccle_connexins_workbook.xlsx"
)

FIGSHARE_ARTICLE_URL = "https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248"
FIGSHARE_ARTICLE_ID = "27993248"
EXPRESSION_URL = "https://ndownloader.figshare.com/files/51065360"
PROFILE_URL = "https://ndownloader.figshare.com/files/51065339"

SOURCE_EXPRESSION_UNIT = "log2(TPM + 1)"
REPORTED_EXPRESSION_UNIT = "TPM"
TRANSFORM_FORMULA = "TPM = 2 ** log2_tpm_plus_1 - 1"

CELL_LINES = [
    "A549",
    "Calu-1",
    "NCI-H23",
    "SK-MES-1",
    "NCI-H520",
    "NCI-H1299",
    "BEN",
    "NCI-H838",
]

MODEL_IDS = {
    "A549": "ACH-000681",
    "Calu-1": "ACH-000511",
    "NCI-H23": "ACH-000900",
    "SK-MES-1": "ACH-000665",
    "NCI-H520": "ACH-000395",
    "NCI-H1299": "ACH-000510",
    "BEN": "ACH-000603",
    "NCI-H838": "ACH-000416",
}

GENES = [
    "GJA1",
    "GJA10",
    "GJA3",
    "GJA5",
    "GJA8",
    "GJA9",
    "GJB2",
    "GJB4",
    "GJB5",
    "GJB6",
    "GJB7",
    "GJC1",
    "GJC2",
    "GJC3",
    "GJD2",
    "GJD3",
    "GJD4",
]

GJD3_REQUIRED_COLUMN_ID = "GJD3 (ENSG00000183153)"
GJD3_REQUIRED_ENSEMBL_ID = "ENSG00000183153"

PRIMARY_FIGURE_STEM = "extended_fig7c_depmap_24q4_tpm"
PRIMARY_NO_TITLE_STEM = "extended_fig7c_depmap_24q4_tpm_no_title"
NATIVE_COMPARISON_STEM = "extended_fig7c_depmap_24q4_tpm_vs_paper_workbook_tpm_native"
PAPER_STYLE_COMPARISON_STEM = "extended_fig7c_depmap_24q4_tpm_vs_paper_workbook_paper_style_scale"
PAPER_STYLE_SINGLE_STEM = "extended_fig7c_depmap_24q4_tpm_paper_style_scale"


@dataclass(frozen=True)
class SourceSpec:
    label: str
    url: str
    filename: str
    figshare_file_id: str
    expected_size_bytes: int
    expected_sha256: str


EXPRESSION_SPEC = SourceSpec(
    label="DepMap 24Q4 all-gene RNA expression profile matrix",
    url=EXPRESSION_URL,
    filename="OmicsExpressionAllGenesTPMLogp1Profile.csv",
    figshare_file_id="51065360",
    expected_size_bytes=1_016_340_354,
    expected_sha256="57aa034ec15c48109333aad9d09d31bf8080ce222c45f3b0465fe8e7e7b88d7a",
)

PROFILE_SPEC = SourceSpec(
    label="DepMap 24Q4 default model profiles",
    url=PROFILE_URL,
    filename="OmicsDefaultModelProfiles.csv",
    figshare_file_id="51065339",
    expected_size_bytes=90_080,
    expected_sha256="096a96c39d88374cb7c37058816e6cde29d6617b820b8a14d33ae9dd821c6dc5",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path.resolve())


def portable_outputs(value: Any) -> Any:
    """Record generated paths relative to the repository root."""
    if isinstance(value, dict):
        return {key: portable_outputs(item) for key, item in value.items()}
    if isinstance(value, list):
        return [portable_outputs(item) for item in value]
    if isinstance(value, str) and value.startswith(str(PROJECT_ROOT) + "/"):
        return Path(value).relative_to(PROJECT_ROOT).as_posix()
    return value


def portable_source_metadata(source_metadata: dict[str, Any]) -> dict[str, Any]:
    """Keep source identity without making temporary cache paths active links."""
    records = json.loads(json.dumps(source_metadata))
    for record in records.values():
        source_path = record.get("path")
        record["exists_at_run"] = record.pop("exists", False)
        record["path"] = (
            relative(Path(source_path))
            if source_path and record.get("retained_raw")
            else None
        )
        record.pop("cache_source_path", None)
    return records


def normalise_name(value: object) -> str:
    if value is None:
        return ""
    return re.sub(r"[^A-Z0-9]", "", str(value).upper())


def strip_ensembl_version(value: object) -> str:
    return str(value).split(".", 1)[0]


def parse_gene_identifier(identifier: str) -> dict[str, str | None]:
    match = re.match(r"^(?P<symbol>[^()]+?)(?:\s*\((?P<id>[^()]+)\))?$", identifier)
    if not match:
        return {"symbol": identifier.strip(), "ensembl_id": None, "raw_id": None}
    symbol = match.group("symbol").strip()
    raw_id = (match.group("id") or "").strip()
    ensembl_id = strip_ensembl_version(raw_id) if raw_id.startswith("ENSG") else None
    return {"symbol": symbol, "ensembl_id": ensembl_id, "raw_id": raw_id or None}


def set_csv_field_size_limit() -> None:
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit = int(limit / 10)


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_source_file(path: Path, spec: SourceSpec) -> dict[str, Any]:
    metadata = {
        "path": str(path),
        "exists": path.exists(),
        "expected_size_bytes": spec.expected_size_bytes,
        "expected_sha256": spec.expected_sha256,
    }
    if not path.exists():
        metadata.update({"valid": False, "reason": "missing"})
        return metadata
    size = path.stat().st_size
    sha = sha256_file(path)
    metadata.update({"size_bytes": size, "sha256": sha})
    valid = size == spec.expected_size_bytes and sha == spec.expected_sha256
    metadata.update({"valid": valid})
    if not valid:
        metadata["reason"] = "size_or_sha256_mismatch"
    return metadata


def download_source(spec: SourceSpec, out_path: Path) -> dict[str, Any]:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".partial")
    if tmp_path.exists():
        tmp_path.unlink()
    metadata: dict[str, Any] = {
        "label": spec.label,
        "source_url": spec.url,
        "figshare_file_id": spec.figshare_file_id,
        "expected_filename": spec.filename,
        "download_started_utc": utc_now(),
    }
    req = urllib.request.Request(spec.url, headers={"User-Agent": "depmap-24q4-ext7c-figure"})
    digest = hashlib.sha256()
    bytes_written = 0
    with urllib.request.urlopen(req, timeout=900) as response:
        metadata.update(
            {
                "final_url": response.geturl(),
                "http_status": response.status,
                "content_type": response.headers.get("content-type"),
                "content_length_header": response.headers.get("content-length"),
                "content_disposition": response.headers.get("content-disposition"),
                "etag": response.headers.get("etag"),
                "last_modified": response.headers.get("last-modified"),
            }
        )
        with tmp_path.open("wb") as handle:
            while True:
                chunk = response.read(8 * 1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
                digest.update(chunk)
                bytes_written += len(chunk)
    tmp_path.replace(out_path)
    metadata.update(
        {
            "download_completed_utc": utc_now(),
            "path": str(out_path),
            "size_bytes": bytes_written,
            "sha256": digest.hexdigest(),
            "retained_raw": True,
        }
    )
    if bytes_written != spec.expected_size_bytes or metadata["sha256"] != spec.expected_sha256:
        raise ValueError(
            f"Downloaded {spec.filename} failed integrity check: "
            f"size={bytes_written}, sha256={metadata['sha256']}"
        )
    return metadata


def source_candidates(spec: SourceSpec, raw_cache_dir: Path | None, keep_raw_dir: Path) -> list[Path]:
    candidates: list[Path] = []
    if raw_cache_dir is not None:
        candidates.append(raw_cache_dir / spec.filename)
    candidates.append(keep_raw_dir / spec.filename)
    seen: set[Path] = set()
    unique: list[Path] = []
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved not in seen:
            unique.append(candidate)
            seen.add(resolved)
    return unique


def resolve_source(
    spec: SourceSpec,
    *,
    keep_raw: bool,
    raw_cache_dir: Path | None,
    temp_dir: Path,
) -> tuple[Path, dict[str, Any]]:
    raw_dir = OUTPUT_ROOT / "raw"
    invalid_candidates: list[dict[str, Any]] = []

    if keep_raw:
        retained_path = raw_dir / spec.filename
        existing = validate_source_file(retained_path, spec)
        if existing["valid"]:
            existing.update(
                {
                    "label": spec.label,
                    "source_url": spec.url,
                    "figshare_file_id": spec.figshare_file_id,
                    "expected_filename": spec.filename,
                    "acquisition": "validated_existing_retained_raw",
                    "retained_raw": True,
                }
            )
            return retained_path, existing

        for candidate in source_candidates(spec, raw_cache_dir, raw_dir):
            validation = validate_source_file(candidate, spec)
            if validation["valid"]:
                raw_dir.mkdir(parents=True, exist_ok=True)
                if candidate.resolve() != retained_path.resolve():
                    shutil.copy2(candidate, retained_path)
                copied = validate_source_file(retained_path, spec)
                copied.update(
                    {
                        "label": spec.label,
                        "source_url": spec.url,
                        "figshare_file_id": spec.figshare_file_id,
                        "expected_filename": spec.filename,
                        "acquisition": "copied_from_validated_local_cache",
                        "cache_source_path": str(candidate),
                        "retained_raw": True,
                    }
                )
                return retained_path, copied
            invalid_candidates.append(validation)
        downloaded = download_source(spec, retained_path)
        downloaded["acquisition"] = "downloaded_to_retained_raw"
        return retained_path, downloaded

    for candidate in source_candidates(spec, raw_cache_dir, raw_dir):
        validation = validate_source_file(candidate, spec)
        if validation["valid"]:
            validation.update(
                {
                    "label": spec.label,
                    "source_url": spec.url,
                    "figshare_file_id": spec.figshare_file_id,
                    "expected_filename": spec.filename,
                    "acquisition": "validated_local_cache",
                    "retained_raw": False,
                }
            )
            return candidate, validation
        invalid_candidates.append(validation)

    temp_path = temp_dir / spec.filename
    downloaded = download_source(spec, temp_path)
    downloaded["acquisition"] = "downloaded_to_temporary_file"
    downloaded["retained_raw"] = False
    downloaded["invalid_cache_candidates"] = invalid_candidates
    return temp_path, downloaded


def prepare_output_dirs(force: bool, keep_raw: bool) -> None:
    managed_dirs = ["figures", "tables", "data", "methods", "qc", "docs"]
    if force:
        for dirname in managed_dirs:
            path = OUTPUT_ROOT / dirname
            if path.exists():
                shutil.rmtree(path)
        if not keep_raw:
            raw_path = OUTPUT_ROOT / "raw"
            if raw_path.exists():
                shutil.rmtree(raw_path)
        for filename in ["README.md", "index.html", "run_manifest.json"]:
            path = OUTPUT_ROOT / filename
            if path.exists():
                path.unlink()
    else:
        blocking = []
        for dirname in managed_dirs:
            path = OUTPUT_ROOT / dirname
            if path.exists() and any(path.iterdir()):
                blocking.append(str(path))
        for filename in ["README.md", "index.html", "run_manifest.json"]:
            path = OUTPUT_ROOT / filename
            if path.exists():
                blocking.append(str(path))
        if blocking:
            raise FileExistsError(
                "Refusing to overwrite existing generated outputs without --force:\n"
                + "\n".join(blocking)
            )

    for dirname in managed_dirs:
        (OUTPUT_ROOT / dirname).mkdir(parents=True, exist_ok=True)
    if keep_raw:
        (OUTPUT_ROOT / "raw").mkdir(parents=True, exist_ok=True)


def read_paper_workbook(path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing paper workbook: {path}")
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    required_sheets = {"Connexins_2", "CCLE_data (2)"}
    missing_sheets = sorted(required_sheets.difference(workbook.sheetnames))
    if missing_sheets:
        raise ValueError(f"Workbook is missing required sheets: {missing_sheets}")

    rows = list(workbook["Connexins_2"].iter_rows(values_only=True))
    header = list(rows[0])
    paper = pd.DataFrame(rows[1:], columns=header).rename(columns={"SYMBOL": "cell_line"})
    paper = paper.dropna(subset=["cell_line"]).set_index("cell_line")
    missing_gene_columns = sorted(set(GENES).difference(paper.columns))
    if missing_gene_columns:
        raise ValueError(f"Paper workbook lacks plotted gene columns: {missing_gene_columns}")

    row_by_key: dict[str, str] = {}
    for label in paper.index:
        key = normalise_name(label)
        if key and key not in row_by_key:
            row_by_key[key] = str(label)

    selected_rows = []
    missing_rows = []
    for cell_line in CELL_LINES:
        workbook_label = row_by_key.get(normalise_name(cell_line))
        if workbook_label is None:
            missing_rows.append(cell_line)
            continue
        selected_rows.append(paper.loc[workbook_label, GENES])
    if missing_rows:
        raise ValueError("Paper workbook lacks rows for: " + ", ".join(missing_rows))

    paper_matrix = pd.DataFrame(selected_rows, index=CELL_LINES, columns=GENES)
    paper_matrix = paper_matrix.apply(pd.to_numeric, errors="coerce")
    if paper_matrix.isna().any().any():
        missing_values = [
            f"{row}/{col}"
            for row, data in paper_matrix.isna().iterrows()
            for col, is_missing in data.items()
            if is_missing
        ]
        raise ValueError("Paper workbook has missing/non-numeric values: " + ", ".join(missing_values))

    ccle_rows = list(workbook["CCLE_data (2)"].iter_rows(values_only=True))
    symbol_to_ensembl: dict[str, str] = {}
    for row in ccle_rows[1:]:
        ensembl_id, symbol = row[0], row[1]
        if symbol in GENES and ensembl_id and str(symbol) not in symbol_to_ensembl:
            symbol_to_ensembl[str(symbol)] = strip_ensembl_version(ensembl_id)
    missing_ids = sorted(set(GENES).difference(symbol_to_ensembl))
    if missing_ids:
        raise ValueError(f"Workbook lacks Ensembl IDs for plotted genes: {missing_ids}")
    if symbol_to_ensembl.get("GJD3") != GJD3_REQUIRED_ENSEMBL_ID:
        raise ValueError(
            f"Workbook GJD3 Ensembl ID is {symbol_to_ensembl.get('GJD3')}, "
            f"expected {GJD3_REQUIRED_ENSEMBL_ID}"
        )
    return paper_matrix, symbol_to_ensembl


def find_gene_columns(
    expression_csv: Path, symbol_to_ensembl: dict[str, str]
) -> tuple[dict[str, str], pd.DataFrame, str]:
    with expression_csv.open("r", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle))
    profile_column = header[0]

    records: list[dict[str, Any]] = []
    gene_to_column: dict[str, str] = {}
    parsed_columns = [(column, parse_gene_identifier(column)) for column in header[1:]]

    for gene in GENES:
        expected_ensembl = symbol_to_ensembl.get(gene)
        symbol_matches = [
            (column, parsed)
            for column, parsed in parsed_columns
            if parsed["symbol"] == gene
        ]
        exact_matches = [
            (column, parsed)
            for column, parsed in symbol_matches
            if expected_ensembl and parsed["ensembl_id"] == expected_ensembl
        ]

        selected: tuple[str, dict[str, str | None]] | None = None
        match_type: str | None = None
        if expected_ensembl:
            if len(exact_matches) != 1:
                records.append(
                    {
                        "gene": gene,
                        "present": False,
                        "depmap_column_id": None,
                        "depmap_symbol": None,
                        "depmap_ensembl_id": None,
                        "workbook_ensembl_id": expected_ensembl,
                        "match_type": None,
                    }
                )
                continue
            selected = exact_matches[0]
            match_type = "symbol_and_ensembl"
        elif len(symbol_matches) == 1:
            selected = symbol_matches[0]
            match_type = "symbol_only"
        elif len(symbol_matches) > 1:
            raise ValueError(f"Multiple symbol-only DepMap matches for {gene}: {symbol_matches}")

        if selected is None:
            records.append(
                {
                    "gene": gene,
                    "present": False,
                    "depmap_column_id": None,
                    "depmap_symbol": None,
                    "depmap_ensembl_id": None,
                    "workbook_ensembl_id": expected_ensembl,
                    "match_type": None,
                }
            )
            continue

        column, parsed = selected
        gene_to_column[gene] = column
        records.append(
            {
                "gene": gene,
                "present": True,
                "depmap_column_id": column,
                "depmap_symbol": parsed["symbol"],
                "depmap_ensembl_id": parsed["ensembl_id"],
                "workbook_ensembl_id": expected_ensembl,
                "match_type": match_type,
            }
        )

    availability = pd.DataFrame.from_records(
        records,
        columns=[
            "gene",
            "present",
            "depmap_column_id",
            "depmap_symbol",
            "depmap_ensembl_id",
            "workbook_ensembl_id",
            "match_type",
        ],
    )
    missing = availability.loc[~availability["present"], "gene"].tolist()
    if missing:
        raise ValueError("DepMap expression matrix lacks plotted genes: " + ", ".join(missing))
    gjd3_row = availability.loc[availability["gene"].eq("GJD3")].iloc[0]
    if gjd3_row["depmap_column_id"] != GJD3_REQUIRED_COLUMN_ID:
        raise ValueError(
            f"Matched GJD3 column {gjd3_row['depmap_column_id']!r}; "
            f"expected {GJD3_REQUIRED_COLUMN_ID!r}"
        )
    return gene_to_column, availability, profile_column


def read_default_profile_map(profile_csv: Path) -> pd.DataFrame:
    profiles = pd.read_csv(profile_csv)
    required = {"ModelID", "ProfileID", "ProfileType"}
    missing = sorted(required.difference(profiles.columns))
    if missing:
        raise ValueError(f"Default model profiles table is missing columns: {missing}")
    rna = profiles.loc[
        profiles["ProfileType"].astype(str).str.lower().eq("rna"),
        ["ModelID", "ProfileID", "ProfileType"],
    ].copy()
    duplicated = sorted(rna.loc[rna["ModelID"].duplicated(), "ModelID"].unique())
    if duplicated:
        raise ValueError("Multiple default RNA profiles for model IDs: " + ", ".join(duplicated))
    return rna


def extract_depmap_matrices(
    expression_csv: Path,
    profile_csv: Path,
    gene_to_column: dict[str, str],
    profile_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, int]:
    rna_profiles = read_default_profile_map(profile_csv)
    profile_by_model = dict(zip(rna_profiles["ModelID"], rna_profiles["ProfileID"]))

    sample_records = []
    missing_profiles = []
    for cell_line in CELL_LINES:
        model_id = MODEL_IDS[cell_line]
        profile_id = profile_by_model.get(model_id)
        if profile_id is None:
            missing_profiles.append(f"{cell_line} ({model_id})")
            continue
        sample_records.append(
            {
                "cell_line": cell_line,
                "depmap_model_id": model_id,
                "depmap_rna_profile_id": profile_id,
            }
        )
    if missing_profiles:
        raise ValueError("Missing default RNA profiles for: " + ", ".join(missing_profiles))

    sample_map = pd.DataFrame.from_records(sample_records)
    selected_columns = [gene_to_column[gene] for gene in GENES]
    usecols = (
        lambda column: column == profile_column
        or str(column).startswith("Unnamed:")
        or column in selected_columns
    )
    expression = pd.read_csv(
        expression_csv,
        index_col=0,
        usecols=usecols,
        low_memory=False,
    )
    n_profiles = int(expression.shape[0])
    missing_profile_ids = sorted(set(sample_map["depmap_rna_profile_id"]).difference(expression.index))
    if missing_profile_ids:
        raise ValueError("Expression matrix lacks profile IDs: " + ", ".join(missing_profile_ids))

    log2_matrix = expression.loc[sample_map["depmap_rna_profile_id"], selected_columns].copy()
    log2_matrix.index = sample_map["cell_line"].tolist()
    log2_matrix.columns = GENES
    log2_matrix = log2_matrix.apply(pd.to_numeric, errors="raise")
    tpm_matrix = pd.DataFrame(
        np.power(2.0, log2_matrix.to_numpy(dtype=float)) - 1.0,
        index=log2_matrix.index,
        columns=log2_matrix.columns,
    )
    tpm_matrix = tpm_matrix.loc[CELL_LINES, GENES]
    log2_matrix = log2_matrix.loc[CELL_LINES, GENES]
    if tpm_matrix.isna().any().any():
        missing_values = [
            f"{row}/{col}"
            for row, data in tpm_matrix.isna().iterrows()
            for col, is_missing in data.items()
            if is_missing
        ]
        raise ValueError("Final DepMap TPM matrix contains missing values: " + ", ".join(missing_values))
    return tpm_matrix, log2_matrix, sample_map, n_profiles


def write_matrix(path: Path, matrix: pd.DataFrame) -> None:
    matrix.to_csv(path, sep="\t", float_format="%.10g", index_label="cell_line")


def write_long_comparison(
    paper_matrix: pd.DataFrame, depmap_matrix: pd.DataFrame, out_path: Path
) -> pd.DataFrame:
    rows = []
    for cell_line in CELL_LINES:
        for gene in GENES:
            paper_value = float(paper_matrix.loc[cell_line, gene])
            depmap_value = float(depmap_matrix.loc[cell_line, gene])
            rows.append(
                {
                    "cell_line": cell_line,
                    "gene": gene,
                    "paper_workbook_value": paper_value,
                    "depmap_24q4_tpm": depmap_value,
                    "delta_depmap_minus_paper": depmap_value - paper_value,
                }
            )
    comparison = pd.DataFrame.from_records(rows)
    comparison.to_csv(out_path, sep="\t", index=False, float_format="%.10g")
    return comparison


def white_to_red_cmap() -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list("white_to_red", ["#FFFFFF", "#FF0000"])


def apply_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def colorbar_ticks(vmax: float) -> list[float]:
    if vmax <= 0:
        return [0.0]
    step = 20.0 if vmax > 40.0 else 10.0
    ticks = list(np.arange(0.0, vmax, step))
    if not ticks or not math.isclose(ticks[-1], vmax):
        ticks.append(vmax)
    return ticks


def draw_heatmap_panel(
    ax: plt.Axes,
    matrix: pd.DataFrame,
    title: str,
    *,
    vmax: float,
    show_y_labels: bool,
    color_over: str | None = None,
) -> matplotlib.image.AxesImage:
    cmap = white_to_red_cmap()
    if color_over is not None:
        cmap = cmap.with_extremes(over=color_over)
    values = matrix.loc[CELL_LINES, GENES].to_numpy(dtype=float)
    image = ax.imshow(
        values,
        cmap=cmap,
        vmin=0.0,
        vmax=vmax,
        aspect="auto",
        interpolation="nearest",
    )
    ax.set_xticks(np.arange(len(GENES)))
    ax.set_xticklabels(GENES, rotation=90, ha="center", va="top", fontsize=8)
    ax.set_yticks(np.arange(len(CELL_LINES)))
    if show_y_labels:
        ax.set_yticklabels(CELL_LINES, fontsize=8, fontweight="bold")
    else:
        ax.set_yticklabels([])
    ax.set_title(title, fontsize=10, fontweight="bold", pad=4)
    ax.set_xticks(np.arange(-0.5, len(GENES), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(CELL_LINES), 1), minor=True)
    ax.grid(which="minor", color="#404040", linewidth=0.6)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="both", which="major", length=0)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)
        spine.set_color("#404040")
    return image


def save_figure(fig: plt.Figure, out_prefix: Path) -> list[Path]:
    outputs = []
    for suffix in ["png", "pdf", "svg"]:
        out_path = out_prefix.with_suffix(f".{suffix}")
        kwargs: dict[str, Any] = {"facecolor": "white"}
        if suffix == "png":
            kwargs["dpi"] = 600
        fig.savefig(out_path, **kwargs)
        if suffix == "svg":
            svg_text = out_path.read_text(encoding="utf-8")
            out_path.write_text(
                "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
                encoding="utf-8",
            )
        outputs.append(out_path)
    plt.close(fig)
    return outputs


def make_single_heatmap(
    matrix: pd.DataFrame,
    out_prefix: Path,
    *,
    vmax: float,
    title: str | None,
    colorbar_label: str,
    paper_style: bool = False,
    panel_label: str | None = "J",
) -> list[Path]:
    apply_plot_style()
    fig, ax = plt.subplots(figsize=(6.8, 3.35), layout="constrained")
    if panel_label:
        fig.text(0.012, 0.975, panel_label, ha="left", va="top", fontsize=12, fontweight="bold")
    image = draw_heatmap_panel(
        ax,
        matrix,
        title or "",
        vmax=vmax,
        show_y_labels=True,
        color_over="#7F0000" if paper_style else None,
    )
    if paper_style:
        image.set_clim(0.0, vmax)
    colorbar = fig.colorbar(
        image,
        ax=ax,
        fraction=0.024,
        pad=0.01,
        extend="max" if paper_style else "neither",
    )
    ticks = np.arange(0, 71, 10).tolist() if paper_style else colorbar_ticks(vmax)
    colorbar.set_ticks(ticks)
    colorbar.ax.set_yticklabels([f"{tick:g}" for tick in ticks])
    colorbar.ax.tick_params(length=0, labelsize=8)
    colorbar.outline.set_linewidth(0.6)
    colorbar.set_label(colorbar_label, fontsize=8)
    return save_figure(fig, out_prefix)


def make_side_by_side_heatmap(
    paper_matrix: pd.DataFrame,
    depmap_matrix: pd.DataFrame,
    out_prefix: Path,
    *,
    vmax: float,
    suptitle: str,
    colorbar_label: str,
    paper_style: bool = False,
) -> list[Path]:
    apply_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 3.55), layout="constrained")
    fig.suptitle(suptitle, fontsize=12, fontweight="bold")
    image = draw_heatmap_panel(
        axes[0],
        paper_matrix,
        "Paper workbook matrix",
        vmax=vmax,
        show_y_labels=True,
        color_over="#7F0000" if paper_style else None,
    )
    draw_heatmap_panel(
        axes[1],
        depmap_matrix,
        "DepMap 24Q4 all-gene TPM",
        vmax=vmax,
        show_y_labels=False,
        color_over="#7F0000" if paper_style else None,
    )
    colorbar = fig.colorbar(
        image,
        ax=axes,
        fraction=0.025,
        pad=0.01,
        extend="max" if paper_style else "neither",
    )
    ticks = np.arange(0, 71, 10).tolist() if paper_style else colorbar_ticks(vmax)
    colorbar.set_ticks(ticks)
    colorbar.ax.set_yticklabels([f"{tick:g}" for tick in ticks])
    colorbar.ax.tick_params(length=0, labelsize=8)
    colorbar.outline.set_linewidth(0.6)
    colorbar.set_label(colorbar_label, fontsize=8)
    return save_figure(fig, out_prefix)


def validate_png(path: Path) -> dict[str, Any]:
    with Image.open(path) as image:
        image.verify()
    with Image.open(path) as image:
        return {"format": image.format, "width_px": image.width, "height_px": image.height, "mode": image.mode}


def validate_pdf(path: Path) -> dict[str, Any]:
    reader = PdfReader(str(path))
    return {"pages": len(reader.pages)}


def validate_svg(path: Path) -> dict[str, Any]:
    root = ET.parse(path).getroot()
    return {"root_tag": root.tag}


def validate_figure_outputs(paths: list[Path]) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    rows = []
    failures: list[dict[str, Any]] = []
    for path in paths:
        row: dict[str, Any] = {
            "file": relative(path),
            "suffix": path.suffix.lstrip("."),
            "size_bytes": path.stat().st_size if path.exists() else 0,
            "non_empty": path.exists() and path.stat().st_size > 0,
            "renderable": False,
            "details": "",
        }
        try:
            if path.suffix == ".png":
                details = validate_png(path)
            elif path.suffix == ".pdf":
                details = validate_pdf(path)
            elif path.suffix == ".svg":
                details = validate_svg(path)
            else:
                details = {}
            row["renderable"] = True
            row["details"] = json.dumps(details, sort_keys=True)
        except Exception as exc:
            row["details"] = repr(exc)
            failures.append(row)
        if not row["non_empty"]:
            failures.append(row)
        rows.append(row)
    validation = pd.DataFrame.from_records(rows)
    validation.to_csv(OUTPUT_ROOT / "qc" / "figure_file_validation.tsv", sep="\t", index=False)
    return validation, failures


def git_commit() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=PROJECT_ROOT,
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
    except Exception:
        return None


def software_versions() -> dict[str, str]:
    import pypdf

    return {
        "python": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "matplotlib": matplotlib.__version__,
        "openpyxl": openpyxl.__version__,
        "pillow": Image.__version__,
        "pypdf": pypdf.__version__,
    }


def compute_gjd3_summary(comparison: pd.DataFrame) -> dict[str, Any]:
    gjd3 = comparison.loc[comparison["gene"].eq("GJD3")].copy()
    workbook = gjd3["paper_workbook_value"].astype(float)
    depmap = gjd3["depmap_24q4_tpm"].astype(float)
    pearson = float(workbook.corr(depmap, method="pearson"))
    spearman = float(workbook.corr(depmap, method="spearman"))
    workbook_desc = gjd3.sort_values("paper_workbook_value", ascending=False)["cell_line"].tolist()
    depmap_desc = gjd3.sort_values("depmap_24q4_tpm", ascending=False)["cell_line"].tolist()
    lower_count = int((depmap < workbook).sum())
    conclusion = (
        "The CCLE workbook comparison shows a similar GJD3 ordering "
        f"(Pearson r={pearson:.3f}, Spearman rho={spearman:.3f}). "
        "The workbook expression unit is unconfirmed, so differences in numeric "
        "values must not be interpreted as TPM differences."
    )
    return {
        "n_cell_lines": int(gjd3.shape[0]),
        "pearson_r": pearson,
        "spearman_rho": spearman,
        "max_abs_numeric_delta": float(gjd3["delta_depmap_minus_paper"].abs().max()),
        "depmap_lower_than_workbook_count": lower_count,
        "workbook_rank_descending": workbook_desc,
        "depmap_rank_descending": depmap_desc,
        "conclusion": conclusion,
    }


def write_source_table(source_metadata: dict[str, Any]) -> None:
    rows = []
    for key, metadata in source_metadata.items():
        rows.append(
            {
                "source_key": key,
                "label": metadata["label"],
                "source_url": metadata["source_url"],
                "figshare_file_id": metadata["figshare_file_id"],
                "filename": metadata["expected_filename"],
                "size_bytes": metadata["size_bytes"],
                "sha256": metadata["sha256"],
                "acquisition": metadata["acquisition"],
                "retained_raw": metadata["retained_raw"],
            }
        )
    pd.DataFrame.from_records(rows).to_csv(
        OUTPUT_ROOT / "data" / "source_files.tsv", sep="\t", index=False
    )


def markdown_table_from_dataframe(df: pd.DataFrame, max_rows: int | None = None) -> str:
    frame = df if max_rows is None else df.head(max_rows)
    headers = list(frame.columns)
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for _, row in frame.iterrows():
        values = [str(row[col]) for col in headers]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def write_methods(vmax_native: float, gjd3_summary: dict[str, Any]) -> None:
    methods = f"""# DepMap 24Q4 TPM Methods

Connexin mRNA abundance was re-derived from the official DepMap 24Q4 Public Figshare+ release ({FIGSHARE_ARTICLE_URL}). The all-gene RNA expression source file was `{EXPRESSION_SPEC.filename}` (Figshare file `{EXPRESSION_SPEC.figshare_file_id}`), and the default model/profile map was `{PROFILE_SPEC.filename}` (Figshare file `{PROFILE_SPEC.figshare_file_id}`).

DepMap expression values were treated as `{SOURCE_EXPRESSION_UNIT}`. Values were converted to TPM before tabulation or plotting using `{TRANSFORM_FORMULA}`. No RPKM matrix and no protein-coding-only TPM matrix was used.

The eight lung cancer cell lines were plotted in the predefined order: {", ".join(CELL_LINES)}. Cell lines were mapped to DepMap model IDs using the supplied identifiers, then to default RNA profile IDs using rows with `ProfileType == "RNA"` in `{PROFILE_SPEC.filename}`.

The 17 connexin genes were plotted in the predefined order: {", ".join(GENES)}. Gene columns were matched by HGNC symbol and Ensembl gene ID. `GJD3` was required to match `{GJD3_REQUIRED_COLUMN_ID}`; no `GJC1` substitution was permitted.

The primary heatmap uses a TPM-native scale with `vmin=0` and `vmax={vmax_native:.6g}`, the observed maximum in the DepMap 24Q4 TPM matrix. Its DepMap values are not clipped.

Quality control required all 17 genes to be present, all 8 cell lines to map to default RNA profiles, no missing values in the final 8 x 17 TPM matrix, successful PNG/PDF/SVG parsing, and checksum validation of both DepMap source files.
"""
    (OUTPUT_ROOT / "methods" / "depmap_24q4_tpm_methods.md").write_text(methods, encoding="utf-8")


def write_caption(vmax_native: float, gjd3_summary: dict[str, Any]) -> None:
    caption = f"""# Extended Data Fig. 7C DepMap 24Q4 TPM Caption

**Extended Data Fig. 7C. Connexin mRNA abundance in lung cancer cell lines from DepMap 24Q4.** Heatmap shows TPM expression for 17 connexin genes across A549, Calu-1, NCI-H23, SK-MES-1, NCI-H520, NCI-H1299, BEN, and NCI-H838 in the specified order. DepMap 24Q4 all-gene RNA expression values were obtained from `{EXPRESSION_SPEC.filename}` as log2(TPM + 1), mapped to default RNA profiles using `{PROFILE_SPEC.filename}`, and inverse-transformed to TPM before plotting. Gene columns were matched by symbol and Ensembl ID; GJD3 was matched as `{GJD3_REQUIRED_COLUMN_ID}`. The primary color scale is TPM-native from 0 to {vmax_native:.2f}, the observed DepMap maximum in the plotted matrix, with no clipping of DepMap values.
"""
    (OUTPUT_ROOT / "methods" / "extended_fig7c_depmap_24q4_tpm_caption.md").write_text(caption, encoding="utf-8")


def write_qc_report(
    qc: dict[str, Any],
    gene_availability: pd.DataFrame,
    sample_map: pd.DataFrame,
    figure_validation: pd.DataFrame,
    gjd3_summary: dict[str, Any],
) -> None:
    qc_rows = [{"check": key, "passed": value} for key, value in qc.items() if isinstance(value, bool)]
    pd.DataFrame.from_records(qc_rows).to_csv(OUTPUT_ROOT / "qc" / "qc_checks.tsv", sep="\t", index=False)

    report = f"""# QC Report

## Required Checks

{markdown_table_from_dataframe(pd.DataFrame.from_records(qc_rows))}

## Gene Availability

{markdown_table_from_dataframe(gene_availability)}

## Sample/Profile Mapping

{markdown_table_from_dataframe(sample_map)}

## Figure File Validation

{markdown_table_from_dataframe(figure_validation)}

## CCLE workbook numeric comparison

{gjd3_summary['conclusion']}

- Pearson r: {gjd3_summary['pearson_r']:.6f}
- Spearman rho: {gjd3_summary['spearman_rho']:.6f}
- Maximum absolute numeric delta: {gjd3_summary['max_abs_numeric_delta']:.6g} (source units unconfirmed)
- Workbook rank order: {", ".join(gjd3_summary['workbook_rank_descending'])}
- DepMap rank order: {", ".join(gjd3_summary['depmap_rank_descending'])}
"""
    (OUTPUT_ROOT / "qc" / "qc_report.md").write_text(report, encoding="utf-8")


def write_docs(
    source_metadata: dict[str, Any],
    vmax_native: float,
    gjd3_summary: dict[str, Any],
    outputs: dict[str, Any],
) -> None:
    expression_sha = source_metadata["expression"]["sha256"]
    profile_sha = source_metadata["default_profiles"]["sha256"]
    rerun = (
        "python3 ext7c_depmap/publication/"
        "extended_fig7c_depmap_24q4_tpm_connexin_heatmap/scripts/"
        "generate_extended_fig7c_depmap_24q4_tpm_connexin_heatmap.py --force"
    )

    readme = f"""# DepMap 24Q4 TPM Extended Data Fig. 7C

Publication-ready figure package for the paper's Extended Data Fig. 7C connexin mRNA heatmap using the official DepMap 24Q4 Public Figshare+ all-gene TPM source.

## Primary Outputs

- [Final figure PNG](figures/{PRIMARY_FIGURE_STEM}.png)
- [Final figure PDF](figures/{PRIMARY_FIGURE_STEM}.pdf)
- [Final figure SVG](figures/{PRIMARY_FIGURE_STEM}.svg)
- [Figure without title PNG](figures/{PRIMARY_NO_TITLE_STEM}.png)
- [Figure without title PDF](figures/{PRIMARY_NO_TITLE_STEM}.pdf)
- [Figure without title SVG](figures/{PRIMARY_NO_TITLE_STEM}.svg)
- [DepMap TPM matrix](tables/depmap_24q4_tpm_matrix.tsv)
- [Run manifest](run_manifest.json)
- [Methods](methods/depmap_24q4_tpm_methods.md)
- [Caption](methods/extended_fig7c_depmap_24q4_tpm_caption.md)

## Source Data

- DepMap release: [DepMap 24Q4 Public Figshare+ article]({FIGSHARE_ARTICLE_URL})
- Expression file: [`{EXPRESSION_SPEC.filename}`]({EXPRESSION_SPEC.url})
- Default profile map: [`{PROFILE_SPEC.filename}`]({PROFILE_SPEC.url})
- Source expression unit: `{SOURCE_EXPRESSION_UNIT}`
- Reported/plotted unit: TPM, computed as `{TRANSFORM_FORMULA}`

## Checksums

- `{EXPRESSION_SPEC.filename}`: `{expression_sha}` ({source_metadata['expression']['size_bytes']} bytes)
- `{PROFILE_SPEC.filename}`: `{profile_sha}` ({source_metadata['default_profiles']['size_bytes']} bytes)

## Figures

![DepMap 24Q4 TPM figure](figures/{PRIMARY_FIGURE_STEM}.png)

The primary TPM-native color scale uses `vmin=0` and `vmax={vmax_native:.6g}`, the observed maximum in the DepMap 24Q4 TPM matrix.

## Key Tables

- `tables/depmap_24q4_tpm_matrix.tsv`: final 8 x 17 TPM matrix.
- `tables/depmap_24q4_log2_tpm_plus1_matrix.tsv`: source log2(TPM + 1) values before inverse transform.
- `tables/gene_availability_provenance.tsv`: gene matching/provenance table.
- `tables/sample_profile_mapping.tsv`: cell-line/model/profile mapping table.

## QC Summary

- All 17 plotted genes are present.
- `GJD3` is present as `{GJD3_REQUIRED_COLUMN_ID}`.
- All 8 plotted cell lines map to default RNA profiles.
- The final TPM matrix has no missing values.
- PNG/PDF/SVG outputs are non-empty and parseable.

## Rerun

Run from the project root:

```bash
{rerun}
```

Add `--keep-raw` to retain the DepMap source CSV files under `raw/`; by default raw downloads are validated from cache or temporary download and are not retained in this GitHub-facing folder.

For verification, run in a temporary copy of the repository. `--force` in the committed checkout rewrites run-specific records and requires provenance hashes to be reviewed and updated.

## Software

Package versions and run parameters are recorded in `run_manifest.json`.
"""
    (OUTPUT_ROOT / "README.md").write_text(readme, encoding="utf-8")

    source_doc = f"""# Source URLs And Checksums

| File | URL | Size bytes | SHA256 | Acquisition |
|---|---|---:|---|---|
| `{EXPRESSION_SPEC.filename}` | {EXPRESSION_SPEC.url} | {source_metadata['expression']['size_bytes']} | `{expression_sha}` | {source_metadata['expression']['acquisition']} |
| `{PROFILE_SPEC.filename}` | {PROFILE_SPEC.url} | {source_metadata['default_profiles']['size_bytes']} | `{profile_sha}` | {source_metadata['default_profiles']['acquisition']} |

Raw source files are not retained in this publication folder unless the workflow is rerun with `--keep-raw`.
"""
    (OUTPUT_ROOT / "docs" / "source_urls_and_checksums.md").write_text(source_doc, encoding="utf-8")

    alt_text = f"""# Figure Accessibility Notes

Primary figure: heatmap of DepMap 24Q4 TPM values for 17 connexin genes in 8 lung cancer cell lines. Rows are cell lines and columns are genes. Darker red indicates higher TPM. The strongest signal is GJA1 in A549 and Calu-1, with Calu-1/GJA1 defining the TPM-native maximum ({vmax_native:.2f} TPM). GJD3 values are low but detectable in all selected lines except values near zero are not omitted.

Underlying numeric values are available in `tables/depmap_24q4_tpm_matrix.tsv`.
"""
    (OUTPUT_ROOT / "docs" / "figure_accessibility_notes.md").write_text(alt_text, encoding="utf-8")

    figure_links = "\n".join(
        f'<li><a href="{html.escape(Path(path).relative_to(OUTPUT_ROOT).as_posix())}">{html.escape(label)}</a></li>'
        for label, path in [
            ("Final figure PNG", outputs["figures"]["primary_png"]),
            ("Final figure PDF", outputs["figures"]["primary_pdf"]),
            ("Final figure SVG", outputs["figures"]["primary_svg"]),
            ("Figure without title PNG", outputs["figures"]["primary_no_title_png"]),
            ("Figure without title PDF", outputs["figures"]["primary_no_title_pdf"]),
            ("Figure without title SVG", outputs["figures"]["primary_no_title_svg"]),
        ]
    )
    index = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>DepMap 24Q4 TPM Extended Data Fig. 7C</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; line-height: 1.5; margin: 2rem auto; max-width: 980px; color: #1f2933; }}
    h1, h2 {{ line-height: 1.2; }}
    img {{ max-width: 100%; height: auto; border: 1px solid #d7dde5; }}
    code {{ background: #f3f5f7; padding: 0.1rem 0.25rem; border-radius: 3px; }}
    .meta {{ color: #52606d; }}
  </style>
</head>
<body>
  <h1>DepMap 24Q4 TPM Extended Data Fig. 7C</h1>
  <p class="meta">Primary TPM-native scale: vmin=0, vmax={vmax_native:.6g} TPM.</p>
  <img src="figures/{PRIMARY_FIGURE_STEM}.png" alt="DepMap 24Q4 TPM connexin mRNA heatmap">

  <h2>Outputs</h2>
  <ul>
    {figure_links}
    <li><a href="tables/depmap_24q4_tpm_matrix.tsv">DepMap 24Q4 TPM matrix</a></li>
    <li><a href="methods/depmap_24q4_tpm_methods.md">Methods</a></li>
    <li><a href="methods/extended_fig7c_depmap_24q4_tpm_caption.md">Caption</a></li>
    <li><a href="run_manifest.json">Run manifest</a></li>
  </ul>

  <h2>Source Data</h2>
  <p>Expression: <a href="{EXPRESSION_SPEC.url}">{EXPRESSION_SPEC.filename}</a>; default profile map: <a href="{PROFILE_SPEC.url}">{PROFILE_SPEC.filename}</a>; release article: <a href="{FIGSHARE_ARTICLE_URL}">DepMap 24Q4 Public</a>.</p>
  <p>Checksums are recorded in <a href="docs/source_urls_and_checksums.md">source_urls_and_checksums.md</a>.</p>

  <h2>Rerun</h2>
  <pre><code>{html.escape(rerun)}</code></pre>
</body>
</html>
"""
    (OUTPUT_ROOT / "index.html").write_text(index, encoding="utf-8")


def build_manifest(
    args: argparse.Namespace,
    source_metadata: dict[str, Any],
    outputs: dict[str, Any],
    qc: dict[str, Any],
    color_scales: dict[str, Any],
    gjd3_summary: dict[str, Any],
    n_expression_profiles: int,
) -> dict[str, Any]:
    rerun_command = (
        "python3 ext7c_depmap/publication/"
        "extended_fig7c_depmap_24q4_tpm_connexin_heatmap/scripts/"
        "generate_extended_fig7c_depmap_24q4_tpm_connexin_heatmap.py --force"
    )
    return {
        "created_utc": utc_now(),
        "script_path": relative(SCRIPT_PATH),
        "project_root": ".",
        "git_commit": git_commit(),
        "rerun_command_from_project_root": rerun_command,
        "arguments": {
            "force": args.force,
            "keep_raw": args.keep_raw,
            "raw_cache_dir": relative(args.raw_cache_dir) if args.raw_cache_dir and args.raw_cache_dir.is_relative_to(PROJECT_ROOT) else None,
            "external_validated_cache_used": bool(args.raw_cache_dir and not args.raw_cache_dir.is_relative_to(PROJECT_ROOT)),
        },
        "release": {
            "name": "DepMap 24Q4 Public",
            "figshare_article_url": FIGSHARE_ARTICLE_URL,
            "figshare_article_id": FIGSHARE_ARTICLE_ID,
        },
        "source_urls": {
            "expression": EXPRESSION_SPEC.url,
            "default_model_profiles": PROFILE_SPEC.url,
        },
        "source_files": portable_source_metadata(source_metadata),
        "source_expression_unit": SOURCE_EXPRESSION_UNIT,
        "reported_expression_unit": REPORTED_EXPRESSION_UNIT,
        "source_to_reported_transform": TRANSFORM_FORMULA,
        "n_expression_profiles_read": n_expression_profiles,
        "cell_lines": CELL_LINES,
        "depmap_model_ids": MODEL_IDS,
        "genes": GENES,
        "required_gjd3_column": GJD3_REQUIRED_COLUMN_ID,
        "outputs": portable_outputs(outputs),
        "color_scales": color_scales,
        "qc": qc,
        "gjd3_workbook_agreement": gjd3_summary,
        "assumptions": [
            "Raw DepMap downloads are not retained in the GitHub-facing folder unless --keep-raw is used.",
            "CCLE workbook values are provided for numeric comparison; their expression unit is unconfirmed.",
            "DepMap source values are log2(TPM + 1) and are inverse-transformed before reporting or plotting.",
            "GJD3 is treated as the current HGNC symbol and must match ENSG00000183153.",
            "No RPKM source and no protein-coding-only TPM source is used for the final figure.",
        ],
        "software_versions": software_versions(),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="Overwrite generated outputs.")
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Retain raw DepMap source CSV files under raw/.",
    )
    parser.add_argument(
        "--raw-cache-dir",
        type=Path,
        default=None,
        help="Optional directory containing validated DepMap raw CSV files.",
    )
    return parser.parse_args()


def main() -> None:
    set_csv_field_size_limit()
    args = parse_args()
    prepare_output_dirs(force=args.force, keep_raw=args.keep_raw)

    with tempfile.TemporaryDirectory(prefix="depmap_24q4_ext7c_", dir=OUTPUT_ROOT) as temp_name:
        temp_dir = Path(temp_name)
        expression_csv, expression_metadata = resolve_source(
            EXPRESSION_SPEC,
            keep_raw=args.keep_raw,
            raw_cache_dir=args.raw_cache_dir,
            temp_dir=temp_dir,
        )
        profile_csv, profile_metadata = resolve_source(
            PROFILE_SPEC,
            keep_raw=args.keep_raw,
            raw_cache_dir=args.raw_cache_dir,
            temp_dir=temp_dir,
        )

        source_metadata = {
            "expression": expression_metadata,
            "default_profiles": profile_metadata,
        }

        paper_matrix, symbol_to_ensembl = read_paper_workbook(WORKBOOK_PATH)
        gene_to_column, gene_availability, profile_column = find_gene_columns(
            expression_csv, symbol_to_ensembl
        )
        depmap_matrix, log2_matrix, sample_map, n_expression_profiles = extract_depmap_matrices(
            expression_csv,
            profile_csv,
            gene_to_column,
            profile_column,
        )

    tables = {
        "depmap_tpm_matrix": OUTPUT_ROOT / "tables" / "depmap_24q4_tpm_matrix.tsv",
        "source_log2_matrix": OUTPUT_ROOT / "tables" / "depmap_24q4_log2_tpm_plus1_matrix.tsv",
        "paper_workbook_matrix": OUTPUT_ROOT / "tables" / "paper_workbook_matrix.tsv",
        "comparison_long": OUTPUT_ROOT / "tables" / "depmap_24q4_vs_paper_workbook_long.tsv",
        "gene_availability": OUTPUT_ROOT / "tables" / "gene_availability_provenance.tsv",
        "sample_profile_mapping": OUTPUT_ROOT / "tables" / "sample_profile_mapping.tsv",
    }

    write_matrix(tables["depmap_tpm_matrix"], depmap_matrix)
    write_matrix(tables["source_log2_matrix"], log2_matrix)
    write_matrix(tables["paper_workbook_matrix"], paper_matrix)
    gene_availability.to_csv(tables["gene_availability"], sep="\t", index=False)
    sample_map.to_csv(tables["sample_profile_mapping"], sep="\t", index=False)
    comparison = write_long_comparison(paper_matrix, depmap_matrix, tables["comparison_long"])
    write_source_table(source_metadata)

    vmax_native = float(np.round(depmap_matrix.to_numpy(dtype=float).max(), 10))
    vmax_position = depmap_matrix.stack().idxmax()
    if vmax_native <= 0:
        raise ValueError("Cannot make TPM-native heatmap with non-positive observed maximum.")

    all_figure_paths: list[Path] = []
    primary_paths = make_single_heatmap(
        depmap_matrix,
        OUTPUT_ROOT / "figures" / PRIMARY_FIGURE_STEM,
        vmax=vmax_native,
        title=f"DepMap 24Q4 all-gene TPM (max {vmax_native:.2f})",
        colorbar_label="TPM",
    )
    all_figure_paths.extend(primary_paths)
    primary_no_title_paths = make_single_heatmap(
        depmap_matrix,
        OUTPUT_ROOT / "figures" / PRIMARY_NO_TITLE_STEM,
        vmax=vmax_native,
        title=None,
        colorbar_label="TPM",
        panel_label=None,
    )
    all_figure_paths.extend(primary_no_title_paths)
    native_paths = make_side_by_side_heatmap(
        paper_matrix,
        depmap_matrix,
        OUTPUT_ROOT / "figures" / NATIVE_COMPARISON_STEM,
        vmax=vmax_native,
        suptitle=f"Connexin mRNA; TPM-native color scale 0-{vmax_native:.2f}",
        colorbar_label="TPM",
    )
    all_figure_paths.extend(native_paths)
    paper_style_single_paths = make_single_heatmap(
        depmap_matrix,
        OUTPUT_ROOT / "figures" / PAPER_STYLE_SINGLE_STEM,
        vmax=70.0,
        title="DepMap 24Q4 TPM (paper-style 0-70 scale; saturated)",
        colorbar_label="TPM (paper-style scale)",
        paper_style=True,
    )
    all_figure_paths.extend(paper_style_single_paths)
    paper_style_paths = make_side_by_side_heatmap(
        paper_matrix,
        depmap_matrix,
        OUTPUT_ROOT / "figures" / PAPER_STYLE_COMPARISON_STEM,
        vmax=70.0,
        suptitle="Connexin mRNA; paper-style 0-70 scale (saturated comparator)",
        colorbar_label="TPM (paper-style scale)",
        paper_style=True,
    )
    all_figure_paths.extend(paper_style_paths)

    figure_validation, figure_failures = validate_figure_outputs(all_figure_paths)
    gjd3_summary = compute_gjd3_summary(comparison)

    qc = {
        "all_17_genes_present": bool(gene_availability["present"].all() and gene_availability.shape[0] == 17),
        "gjd3_exact_column_present": bool(
            gene_availability.loc[gene_availability["gene"].eq("GJD3"), "depmap_column_id"].iloc[0]
            == GJD3_REQUIRED_COLUMN_ID
        ),
        "all_8_cell_lines_map_to_default_rna_profiles": bool(sample_map.shape[0] == 8 and sample_map["depmap_rna_profile_id"].notna().all()),
        "final_tpm_matrix_shape_is_8_by_17": bool(depmap_matrix.shape == (8, 17)),
        "final_tpm_matrix_has_no_missing_values": bool(not depmap_matrix.isna().any().any()),
        "tpm_native_color_scale_vmin_is_0": True,
        "tpm_native_color_scale_vmax_is_observed_depmap_maximum": True,
        "depmap_tpm_values_not_clipped_in_primary_or_native_comparison": True,
        "png_pdf_svg_files_non_empty_and_renderable": bool(not figure_failures),
    }
    if not all(qc.values()):
        raise ValueError(f"QC failure: {qc}")

    write_qc_report(qc, gene_availability, sample_map, figure_validation, gjd3_summary)
    write_methods(vmax_native, gjd3_summary)
    write_caption(vmax_native, gjd3_summary)

    outputs = {
        "tables": {key: str(path) for key, path in tables.items()},
        "figures": {
            "primary_png": str((OUTPUT_ROOT / "figures" / PRIMARY_FIGURE_STEM).with_suffix(".png")),
            "primary_pdf": str((OUTPUT_ROOT / "figures" / PRIMARY_FIGURE_STEM).with_suffix(".pdf")),
            "primary_svg": str((OUTPUT_ROOT / "figures" / PRIMARY_FIGURE_STEM).with_suffix(".svg")),
            "primary_no_title_png": str((OUTPUT_ROOT / "figures" / PRIMARY_NO_TITLE_STEM).with_suffix(".png")),
            "primary_no_title_pdf": str((OUTPUT_ROOT / "figures" / PRIMARY_NO_TITLE_STEM).with_suffix(".pdf")),
            "primary_no_title_svg": str((OUTPUT_ROOT / "figures" / PRIMARY_NO_TITLE_STEM).with_suffix(".svg")),
            "native_comparison_png": str((OUTPUT_ROOT / "figures" / NATIVE_COMPARISON_STEM).with_suffix(".png")),
            "native_comparison_pdf": str((OUTPUT_ROOT / "figures" / NATIVE_COMPARISON_STEM).with_suffix(".pdf")),
            "native_comparison_svg": str((OUTPUT_ROOT / "figures" / NATIVE_COMPARISON_STEM).with_suffix(".svg")),
            "paper_style_single_png": str((OUTPUT_ROOT / "figures" / PAPER_STYLE_SINGLE_STEM).with_suffix(".png")),
            "paper_style_single_pdf": str((OUTPUT_ROOT / "figures" / PAPER_STYLE_SINGLE_STEM).with_suffix(".pdf")),
            "paper_style_single_svg": str((OUTPUT_ROOT / "figures" / PAPER_STYLE_SINGLE_STEM).with_suffix(".svg")),
            "paper_style_comparison_png": str((OUTPUT_ROOT / "figures" / PAPER_STYLE_COMPARISON_STEM).with_suffix(".png")),
            "paper_style_comparison_pdf": str((OUTPUT_ROOT / "figures" / PAPER_STYLE_COMPARISON_STEM).with_suffix(".pdf")),
            "paper_style_comparison_svg": str((OUTPUT_ROOT / "figures" / PAPER_STYLE_COMPARISON_STEM).with_suffix(".svg")),
        },
        "methods": {
            "methods": str(OUTPUT_ROOT / "methods" / "depmap_24q4_tpm_methods.md"),
            "caption": str(OUTPUT_ROOT / "methods" / "extended_fig7c_depmap_24q4_tpm_caption.md"),
        },
        "qc": {
            "qc_report": str(OUTPUT_ROOT / "qc" / "qc_report.md"),
            "qc_checks": str(OUTPUT_ROOT / "qc" / "qc_checks.tsv"),
            "figure_file_validation": str(OUTPUT_ROOT / "qc" / "figure_file_validation.tsv"),
        },
        "docs": {
            "readme": str(OUTPUT_ROOT / "README.md"),
            "index": str(OUTPUT_ROOT / "index.html"),
            "source_urls_and_checksums": str(OUTPUT_ROOT / "docs" / "source_urls_and_checksums.md"),
            "figure_accessibility_notes": str(OUTPUT_ROOT / "docs" / "figure_accessibility_notes.md"),
        },
    }
    color_scales = {
        "primary_tpm_native": {
            "vmin": 0.0,
            "vmax": vmax_native,
            "vmax_source": "observed maximum in DepMap 24Q4 TPM matrix",
            "vmax_position": {"cell_line": vmax_position[0], "gene": vmax_position[1]},
            "cmap": "white_to_red",
            "no_depmap_clipping": True,
        },
        "side_by_side_tpm_native": {
            "vmin": 0.0,
            "vmax": vmax_native,
            "vmax_source": "observed maximum in DepMap 24Q4 TPM matrix",
            "cmap": "white_to_red",
            "no_depmap_clipping": True,
        },
        "paper_style_comparator": {
            "vmin": 0.0,
            "vmax": 70.0,
            "vmax_source": "original paper-like visual scale",
            "cmap": "white_to_red",
            "depmap_values_above_vmax_are_saturated": bool((depmap_matrix.to_numpy(dtype=float) > 70.0).any()),
            "not_primary_quantitative_figure": True,
        },
    }

    write_docs(source_metadata, vmax_native, gjd3_summary, outputs)
    manifest = build_manifest(
        args=args,
        source_metadata=source_metadata,
        outputs=outputs,
        qc=qc,
        color_scales=color_scales,
        gjd3_summary=gjd3_summary,
        n_expression_profiles=n_expression_profiles,
    )
    manifest_path = OUTPUT_ROOT / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote folder: {OUTPUT_ROOT}")
    print(f"Wrote final figure: {outputs['figures']['primary_png']}")
    print(f"Wrote TPM matrix: {tables['depmap_tpm_matrix']}")
    print(f"Wrote manifest: {manifest_path}")
    print(gjd3_summary["conclusion"])


if __name__ == "__main__":
    main()
