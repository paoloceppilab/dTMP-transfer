#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import socket
import sys
import time
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import chain
from pathlib import Path
from typing import Iterable, Iterator, TextIO


GENETIC_CODE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

KNOWN_ONCOGENIC_AA = {
    12: set("ACDRSV"),
    13: set("ACDRSV"),
    61: set("HKLR"),
}

CODON_START_DEFAULTS = {
    12: 21,
    13: 24,
}


@dataclass(frozen=True)
class Sample:
    specimen_id: str
    cohort: str
    analysis_name: str


@dataclass
class AlleleCall:
    observed_by_pos: dict[int, str]
    insertions_by_anchor: dict[int, list[str]]
    window_observed: str
    event_notation: str
    has_window_indel: bool
    has_window_substitution: bool
    affects_codon12: bool
    affects_codon13: bool
    affects_codon61: bool
    codon12_nt: str
    codon12_aa: str
    codon12_effect: str
    codon13_nt: str
    codon13_aa: str
    codon13_effect: str
    known_oncogenic_residue_flag: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize Kras codon-centered allele spectra from existing CRISPResso "
            "prefix-amplicon allele tables. Intended for codon 12/13 windows that are "
            "inside the 67 bp forward-prefix Kras amplicon."
        )
    )
    parser.add_argument("--analysis-manifest", required=True)
    parser.add_argument("--amplicons", required=True)
    parser.add_argument("--editing-summary", required=True)
    parser.add_argument("--crispresso-dir", required=True)
    parser.add_argument("--locus", default="Kras")
    parser.add_argument("--codon-number", type=int, default=13)
    parser.add_argument(
        "--codon-start",
        type=int,
        default=None,
        help="1-based amplicon start for the target codon. Defaults are known for Kras codons 12 and 13.",
    )
    parser.add_argument("--flank-bp", type=int, default=20)
    parser.add_argument("--min-classifiable-reads", type=int, default=10000)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def open_text_maybe_gzip(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def iter_fastq(path: Path) -> Iterator[tuple[str, str, str, str]]:
    with open_text_maybe_gzip(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip().upper()
            plus = handle.readline()
            qual = handle.readline().strip()
            if not qual:
                raise ValueError(f"Malformed FASTQ record in {path}")
            yield header, seq, plus, qual


def translate(codon: str) -> str:
    if len(codon) != 3 or any(base not in "ACGT" for base in codon):
        return "NA"
    return GENETIC_CODE.get(codon, "NA")


def codon_effect(codon_number: int, ref_codon: str, observed: str) -> tuple[str, str]:
    aa = translate(observed)
    ref_aa = translate(ref_codon)
    if aa == "NA":
        return aa, "uncallable_indel_or_ambiguous"
    if observed == ref_codon:
        return aa, "exact_wt_codon"
    if aa == ref_aa:
        return aa, "synonymous_noncanonical_codon"
    if aa == "*":
        return aa, "stop_gain"
    if aa in KNOWN_ONCOGENIC_AA.get(codon_number, set()):
        return aa, "known_oncogenic_residue_change"
    return aa, "missense_unknown_significance"


def find_samples(manifest_rows: list[dict[str, str]], locus: str) -> list[Sample]:
    samples = []
    for row in manifest_rows:
        if row.get("locus") != locus:
            continue
        if row.get("analysis_include", "yes").lower() != "yes":
            continue
        samples.append(
            Sample(
                specimen_id=row["specimen_id"],
                cohort=row["cohort"],
                analysis_name=row["analysis_name"],
            )
        )
    if not samples:
        raise ValueError(f"No included samples found for locus={locus!r}")
    return samples


def pick_count_column(fieldnames: list[str]) -> str:
    candidates = ["#Reads", "Reads", "read_count", "n_reads", "Read_Count", "Count"]
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    raise ValueError(f"Could not find read-count column in CRISPResso table header: {fieldnames}")


def find_crispresso_allele_table(crispresso_dir: Path, analysis_name: str) -> Path:
    sample_dirs = [p for p in crispresso_dir.rglob(f"*{analysis_name}*") if p.is_dir()]
    search_roots = sample_dirs or [crispresso_dir]
    candidates: list[Path] = []
    for root in search_roots:
        candidates.extend(root.rglob("Alleles_frequency_table*.txt"))
        candidates.extend(root.rglob("Alleles_frequency_table*.tsv"))
        candidates.extend(root.rglob("Alleles_frequency_table*.txt.gz"))
        candidates.extend(root.rglob("Alleles_frequency_table*.zip"))
    candidates = sorted(set(candidates))
    if not candidates:
        raise FileNotFoundError(
            f"No Alleles_frequency_table file found for {analysis_name} under {crispresso_dir}"
        )
    scored = []
    for path in candidates:
        name = path.name
        score = 0
        if path.parent.name == analysis_name or analysis_name in str(path.parent):
            score += 10
        if name == "Alleles_frequency_table.txt":
            score += 5
        if path.suffix == ".zip":
            score += 2
        scored.append((score, len(str(path)), path))
    return sorted(scored, key=lambda item: (-item[0], item[1]))[0][2]


def iter_crispresso_rows(path: Path) -> Iterator[dict[str, str]]:
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as zf:
            members = [
                name
                for name in zf.namelist()
                if Path(name).name.startswith("Alleles_frequency_table")
                and Path(name).suffix in {".txt", ".tsv"}
            ]
            if not members:
                raise ValueError(f"No allele table member found in {path}")
            with zf.open(sorted(members)[0]) as raw:
                text = (line.decode() for line in raw)
                yield from csv.DictReader(text, delimiter="\t")
        return
    with open_text_maybe_gzip(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def reference_name(row: dict[str, str]) -> str:
    for key in ["Reference_Name", "Reference name", "ref_name", "Reference"]:
        if key in row:
            return row[key]
    return ""


def aligned_sequence_fields(row: dict[str, str]) -> tuple[str, str]:
    aligned_keys = ["Aligned_Sequence", "Aligned Sequence", "aligned_sequence"]
    reference_keys = ["Reference_Sequence", "Reference Sequence", "reference_sequence"]
    aligned = next((row[key] for key in aligned_keys if key in row), "")
    reference = next((row[key] for key in reference_keys if key in row), "")
    if not aligned or not reference:
        raise ValueError(f"Missing aligned/reference sequence columns in row keys: {list(row)}")
    return aligned.upper(), reference.upper()


def build_call(
    aligned_seq: str,
    aligned_ref: str,
    wt_seq: str,
    hdr_seq: str,
    window_start: int,
    window_end: int,
    codon12_start: int,
    codon13_start: int,
) -> AlleleCall:
    observed_by_pos: dict[int, str] = {}
    insertions_by_anchor: dict[int, list[str]] = defaultdict(list)
    ref_pos = 0
    for obs, ref in zip(aligned_seq, aligned_ref):
        if ref != "-":
            ref_pos += 1
            observed_by_pos[ref_pos] = obs
        elif obs != "-":
            insertions_by_anchor[ref_pos].append(obs)

    events: list[str] = []
    has_window_indel = False
    has_window_substitution = False
    for pos in range(window_start, window_end + 1):
        obs = observed_by_pos.get(pos, "N")
        ref = wt_seq[pos - 1]
        if obs == "-":
            has_window_indel = True
            events.append(f"{pos}del{ref}")
        elif obs not in "ACGT":
            has_window_substitution = True
            events.append(f"{pos}{ref}>{obs}")
        elif obs != ref:
            has_window_substitution = True
            events.append(f"{pos}{ref}>{obs}")
    for anchor, inserted in sorted(insertions_by_anchor.items()):
        if window_start - 1 <= anchor <= window_end:
            has_window_indel = True
            events.append(f"{anchor}_{anchor + 1}ins{''.join(inserted)}")

    def codon_has_indel(start: int) -> bool:
        codon_positions = range(start, start + 3)
        has_deletion = any(observed_by_pos.get(pos, "N") == "-" for pos in codon_positions)
        has_insertion = any(start - 1 <= anchor <= start + 2 for anchor in insertions_by_anchor)
        return has_deletion or has_insertion

    def codon_nt(start: int) -> str:
        if codon_has_indel(start):
            return "INDEL"
        bases = [observed_by_pos.get(pos, "N") for pos in range(start, start + 3)]
        return "".join(bases)

    codon12_ref = wt_seq[codon12_start - 1 : codon12_start + 2]
    codon13_ref = wt_seq[codon13_start - 1 : codon13_start + 2]
    codon12_nt = codon_nt(codon12_start)
    codon13_nt = codon_nt(codon13_start)
    codon12_aa, codon12_eff = codon_effect(12, codon12_ref, codon12_nt)
    codon13_aa, codon13_eff = codon_effect(13, codon13_ref, codon13_nt)

    affects_codon12 = any(
        observed_by_pos.get(pos, "N") != wt_seq[pos - 1]
        for pos in range(codon12_start, codon12_start + 3)
    ) or any(codon12_start - 1 <= anchor <= codon12_start + 2 for anchor in insertions_by_anchor)
    affects_codon13 = any(
        observed_by_pos.get(pos, "N") != wt_seq[pos - 1]
        for pos in range(codon13_start, codon13_start + 3)
    ) or any(codon13_start - 1 <= anchor <= codon13_start + 2 for anchor in insertions_by_anchor)

    observed_window_parts: list[str] = []
    for pos in range(window_start, window_end + 1):
        for inserted in insertions_by_anchor.get(pos - 1, []):
            observed_window_parts.append(f"+{inserted}")
        observed_window_parts.append(observed_by_pos.get(pos, "N"))
    for inserted in insertions_by_anchor.get(window_end, []):
        observed_window_parts.append(f"+{inserted}")

    flags = []
    if codon12_eff == "known_oncogenic_residue_change":
        flags.append(f"codon12_{codon12_aa}")
    if codon13_eff == "known_oncogenic_residue_change":
        flags.append(f"codon13_{codon13_aa}")
    known_oncogenic_residue_flag = ";".join(flags) if flags else "none_detected"

    return AlleleCall(
        observed_by_pos=observed_by_pos,
        insertions_by_anchor=dict(insertions_by_anchor),
        window_observed="".join(observed_window_parts),
        event_notation=";".join(events) if events else "none",
        has_window_indel=has_window_indel,
        has_window_substitution=has_window_substitution,
        affects_codon12=affects_codon12,
        affects_codon13=affects_codon13,
        affects_codon61=False,
        codon12_nt=codon12_nt,
        codon12_aa=codon12_aa,
        codon12_effect=codon12_eff,
        codon13_nt=codon13_nt,
        codon13_aa=codon13_aa,
        codon13_effect=codon13_eff,
        known_oncogenic_residue_flag=known_oncogenic_residue_flag,
    )


def strict_hdr_barcode(call: AlleleCall, wt_seq: str, hdr_seq: str, prefix_len: int = 67) -> bool:
    diff_positions = [
        pos
        for pos in range(1, min(prefix_len, len(wt_seq), len(hdr_seq)) + 1)
        if wt_seq[pos - 1] != hdr_seq[pos - 1]
    ]
    return all(call.observed_by_pos.get(pos) == hdr_seq[pos - 1] for pos in diff_positions)


def window_exact(call: AlleleCall, seq: str, window_start: int, window_end: int) -> bool:
    if any(window_start - 1 <= anchor <= window_end for anchor in call.insertions_by_anchor):
        return False
    return all(call.observed_by_pos.get(pos) == seq[pos - 1] for pos in range(window_start, window_end + 1))


def window_class(call: AlleleCall, wt_seq: str, hdr_seq: str, window_start: int, window_end: int) -> str:
    if window_exact(call, wt_seq, window_start, window_end):
        return "window_exact_wt"
    if window_exact(call, hdr_seq, window_start, window_end):
        return "window_exact_hdr"
    if call.has_window_indel and call.has_window_substitution:
        return "window_complex_indel_plus_substitution"
    if call.has_window_indel:
        return "window_indel"
    if call.has_window_substitution:
        return "window_substitution_only"
    return "window_other_uncallable"


def load_existing_hdr_by_sample(path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(path)
    return {row["specimen_id"]: row for row in rows if row.get("kras_analysis_name")}


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def pct(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return "NA"
    return f"{100.0 * numerator / denominator:.6f}"


def main() -> None:
    args = parse_args()
    start_time = time.time()
    outdir = Path(args.outdir)
    summary_dir = outdir / "summary"
    per_allele_dir = outdir / "per_allele"
    figures_dir = outdir / "figures"
    metadata_dir = outdir / "metadata"
    for directory in [summary_dir, per_allele_dir, figures_dir, metadata_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    analysis_manifest = Path(args.analysis_manifest)
    amplicons = Path(args.amplicons)
    editing_summary = Path(args.editing_summary)
    crispresso_dir = Path(args.crispresso_dir)

    amp_rows = read_tsv(amplicons)
    amp = next((row for row in amp_rows if row["locus"] == args.locus), None)
    if amp is None:
        raise ValueError(f"Locus {args.locus!r} not found in {amplicons}")
    wt_seq = amp["wt_amplicon_seq"].upper()
    hdr_seq = amp["expected_hdr_amplicon_seq"].upper()
    if not hdr_seq or hdr_seq == "NA":
        hdr_seq = wt_seq

    codon_start = args.codon_start or CODON_START_DEFAULTS.get(args.codon_number)
    if codon_start is None:
        raise ValueError("--codon-start is required unless codon-number is 12 or 13")
    codon_end = codon_start + 2
    window_start = max(1, codon_start - args.flank_bp)
    window_end = min(len(wt_seq), codon_end + args.flank_bp)

    codon12_start = CODON_START_DEFAULTS[12]
    codon13_start = CODON_START_DEFAULTS[13]
    if window_end > 67:
        raise ValueError(
            f"Requested window {window_start}-{window_end} extends beyond the existing 67 bp prefix. "
            "Use raw reads or a longer amplicon analysis for this window."
        )

    samples = find_samples(read_tsv(analysis_manifest), args.locus)
    existing_hdr = load_existing_hdr_by_sample(editing_summary)

    per_allele_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    consequence_counter: dict[tuple[str, str, str, str], int] = Counter()
    figure_rows: list[dict[str, object]] = []
    allele_table_paths: dict[str, str] = {}

    for sample in samples:
        table_path = find_crispresso_allele_table(crispresso_dir, sample.analysis_name)
        allele_table_paths[sample.analysis_name] = str(table_path)
        rows_iter = iter_crispresso_rows(table_path)
        first_row = next(rows_iter, None)
        if first_row is None:
            raise ValueError(f"Empty CRISPResso allele table for {sample.analysis_name}: {table_path}")
        fieldnames = list(first_row)
        count_col = pick_count_column(fieldnames)

        counters = Counter()
        counters["total_allele_table_reads"] = 0
        allele_rank_rows: list[dict[str, object]] = []

        for row in chain([first_row], rows_iter):
            try:
                read_count = int(float(row[count_col]))
            except ValueError as exc:
                raise ValueError(f"Invalid read count {row[count_col]!r} in {table_path}") from exc
            if read_count <= 0:
                continue
            aligned_seq, aligned_ref = aligned_sequence_fields(row)
            call = build_call(
                aligned_seq=aligned_seq,
                aligned_ref=aligned_ref,
                wt_seq=wt_seq,
                hdr_seq=hdr_seq,
                window_start=window_start,
                window_end=window_end,
                codon12_start=codon12_start,
                codon13_start=codon13_start,
            )
            counters["total_allele_table_reads"] += read_count
            counters["window_classifiable_reads"] += read_count
            wclass = window_class(call, wt_seq, hdr_seq, window_start, window_end)
            counters[f"{wclass}_reads"] += read_count
            counters[f"codon13_{call.codon13_effect}_reads"] += read_count
            counters[f"codon12_{call.codon12_effect}_reads"] += read_count
            if strict_hdr_barcode(call, wt_seq, hdr_seq):
                counters["strict_hdr_barcode_reads"] += read_count
            if call.known_oncogenic_residue_flag != "none_detected":
                counters["known_oncogenic_residue_change_reads"] += read_count
            if call.has_window_indel:
                counters["window_any_indel_reads"] += read_count
            if call.affects_codon13:
                counters["affects_codon13_reads"] += read_count

            consequence_counter[
                (
                    sample.specimen_id,
                    sample.cohort,
                    call.codon13_effect,
                    call.known_oncogenic_residue_flag,
                )
            ] += read_count

            allele_rank_rows.append(
                {
                    "specimen_id": sample.specimen_id,
                    "cohort": sample.cohort,
                    "analysis_name": sample.analysis_name,
                    "reference_name": reference_name(row),
                    "read_count": read_count,
                    "window_start": window_start,
                    "window_end": window_end,
                    "target_codon_number": args.codon_number,
                    "target_codon_start": codon_start,
                    "target_codon_end": codon_end,
                    "window_class": wclass,
                    "strict_hdr_barcode": "yes" if strict_hdr_barcode(call, wt_seq, hdr_seq) else "no",
                    "codon12_nt": call.codon12_nt,
                    "codon12_aa": call.codon12_aa,
                    "codon12_effect": call.codon12_effect,
                    "codon13_nt": call.codon13_nt,
                    "codon13_aa": call.codon13_aa,
                    "codon13_effect": call.codon13_effect,
                    "affects_codon12": "yes" if call.affects_codon12 else "no",
                    "affects_codon13": "yes" if call.affects_codon13 else "no",
                    "affects_codon61": "no_not_assayed",
                    "known_oncogenic_residue_flag": call.known_oncogenic_residue_flag,
                    "event_notation": call.event_notation,
                    "observed_window": call.window_observed,
                }
            )

        classifiable = counters["window_classifiable_reads"]
        hdr_row = existing_hdr.get(sample.specimen_id, {})
        qc_flags = []
        if classifiable < args.min_classifiable_reads:
            qc_flags.append("low_classifiable_reads")
        if not table_path.exists():
            qc_flags.append("missing_allele_table")
        qc_flag = "pass" if not qc_flags else ";".join(qc_flags)

        summary_row = {
            "specimen_id": sample.specimen_id,
            "cohort": sample.cohort,
            "analysis_name": sample.analysis_name,
            "target_codon_number": args.codon_number,
            "target_codon_start": codon_start,
            "target_codon_end": codon_end,
            "window_start": window_start,
            "window_end": window_end,
            "window_sequence_wt": wt_seq[window_start - 1 : window_end],
            "window_sequence_hdr": hdr_seq[window_start - 1 : window_end],
            "existing_kras_hdr_exact_pct_prefix67": hdr_row.get("kras_hdr_exact_pct", "NA"),
            "existing_kras_hdr_exact_reads_prefix67": hdr_row.get("kras_hdr_exact_reads", "NA"),
            "existing_kras_informative_reads_prefix67": hdr_row.get("kras_informative_reads", "NA"),
            "total_allele_table_reads": counters["total_allele_table_reads"],
            "window_classifiable_reads": classifiable,
            "window_exact_wt_reads": counters["window_exact_wt_reads"],
            "window_exact_wt_pct": pct(counters["window_exact_wt_reads"], classifiable),
            "window_exact_hdr_reads": counters["window_exact_hdr_reads"],
            "window_exact_hdr_pct": pct(counters["window_exact_hdr_reads"], classifiable),
            "strict_hdr_barcode_reads_observed": counters["strict_hdr_barcode_reads"],
            "strict_hdr_barcode_pct_observed": pct(counters["strict_hdr_barcode_reads"], classifiable),
            "codon13_exact_wt_reads": counters["codon13_exact_wt_codon_reads"],
            "codon13_exact_wt_pct": pct(counters["codon13_exact_wt_codon_reads"], classifiable),
            "codon13_synonymous_noncanonical_reads": counters[
                "codon13_synonymous_noncanonical_codon_reads"
            ],
            "codon13_synonymous_noncanonical_pct": pct(
                counters["codon13_synonymous_noncanonical_codon_reads"], classifiable
            ),
            "codon13_known_oncogenic_residue_reads": counters[
                "codon13_known_oncogenic_residue_change_reads"
            ],
            "codon13_known_oncogenic_residue_pct": pct(
                counters["codon13_known_oncogenic_residue_change_reads"], classifiable
            ),
            "codon13_missense_unknown_reads": counters[
                "codon13_missense_unknown_significance_reads"
            ],
            "codon13_missense_unknown_pct": pct(
                counters["codon13_missense_unknown_significance_reads"], classifiable
            ),
            "codon13_stop_gain_reads": counters["codon13_stop_gain_reads"],
            "codon13_stop_gain_pct": pct(counters["codon13_stop_gain_reads"], classifiable),
            "codon13_uncallable_reads": counters["codon13_uncallable_indel_or_ambiguous_reads"],
            "codon13_uncallable_pct": pct(
                counters["codon13_uncallable_indel_or_ambiguous_reads"], classifiable
            ),
            "window_any_indel_reads": counters["window_any_indel_reads"],
            "window_any_indel_pct": pct(counters["window_any_indel_reads"], classifiable),
            "window_substitution_only_reads": counters["window_substitution_only_reads"],
            "window_substitution_only_pct": pct(counters["window_substitution_only_reads"], classifiable),
            "window_complex_indel_plus_substitution_reads": counters[
                "window_complex_indel_plus_substitution_reads"
            ],
            "window_complex_indel_plus_substitution_pct": pct(
                counters["window_complex_indel_plus_substitution_reads"], classifiable
            ),
            "affects_codon13_reads": counters["affects_codon13_reads"],
            "affects_codon13_pct": pct(counters["affects_codon13_reads"], classifiable),
            "known_oncogenic_residue_change_reads": counters["known_oncogenic_residue_change_reads"],
            "known_oncogenic_residue_change_pct": pct(
                counters["known_oncogenic_residue_change_reads"], classifiable
            ),
            "codon61_status": "not_assayed_by_67bp_prefix_or_134bp_amplicon",
            "qc_flag": qc_flag,
            "crispresso_allele_table": str(table_path),
        }
        summary_rows.append(summary_row)

        for category, value in [
            ("exact_wt_window", counters["window_exact_wt_reads"]),
            ("exact_hdr_window", counters["window_exact_hdr_reads"]),
            ("codon13_known_oncogenic", counters["codon13_known_oncogenic_residue_change_reads"]),
            ("codon13_synonymous_noncanonical", counters["codon13_synonymous_noncanonical_codon_reads"]),
            ("window_indel_any", counters["window_any_indel_reads"]),
            ("window_substitution_only", counters["window_substitution_only_reads"]),
        ]:
            figure_rows.append(
                {
                    "specimen_id": sample.specimen_id,
                    "cohort": sample.cohort,
                    "category": category,
                    "read_count": value,
                    "denominator": classifiable,
                    "pct": pct(value, classifiable),
                }
            )

        per_allele_rows.extend(
            sorted(allele_rank_rows, key=lambda row: int(row["read_count"]), reverse=True)
        )

    sample_summary_name = f"kras_codon{args.codon_number}_window_sample_summary.tsv"
    allele_name = f"kras_codon{args.codon_number}_window_allele_counts.tsv"
    consequence_name = f"kras_codon{args.codon_number}_window_consequence_summary.tsv"
    figure_name = f"kras_codon{args.codon_number}_window_figure_ready.tsv"

    summary_fields = list(summary_rows[0])
    write_tsv(summary_dir / sample_summary_name, summary_rows, summary_fields)

    allele_fields = list(per_allele_rows[0]) if per_allele_rows else [
        "specimen_id",
        "cohort",
        "analysis_name",
        "read_count",
    ]
    write_tsv(per_allele_dir / allele_name, per_allele_rows, allele_fields)

    consequence_rows = [
        {
            "specimen_id": specimen_id,
            "cohort": cohort,
            "codon13_effect": codon13_eff,
            "known_oncogenic_residue_flag": onc_flag,
            "read_count": count,
        }
        for (specimen_id, cohort, codon13_eff, onc_flag), count in sorted(consequence_counter.items())
    ]
    write_tsv(
        summary_dir / consequence_name,
        consequence_rows,
        [
            "specimen_id",
            "cohort",
            "codon13_effect",
            "known_oncogenic_residue_flag",
            "read_count",
        ],
    )
    write_tsv(figures_dir / figure_name, figure_rows, list(figure_rows[0]))

    run_info = {
        "script": str(Path(__file__).resolve()),
        "argv": sys.argv,
        "start_time_utc": datetime.fromtimestamp(start_time, timezone.utc).isoformat(),
        "end_time_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(time.time() - start_time, 3),
        "python": sys.version,
        "platform": platform.platform(),
        "hostname": socket.gethostname(),
        "parameters": {
            "locus": args.locus,
            "codon_number": args.codon_number,
            "codon_start": codon_start,
            "flank_bp": args.flank_bp,
            "window_start": window_start,
            "window_end": window_end,
            "min_classifiable_reads": args.min_classifiable_reads,
        },
        "inputs": {
            "analysis_manifest": str(analysis_manifest),
            "analysis_manifest_sha256": sha256(analysis_manifest),
            "amplicons": str(amplicons),
            "amplicons_sha256": sha256(amplicons),
            "editing_summary": str(editing_summary),
            "editing_summary_sha256": sha256(editing_summary),
            "crispresso_dir": str(crispresso_dir),
            "allele_table_paths": allele_table_paths,
        },
        "outputs": {
            "sample_summary": str(summary_dir / sample_summary_name),
            "allele_counts": str(per_allele_dir / allele_name),
            "consequence_summary": str(summary_dir / consequence_name),
            "figure_ready": str(figures_dir / figure_name),
        },
    }
    (metadata_dir / "run_info.json").write_text(json.dumps(run_info, indent=2) + "\n")
    (metadata_dir / "run_manifest.yaml").write_text(
        "\n".join(
            [
                f"script: {run_info['script']}",
                f"created_utc: {run_info['end_time_utc']}",
                f"locus: {args.locus}",
                f"target_codon_number: {args.codon_number}",
                f"target_codon_start_1based: {codon_start}",
                f"target_codon_end_1based: {codon_end}",
                f"flank_bp: {args.flank_bp}",
                f"window_start_1based: {window_start}",
                f"window_end_1based: {window_end}",
                f"analysis_manifest: {analysis_manifest}",
                f"amplicons: {amplicons}",
                f"editing_summary: {editing_summary}",
                f"crispresso_dir: {crispresso_dir}",
                "method: existing CRISPResso prefix67 allele-table parsing; valid because requested window is within bases 1-67",
                "codon61_status: not_assayed_by_this_amplicon",
                "",
            ]
        )
    )
    (outdir / "README.md").write_text(
        "\n".join(
            [
                f"# Kras Codon {args.codon_number} Window Spectrum",
                "",
                f"Target codon: Kras codon {args.codon_number}, amplicon bases {codon_start}-{codon_end}.",
                f"Window: bases {window_start}-{window_end} ({args.flank_bp} bp flanking the codon, bounded by the amplicon).",
                "",
                "This analysis parses existing Kras CRISPResso prefix67 allele tables. The approach is appropriate here because the requested codon-centered window is fully contained within the 67 bp forward-prefix amplicon. It should not be used for the base-91 cut-site question.",
                "",
                "Primary HDR KrasG12D percentages should still be taken from the existing `editing_summary.tsv` prefix67 exact barcode metric. Codon13 calls are reported separately as WT, synonymous noncanonical codons, known oncogenic-residue changes, unknown missense, stop/uncallable, and window indel/substitution classes.",
                "",
            ]
        )
    )


if __name__ == "__main__":
    main()
