#!/usr/bin/env python3
"""Analyze Kras cut-site allele spectra from paired amplicon FASTQs.

This script is intentionally dependency-light and deterministic. It uses raw paired
FASTQs, selects reverse-primer-side reads that cover the Kras cut site, aligns only
observed read sequence to the amplicon suffix, and writes auditable per-pair,
per-allele, per-sample, and consequence summaries.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import os
import platform
import socket
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import zip_longest
from pathlib import Path
from typing import Iterable

BASES = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass(frozen=True)
class Candidate:
    source_read: str
    raw_reverse_seq: str
    raw_reverse_qual: str
    primer_offset: int
    trimmed_len: int
    nominal_start_1based: int
    nominal_window_minq: int | None


@dataclass
class AlignmentResult:
    status: str
    reason: str
    score: int | None
    ref_start_1based: int | None
    ref_end_1based: int | None
    aligned_ref: str
    aligned_read: str
    aligned_qual: list[int | None]
    ref_positions: list[int | None]
    n_inserted: int
    n_deleted: int
    n_substituted: int
    n_inserted_cut: int
    n_deleted_cut: int
    n_substituted_cut: int
    net_indel_length_cut: int
    full_event_key: str
    allele_key: str
    cutsite_event_type: str
    frameshift_status: str
    consequence: str
    predicted_lof: str
    codon12_status: str
    codon13_status: str
    codon61_status: str
    known_oncogenic_status: str
    cut_window_min_baseq: int | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Kras cut-site allele-spectrum analysis from raw paired FASTQs."
    )
    parser.add_argument("--analysis-manifest", required=True)
    parser.add_argument("--amplicons", required=True)
    parser.add_argument("--editing-summary", required=True)
    parser.add_argument("--locus", default="Kras")
    parser.add_argument("--cut-site", type=int, default=91)
    parser.add_argument("--cut-window-start", type=int, default=86)
    parser.add_argument("--cut-window-end", type=int, default=96)
    parser.add_argument("--context-window-start", type=int, default=81)
    parser.add_argument("--context-window-end", type=int, default=101)
    parser.add_argument("--primer-kmer-len", type=int, default=10)
    parser.add_argument("--primer-search-window", type=int, default=5)
    parser.add_argument("--prefix-length", type=int, default=67)
    parser.add_argument("--min-window-baseq", type=int, default=20)
    parser.add_argument("--min-classifiable-pairs", type=int, default=10000)
    parser.add_argument("--max-alignment-score", type=int, default=35)
    parser.add_argument("--gap-cost", type=int, default=2)
    parser.add_argument("--mismatch-cost", type=int, default=1)
    parser.add_argument("--alignment-upstream-padding", type=int, default=35)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--max-pairs-per-sample", type=int, default=None,
                        help="Debug/testing only: cap pairs processed per sample.")
    parser.add_argument("--log-level", default="INFO")
    return parser.parse_args()


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S%z",
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def revcomp(seq: str) -> str:
    return seq.translate(BASES)[::-1].upper()


def normalize_read_id(header: str) -> str:
    rid = header.split()[0].lstrip("@")
    if rid.endswith("/1") or rid.endswith("/2"):
        rid = rid[:-2]
    return rid


def iter_fastq(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        while True:
            header = handle.readline().rstrip("\n")
            if not header:
                break
            seq = handle.readline().strip().upper()
            handle.readline()
            qual = handle.readline().strip()
            yield normalize_read_id(header), seq, qual


def find_primer_start(seq: str, primer_kmer: str, search_window: int) -> int | None:
    max_start = min(search_window, max(0, len(seq) - len(primer_kmer)))
    for start in range(max_start + 1):
        if seq[start:start + len(primer_kmer)] == primer_kmer:
            return start
    return None


def phred_min(qual: str) -> int | None:
    if not qual:
        return None
    return min(ord(ch) - 33 for ch in qual)


def reverse_candidates(
    read_label: str,
    seq: str,
    qual: str,
    reverse_primer_kmer: str,
    search_window: int,
    wt_len: int,
    cut_window_start: int,
    cut_window_end: int,
) -> list[Candidate]:
    candidates: list[Candidate] = []
    direct_start = find_primer_start(seq, reverse_primer_kmer, search_window)
    orientations: list[tuple[str, str, int]] = []
    if direct_start is not None:
        orientations.append((seq, qual, direct_start))

    rc_seq = revcomp(seq)
    rc_qual = qual[::-1]
    rc_start = find_primer_start(rc_seq, reverse_primer_kmer, search_window)
    if rc_start is not None:
        orientations.append((rc_seq, rc_qual, rc_start))

    for oriented_seq, oriented_qual, start in orientations:
        trimmed_seq = oriented_seq[start:]
        trimmed_qual = oriented_qual[start:]
        if not trimmed_seq:
            continue
        nominal_start = wt_len - len(trimmed_seq) + 1
        if nominal_start > cut_window_start:
            continue
        s0 = max(0, cut_window_start - nominal_start)
        s1 = min(len(trimmed_qual), cut_window_end - nominal_start + 1)
        minq = phred_min(trimmed_qual[::-1][s0:s1]) if s1 > s0 else None
        candidates.append(
            Candidate(
                source_read=read_label,
                raw_reverse_seq=trimmed_seq,
                raw_reverse_qual=trimmed_qual,
                primer_offset=start,
                trimmed_len=len(trimmed_seq),
                nominal_start_1based=nominal_start,
                nominal_window_minq=minq,
            )
        )
    return candidates


def forward_oriented_sequences(
    seq: str,
    qual: str,
    forward_primer_kmer: str,
    search_window: int,
) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    direct_start = find_primer_start(seq, forward_primer_kmer, search_window)
    if direct_start is not None:
        out.append((seq[direct_start:], qual[direct_start:]))
    rc_seq = revcomp(seq)
    rc_qual = qual[::-1]
    rc_start = find_primer_start(rc_seq, forward_primer_kmer, search_window)
    if rc_start is not None:
        out.append((rc_seq[rc_start:], rc_qual[rc_start:]))
    return out


def classify_hdr_barcode(
    forward_sequences: list[tuple[str, str]],
    informative_positions: list[int],
    wt_key: str,
    hdr_key: str,
) -> str:
    if not forward_sequences:
        return "no_forward_primer_read"
    statuses = []
    min_len = max(informative_positions) + 1
    for seq, _qual in forward_sequences:
        if len(seq) < min_len:
            statuses.append("short")
            continue
        key = "".join(seq[pos] for pos in informative_positions)
        if key == hdr_key:
            return "hdr_exact"
        if key == wt_key:
            statuses.append("wt_exact")
        else:
            statuses.append("other_pattern")
    if "wt_exact" in statuses:
        return "wt_exact"
    if "other_pattern" in statuses:
        return "other_pattern"
    return "short"


def primer_side_status(
    seq: str,
    fwd_kmer: str,
    rev_kmer: str,
    search_window: int,
) -> str:
    if find_primer_start(seq, fwd_kmer, search_window) is not None:
        return "forward_primer"
    if find_primer_start(seq, rev_kmer, search_window) is not None:
        return "reverse_primer"
    rc_seq = revcomp(seq)
    if find_primer_start(rc_seq, fwd_kmer, search_window) is not None:
        return "forward_primer_rc"
    if find_primer_start(rc_seq, rev_kmer, search_window) is not None:
        return "reverse_primer_rc"
    return "unmatched"


def semiglobal_align_suffix(
    read_seq: str,
    read_qual: str,
    ref: str,
    cut_window: tuple[int, int],
    context_window: tuple[int, int],
    codon_map: dict[str, tuple[int | None, int | None]],
    coding_frame_start_1based: int,
    args: argparse.Namespace,
) -> AlignmentResult:
    # Reverse-primer reads are anchored at the amplicon right end. Use only a suffix
    # plus upstream padding, while allowing a free reference prefix within that suffix.
    ref_end = len(ref)
    seg_start = max(1, ref_end - len(read_seq) - args.alignment_upstream_padding + 1)
    ref_seg = ref[seg_start - 1:ref_end]
    r_len = len(ref_seg)
    q_len = len(read_seq)
    gap = args.gap_cost
    mismatch = args.mismatch_cost

    dp = [[0] * (q_len + 1) for _ in range(r_len + 1)]
    op = [[""] * (q_len + 1) for _ in range(r_len + 1)]
    for i in range(1, r_len + 1):
        dp[i][0] = 0
        op[i][0] = "skip"
    for j in range(1, q_len + 1):
        dp[0][j] = j * gap
        op[0][j] = "ins"

    for i in range(1, r_len + 1):
        rb = ref_seg[i - 1]
        prev_row = dp[i - 1]
        row = dp[i]
        op_row = op[i]
        for j in range(1, q_len + 1):
            diag = prev_row[j - 1] + (0 if rb == read_seq[j - 1] else mismatch)
            delete = prev_row[j] + gap
            insert = row[j - 1] + gap
            best = diag
            best_op = "diag"
            if delete < best:
                best = delete
                best_op = "del"
            if insert < best:
                best = insert
                best_op = "ins"
            row[j] = best
            op_row[j] = best_op

    score = dp[r_len][q_len]
    if score > args.max_alignment_score:
        return empty_alignment("alignment_failed", "score_exceeds_threshold", score)

    i, j = r_len, q_len
    aligned_ref_rev: list[str] = []
    aligned_read_rev: list[str] = []
    aligned_qual_rev: list[int | None] = []
    ref_pos_rev: list[int | None] = []
    while j > 0:
        current_op = op[i][j]
        if current_op == "diag":
            aligned_ref_rev.append(ref_seg[i - 1])
            aligned_read_rev.append(read_seq[j - 1])
            aligned_qual_rev.append(ord(read_qual[j - 1]) - 33)
            ref_pos_rev.append(seg_start + i - 1)
            i -= 1
            j -= 1
        elif current_op == "del":
            aligned_ref_rev.append(ref_seg[i - 1])
            aligned_read_rev.append("-")
            aligned_qual_rev.append(None)
            ref_pos_rev.append(seg_start + i - 1)
            i -= 1
        elif current_op == "ins":
            aligned_ref_rev.append("-")
            aligned_read_rev.append(read_seq[j - 1])
            aligned_qual_rev.append(ord(read_qual[j - 1]) - 33)
            ref_pos_rev.append(None)
            j -= 1
        else:
            return empty_alignment("alignment_failed", "traceback_failed", score)

    aligned_ref = "".join(reversed(aligned_ref_rev))
    aligned_read = "".join(reversed(aligned_read_rev))
    aligned_qual = list(reversed(aligned_qual_rev))
    ref_positions = list(reversed(ref_pos_rev))
    covered_positions = [p for p in ref_positions if p is not None]
    if not covered_positions:
        return empty_alignment("alignment_failed", "no_reference_positions_aligned", score)
    ref_start = min(covered_positions)
    ref_end_aligned = max(covered_positions)
    if ref_start > cut_window[0] or ref_end_aligned < cut_window[1]:
        return empty_alignment("not_covering_cut_window", "cut_window_not_aligned", score)

    return summarize_alignment(
        score=score,
        ref_start=ref_start,
        ref_end=ref_end_aligned,
        aligned_ref=aligned_ref,
        aligned_read=aligned_read,
        aligned_qual=aligned_qual,
        ref_positions=ref_positions,
        ref=ref,
        cut_window=cut_window,
        context_window=context_window,
        codon_map=codon_map,
        coding_frame_start_1based=coding_frame_start_1based,
    )


def empty_alignment(status: str, reason: str, score: int | None) -> AlignmentResult:
    return AlignmentResult(
        status=status,
        reason=reason,
        score=score,
        ref_start_1based=None,
        ref_end_1based=None,
        aligned_ref="",
        aligned_read="",
        aligned_qual=[],
        ref_positions=[],
        n_inserted=0,
        n_deleted=0,
        n_substituted=0,
        n_inserted_cut=0,
        n_deleted_cut=0,
        n_substituted_cut=0,
        net_indel_length_cut=0,
        full_event_key="NA",
        allele_key="NA",
        cutsite_event_type="unclassified",
        frameshift_status="not_evaluated",
        consequence="not_evaluated",
        predicted_lof="not_evaluated",
        codon12_status="not_evaluated",
        codon13_status="not_evaluated",
        codon61_status="not_in_amplicon",
        known_oncogenic_status="not_evaluated",
        cut_window_min_baseq=None,
    )


def summarize_alignment(
    score: int,
    ref_start: int,
    ref_end: int,
    aligned_ref: str,
    aligned_read: str,
    aligned_qual: list[int | None],
    ref_positions: list[int | None],
    ref: str,
    cut_window: tuple[int, int],
    context_window: tuple[int, int],
    codon_map: dict[str, tuple[int | None, int | None]],
    coding_frame_start_1based: int,
) -> AlignmentResult:
    insertions: list[dict[str, object]] = []
    deletions: list[dict[str, object]] = []
    substitutions: list[dict[str, object]] = []
    window_qualities: list[int] = []

    previous_ref_pos = ref_start - 1
    idx = 0
    while idx < len(aligned_ref):
        rb = aligned_ref[idx]
        qb = aligned_read[idx]
        pos = ref_positions[idx]
        qv = aligned_qual[idx]
        if pos is not None:
            previous_ref_pos = pos
        if pos is not None and cut_window[0] <= pos <= cut_window[1] and qv is not None:
            window_qualities.append(qv)

        if rb == "-" and qb != "-":
            bases = []
            quals = []
            anchor = previous_ref_pos
            while idx < len(aligned_ref) and aligned_ref[idx] == "-" and aligned_read[idx] != "-":
                bases.append(aligned_read[idx])
                if aligned_qual[idx] is not None:
                    quals.append(aligned_qual[idx])
                idx += 1
            insertions.append({"anchor": anchor, "seq": "".join(bases), "quals": quals})
            if cut_window[0] <= anchor <= cut_window[1]:
                window_qualities.extend(quals)
            continue
        if rb != "-" and qb == "-":
            positions = []
            bases = []
            while idx < len(aligned_ref) and aligned_ref[idx] != "-" and aligned_read[idx] == "-":
                positions.append(ref_positions[idx])
                bases.append(aligned_ref[idx])
                idx += 1
            positions_int = [p for p in positions if p is not None]
            deletions.append({"start": min(positions_int), "end": max(positions_int), "seq": "".join(bases)})
            continue
        if rb != "-" and qb != "-" and rb != qb:
            substitutions.append({"pos": pos, "ref": rb, "alt": qb, "qual": qv})
        idx += 1

    n_inserted = sum(len(e["seq"]) for e in insertions)
    n_deleted = sum(int(e["end"]) - int(e["start"]) + 1 for e in deletions)
    n_substituted = len(substitutions)

    cut_insertions = [e for e in insertions if cut_window[0] <= int(e["anchor"]) <= cut_window[1]]
    cut_deletions = [e for e in deletions if not (int(e["end"]) < cut_window[0] or int(e["start"]) > cut_window[1])]
    cut_substitutions = [e for e in substitutions if e["pos"] is not None and cut_window[0] <= int(e["pos"]) <= cut_window[1]]
    n_inserted_cut = sum(len(e["seq"]) for e in cut_insertions)
    n_deleted_cut = sum(int(e["end"]) - int(e["start"]) + 1 for e in cut_deletions)
    n_substituted_cut = len(cut_substitutions)
    net_indel_length_cut = n_inserted_cut - n_deleted_cut

    context_events = events_to_strings(insertions, deletions, substitutions, context_window)
    full_events = events_to_strings(insertions, deletions, substitutions, (1, len(ref)))
    allele_key = ";".join(context_events) if context_events else "WT_context_81_101"
    full_event_key = ";".join(full_events) if full_events else "WT_observed_region"

    if cut_insertions or cut_deletions:
        cutsite_event_type = "indel"
    elif cut_substitutions:
        cutsite_event_type = "substitution_only"
    elif n_inserted or n_deleted or n_substituted:
        cutsite_event_type = "other_variant_outside_window"
    else:
        cutsite_event_type = "wt"

    frameshift_status = "not_applicable"
    consequence = "no_indel"
    predicted_lof = "no"
    if cut_insertions or cut_deletions:
        frameshift_status = "frameshift" if net_indel_length_cut % 3 else "in_frame"
        if frameshift_status == "frameshift":
            consequence = "predicted_lof"
            predicted_lof = "yes"
        else:
            mutated = apply_events_to_reference(ref, insertions, deletions, substitutions)
            aa = translate_amplicon(mutated, coding_frame_start_1based)
            if "*" in aa[:-1]:
                consequence = "predicted_lof"
                predicted_lof = "yes"
            else:
                consequence = "in_frame_unknown"

    covered_ref_positions = {p for p in ref_positions if p is not None}
    codon12_status = codon_status("12", codon_map, insertions, deletions, substitutions, ref, covered_ref_positions)
    codon13_status = codon_status("13", codon_map, insertions, deletions, substitutions, ref, covered_ref_positions)
    codon61_status = codon_status("61", codon_map, insertions, deletions, substitutions, ref, covered_ref_positions)
    known_oncogenic_status = "not_known_oncogenic"
    if codon12_status in {"G12D", "G12V", "G12C", "G12R", "G12S", "G12A"}:
        known_oncogenic_status = codon12_status
    elif codon13_status in {"G13D", "G13C", "G13R", "G13S", "G13V"}:
        known_oncogenic_status = codon13_status
    elif codon61_status not in {"not_in_amplicon", "not_covered", "WT", "unchanged"}:
        known_oncogenic_status = codon61_status

    return AlignmentResult(
        status="ok",
        reason="ok",
        score=score,
        ref_start_1based=ref_start,
        ref_end_1based=ref_end,
        aligned_ref=aligned_ref,
        aligned_read=aligned_read,
        aligned_qual=aligned_qual,
        ref_positions=ref_positions,
        n_inserted=n_inserted,
        n_deleted=n_deleted,
        n_substituted=n_substituted,
        n_inserted_cut=n_inserted_cut,
        n_deleted_cut=n_deleted_cut,
        n_substituted_cut=n_substituted_cut,
        net_indel_length_cut=net_indel_length_cut,
        full_event_key=full_event_key,
        allele_key=allele_key,
        cutsite_event_type=cutsite_event_type,
        frameshift_status=frameshift_status,
        consequence=consequence,
        predicted_lof=predicted_lof,
        codon12_status=codon12_status,
        codon13_status=codon13_status,
        codon61_status=codon61_status,
        known_oncogenic_status=known_oncogenic_status,
        cut_window_min_baseq=min(window_qualities) if window_qualities else None,
    )


def events_to_strings(
    insertions: list[dict[str, object]],
    deletions: list[dict[str, object]],
    substitutions: list[dict[str, object]],
    window: tuple[int, int],
) -> list[str]:
    events: list[str] = []
    for e in deletions:
        start, end = int(e["start"]), int(e["end"])
        if end < window[0] or start > window[1]:
            continue
        events.append(f"del:{start}-{end}:{e['seq']}")
    for e in insertions:
        anchor = int(e["anchor"])
        if window[0] <= anchor <= window[1]:
            events.append(f"ins_after:{anchor}:{e['seq']}")
    for e in substitutions:
        pos = int(e["pos"])
        if window[0] <= pos <= window[1]:
            events.append(f"sub:{pos}:{e['ref']}>{e['alt']}")
    return sorted(events, key=event_sort_key)


def event_sort_key(event: str) -> tuple[int, str]:
    parts = event.split(":")
    if parts[0] == "del":
        pos = int(parts[1].split("-")[0])
    elif parts[0] == "ins_after":
        pos = int(parts[1])
    else:
        pos = int(parts[1])
    return pos, event


def apply_events_to_reference(
    ref: str,
    insertions: list[dict[str, object]],
    deletions: list[dict[str, object]],
    substitutions: list[dict[str, object]],
) -> str:
    seq = list(ref)
    for e in substitutions:
        pos = int(e["pos"])
        if 1 <= pos <= len(seq):
            seq[pos - 1] = str(e["alt"])
    # Apply deletions and insertions from right to left to preserve coordinates.
    edits: list[tuple[int, str, object]] = []
    for e in deletions:
        edits.append((int(e["start"]), "del", e))
    for e in insertions:
        edits.append((int(e["anchor"]), "ins", e))
    for _pos, kind, e in sorted(edits, key=lambda x: x[0], reverse=True):
        if kind == "del":
            start, end = int(e["start"]), int(e["end"])
            del seq[start - 1:end]
        else:
            anchor = int(e["anchor"])
            ins_seq = list(str(e["seq"]))
            seq[anchor:anchor] = ins_seq
    return "".join(seq)


def translate_amplicon(seq: str, coding_frame_start_1based: int) -> str:
    aas = []
    start = coding_frame_start_1based - 1
    for idx in range(start, len(seq) - 2, 3):
        codon = seq[idx:idx + 3]
        aas.append(CODON_TABLE.get(codon, "X"))
    return "".join(aas)


def codon_status(
    codon: str,
    codon_map: dict[str, tuple[int | None, int | None]],
    insertions: list[dict[str, object]],
    deletions: list[dict[str, object]],
    substitutions: list[dict[str, object]],
    ref: str,
    covered_ref_positions: set[int],
) -> str:
    start, end = codon_map[codon]
    if start is None or end is None:
        return "not_in_amplicon"
    if not all(pos in covered_ref_positions for pos in range(start, end + 1)):
        return "not_covered_by_selected_read"
    overlapping_indel = False
    for e in deletions:
        if not (int(e["end"]) < start or int(e["start"]) > end):
            overlapping_indel = True
    for e in insertions:
        if start <= int(e["anchor"]) <= end:
            overlapping_indel = True
    overlapping_subs = [e for e in substitutions if start <= int(e["pos"]) <= end]
    if overlapping_indel:
        return "codon_disrupted_by_indel"
    if not overlapping_subs:
        return "WT"
    codon_seq = list(ref[start - 1:end])
    for e in overlapping_subs:
        codon_seq[int(e["pos"]) - start] = str(e["alt"])
    wt_aa = CODON_TABLE.get(ref[start - 1:end], "X")
    alt_aa = CODON_TABLE.get("".join(codon_seq), "X")
    if wt_aa == alt_aa:
        return "synonymous"
    return f"{wt_aa}{codon}{alt_aa}"


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            h.update(block)
    return h.hexdigest()


def command_output(cmd: list[str], cwd: Path | None = None) -> str:
    try:
        return subprocess.check_output(cmd, cwd=cwd, text=True, stderr=subprocess.STDOUT).strip()
    except Exception as exc:
        return f"NA ({type(exc).__name__}: {exc})"


def fmt_pct(numer: int | float, denom: int | float) -> str:
    if not denom:
        return "NA"
    return f"{100.0 * float(numer) / float(denom):.6f}"


def fmt_value(value: object) -> str:
    if value is None:
        return "NA"
    return str(value)


def classify_pair(
    aln: AlignmentResult,
    paired_hdr_status: str,
    min_baseq: int,
    conflict: bool,
) -> tuple[str, str]:
    if aln.status != "ok":
        return "ambiguous_low_confidence", aln.reason
    if aln.cut_window_min_baseq is not None and aln.cut_window_min_baseq < min_baseq:
        return "ambiguous_low_confidence", "cut_window_low_base_quality"
    if conflict:
        return "ambiguous_low_confidence", "conflicting_high_quality_reverse_mates"
    if paired_hdr_status == "hdr_exact":
        return "hdr_kras_g12d", "hdr_barcode_exact"
    if aln.cutsite_event_type == "indel":
        return "cutsite_indel", "cutsite_indel"
    if aln.cutsite_event_type == "substitution_only":
        return "substitution_only_cutsite_window", "cutsite_substitution_only"
    if aln.cutsite_event_type == "other_variant_outside_window":
        return "other_variant_outside_window", "variant_outside_cut_window"
    return "wt_unmodified_cutsite", "wt_cutsite"


def build_fieldnames() -> tuple[list[str], list[str], list[str], list[str], list[str]]:
    per_pair = [
        "specimen_id", "cohort", "sample_name", "analysis_name", "read_id",
        "r1_primer_status", "r2_primer_status", "amplicon_matched",
        "paired_hdr_status", "selected_read", "selected_read_length", "primer_offset",
        "nominal_start_1based", "alignment_status", "ambiguity_reason", "allele_class",
        "allele_key", "full_event_key", "cutsite_event_type", "alignment_score",
        "ref_start_1based", "ref_end_1based", "cut_window_min_baseq",
        "n_inserted", "n_deleted", "n_substituted", "n_inserted_cut", "n_deleted_cut",
        "n_substituted_cut", "net_indel_length_cut", "frameshift_status", "consequence",
        "predicted_lof", "codon12_status", "codon13_status", "codon61_status",
        "known_oncogenic_status",
    ]
    per_allele = [
        "specimen_id", "cohort", "sample_name", "analysis_name", "allele_class",
        "allele_key", "full_event_key", "cutsite_event_type", "frameshift_status",
        "consequence", "predicted_lof", "codon12_status", "codon13_status", "codon61_status",
        "known_oncogenic_status", "count", "pct_of_raw_pairs", "pct_of_classifiable_pairs",
    ]
    sample_summary = [
        "specimen_id", "cohort", "sample_name", "analysis_name", "raw_pairs", "r1_records",
        "r2_records", "pair_id_mismatches", "amplicon_matched_pairs", "cutsite_covering_pairs",
        "cutsite_classifiable_pairs", "alignment_failure_pairs", "low_quality_pairs",
        "conflicting_reverse_mate_pairs", "ambiguous_low_confidence_pairs", "wt_unmodified_cutsite_pairs",
        "cutsite_indel_pairs", "substitution_only_cutsite_window_pairs", "other_variant_outside_window_pairs",
        "hdr_kras_g12d_pair_calls", "non_hdr_classifiable_pairs", "non_hdr_wt_unmodified_cutsite_pairs",
        "non_hdr_cutsite_indel_pairs", "non_hdr_substitution_only_cutsite_window_pairs",
        "non_hdr_other_variant_outside_window_pairs", "paired_hdr_exact_pairs", "paired_hdr_wt_exact_pairs",
        "paired_hdr_other_pattern_pairs", "paired_hdr_short_pairs", "paired_hdr_no_forward_primer_pairs",
        "existing_kras_hdr_exact_pct", "existing_kras_hdr_exact_reads", "existing_kras_hdr_denominator",
        "existing_kras_qc_flag", "cutsite_classifiable_pct_raw", "alignment_failure_pct_cutsite_covering",
        "ambiguous_low_confidence_pct_raw", "ambiguous_low_confidence_pct_cutsite_covering",
        "cutsite_indel_pct_classifiable", "non_hdr_cutsite_indel_pct_non_hdr_classifiable",
        "wt_unmodified_cutsite_pct_classifiable", "qc_flag",
    ]
    consequence = [
        "specimen_id", "cohort", "sample_name", "analysis_name", "consequence", "frameshift_status",
        "count", "pct_of_cutsite_indel_pairs", "pct_of_classifiable_pairs",
    ]
    figure = [
        "specimen_id", "cohort", "sample_name", "analysis_name", "panel", "category",
        "count", "denominator", "denominator_name", "pct",
    ]
    return per_pair, per_allele, sample_summary, consequence, figure


def main() -> None:
    args = parse_args()
    configure_logging(args.log_level)
    start_time = time.time()
    outdir = Path(args.outdir)
    per_pair_dir = outdir / "per_pair"
    per_allele_dir = outdir / "per_allele"
    summary_dir = outdir / "summary"
    metadata_dir = outdir / "metadata"
    for directory in [per_pair_dir, per_allele_dir, summary_dir, metadata_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    amplicon_rows = read_tsv(Path(args.amplicons))
    amplicons = {row["locus"]: row for row in amplicon_rows}
    if args.locus not in amplicons:
        raise SystemExit(f"Locus {args.locus!r} not found in {args.amplicons}")
    amp = amplicons[args.locus]
    wt_amplicon = amp["wt_amplicon_seq"].upper()
    hdr_amplicon = amp["expected_hdr_amplicon_seq"].upper()
    fwd_kmer = amp["primer_fwd"].upper()[:args.primer_kmer_len]
    rev_kmer = amp["primer_rev"].upper()[:args.primer_kmer_len]
    wt_prefix = wt_amplicon[:args.prefix_length]
    hdr_prefix = hdr_amplicon[:args.prefix_length]
    informative_positions = [idx for idx, (w, h) in enumerate(zip(wt_prefix, hdr_prefix)) if w != h]
    wt_key = "".join(wt_prefix[pos] for pos in informative_positions)
    hdr_key = "".join(hdr_prefix[pos] for pos in informative_positions)
    cut_window = (args.cut_window_start, args.cut_window_end)
    context_window = (args.context_window_start, args.context_window_end)

    # Explicit Kras coding map for this amplicon: codon 12 = bases 21-23,
    # codon 13 = bases 24-26, codon 61 is outside the 134 bp amplicon.
    codon_map: dict[str, tuple[int | None, int | None]] = {
        "12": (21, 23),
        "13": (24, 26),
        "61": (None, None),
    }
    coding_frame_start_1based = 3

    manifest_rows = [row for row in read_tsv(Path(args.analysis_manifest)) if row["locus"] == args.locus]
    editing_rows = read_tsv(Path(args.editing_summary))
    editing_by_name = {row["kras_analysis_name"]: row for row in editing_rows if row.get("kras_analysis_name")}

    per_pair_fields, per_allele_fields, sample_fields, consequence_fields, figure_fields = build_fieldnames()
    per_pair_path = per_pair_dir / "kras_cutsite_pair_calls.tsv.gz"
    sample_rows: list[dict[str, object]] = []
    allele_counter: Counter[tuple[str, ...]] = Counter()
    allele_meta: dict[tuple[str, ...], dict[str, object]] = {}
    consequence_counter: Counter[tuple[str, str, str]] = Counter()
    figure_rows: list[dict[str, object]] = []

    with gzip.open(per_pair_path, "wt", newline="") as pair_handle:
        pair_writer = csv.DictWriter(pair_handle, fieldnames=per_pair_fields, delimiter="\t")
        pair_writer.writeheader()

        for row in manifest_rows:
            sample_start = time.time()
            analysis_name = row["analysis_name"]
            logging.info("Processing %s", analysis_name)
            r1_path = Path(row["remote_fastq_r1"])
            r2_path = Path(row["remote_fastq_r2"])
            if not r1_path.exists() or not r2_path.exists():
                raise FileNotFoundError(f"Missing FASTQ for {analysis_name}: {r1_path}, {r2_path}")

            stats: Counter[str] = Counter()
            hdr_status_counts: Counter[str] = Counter()
            class_counts: Counter[str] = Counter()
            non_hdr_class_counts: Counter[str] = Counter()
            consequence_counts: Counter[tuple[str, str]] = Counter()
            alignment_cache: dict[str, AlignmentResult] = {}

            r1_iter = iter_fastq(r1_path)
            r2_iter = iter_fastq(r2_path)
            for r1_record, r2_record in zip_longest(r1_iter, r2_iter):
                if r1_record is None:
                    stats["r2_records"] += 1
                    stats["unpaired_r2_records"] += 1
                    continue
                if r2_record is None:
                    stats["r1_records"] += 1
                    stats["unpaired_r1_records"] += 1
                    continue

                if args.max_pairs_per_sample and stats["raw_pairs"] >= args.max_pairs_per_sample:
                    break
                r1_id, r1_seq, r1_qual = r1_record
                r2_id, r2_seq, r2_qual = r2_record
                stats["raw_pairs"] += 1
                stats["r1_records"] += 1
                stats["r2_records"] += 1
                if r1_id != r2_id:
                    stats["pair_id_mismatches"] += 1

                r1_primer = primer_side_status(r1_seq, fwd_kmer, rev_kmer, args.primer_search_window)
                r2_primer = primer_side_status(r2_seq, fwd_kmer, rev_kmer, args.primer_search_window)
                amplicon_matched = r1_primer != "unmatched" or r2_primer != "unmatched"
                if amplicon_matched:
                    stats["amplicon_matched_pairs"] += 1

                forward_reads = []
                forward_reads.extend(forward_oriented_sequences(r1_seq, r1_qual, fwd_kmer, args.primer_search_window))
                forward_reads.extend(forward_oriented_sequences(r2_seq, r2_qual, fwd_kmer, args.primer_search_window))
                paired_hdr_status = classify_hdr_barcode(forward_reads, informative_positions, wt_key, hdr_key)
                hdr_status_counts[paired_hdr_status] += 1

                candidates = []
                candidates.extend(reverse_candidates(
                    "R1", r1_seq, r1_qual, rev_kmer, args.primer_search_window,
                    len(wt_amplicon), args.cut_window_start, args.cut_window_end,
                ))
                candidates.extend(reverse_candidates(
                    "R2", r2_seq, r2_qual, rev_kmer, args.primer_search_window,
                    len(wt_amplicon), args.cut_window_start, args.cut_window_end,
                ))

                selected: Candidate | None = None
                aln = empty_alignment("not_covering_cut_window", "no_reverse_primer_cutsite_read", None)
                conflict = False
                if candidates:
                    stats["cutsite_covering_pairs"] += 1
                    candidates.sort(
                        key=lambda c: (c.trimmed_len, c.nominal_window_minq if c.nominal_window_minq is not None else -1),
                        reverse=True,
                    )
                    selected = candidates[0]
                    forward_seq = revcomp(selected.raw_reverse_seq)
                    forward_qual = selected.raw_reverse_qual[::-1]
                    aln = alignment_cache.get(forward_seq)
                    if aln is None:
                        aln = semiglobal_align_suffix(
                            forward_seq, forward_qual, wt_amplicon, cut_window, context_window,
                            codon_map, coding_frame_start_1based, args,
                        )
                        alignment_cache[forward_seq] = aln

                    if len(candidates) > 1 and aln.status == "ok" and (
                        aln.cut_window_min_baseq is None or aln.cut_window_min_baseq >= args.min_window_baseq
                    ):
                        for other in candidates[1:]:
                            other_seq = revcomp(other.raw_reverse_seq)
                            other_qual = other.raw_reverse_qual[::-1]
                            other_aln = alignment_cache.get(other_seq)
                            if other_aln is None:
                                other_aln = semiglobal_align_suffix(
                                    other_seq, other_qual, wt_amplicon, cut_window, context_window,
                                    codon_map, coding_frame_start_1based, args,
                                )
                                alignment_cache[other_seq] = other_aln
                            if other_aln.status == "ok" and (
                                other_aln.cut_window_min_baseq is None
                                or other_aln.cut_window_min_baseq >= args.min_window_baseq
                            ) and other_aln.allele_key != aln.allele_key:
                                conflict = True
                                break

                allele_class, ambiguity_reason = classify_pair(aln, paired_hdr_status, args.min_window_baseq, conflict)
                codon12_out = "G12D_hdr_barcode" if paired_hdr_status == "hdr_exact" else aln.codon12_status
                codon13_out = aln.codon13_status
                codon61_out = aln.codon61_status
                known_oncogenic_out = "G12D_hdr_barcode" if paired_hdr_status == "hdr_exact" else aln.known_oncogenic_status
                class_counts[allele_class] += 1
                if paired_hdr_status != "hdr_exact" and allele_class != "ambiguous_low_confidence":
                    non_hdr_class_counts[allele_class] += 1

                if allele_class == "ambiguous_low_confidence":
                    stats["ambiguous_low_confidence_pairs"] += 1
                    if aln.status == "alignment_failed":
                        stats["alignment_failure_pairs"] += 1
                    if ambiguity_reason == "cut_window_low_base_quality":
                        stats["low_quality_pairs"] += 1
                    if conflict:
                        stats["conflicting_reverse_mate_pairs"] += 1
                else:
                    stats["cutsite_classifiable_pairs"] += 1

                if allele_class != "ambiguous_low_confidence" and aln.cutsite_event_type == "indel":
                    consequence_counts[(aln.consequence, aln.frameshift_status)] += 1
                    consequence_counter[(analysis_name, aln.consequence, aln.frameshift_status)] += 1

                allele_key_tuple = (
                    analysis_name, allele_class, aln.allele_key, aln.full_event_key,
                    aln.cutsite_event_type, aln.frameshift_status, aln.consequence,
                    aln.predicted_lof, codon12_out, codon13_out, codon61_out,
                    known_oncogenic_out,
                )
                allele_counter[allele_key_tuple] += 1
                if allele_key_tuple not in allele_meta:
                    allele_meta[allele_key_tuple] = {
                        "specimen_id": row["specimen_id"],
                        "cohort": row["cohort"],
                        "sample_name": row["sample_name"],
                        "analysis_name": analysis_name,
                        "allele_class": allele_class,
                        "allele_key": aln.allele_key,
                        "full_event_key": aln.full_event_key,
                        "cutsite_event_type": aln.cutsite_event_type,
                        "frameshift_status": aln.frameshift_status,
                        "consequence": aln.consequence,
                        "predicted_lof": aln.predicted_lof,
                        "codon12_status": codon12_out,
                        "codon13_status": codon13_out,
                        "codon61_status": codon61_out,
                        "known_oncogenic_status": known_oncogenic_out,
                    }

                pair_writer.writerow({
                    "specimen_id": row["specimen_id"],
                    "cohort": row["cohort"],
                    "sample_name": row["sample_name"],
                    "analysis_name": analysis_name,
                    "read_id": r1_id,
                    "r1_primer_status": r1_primer,
                    "r2_primer_status": r2_primer,
                    "amplicon_matched": "yes" if amplicon_matched else "no",
                    "paired_hdr_status": paired_hdr_status,
                    "selected_read": selected.source_read if selected else "NA",
                    "selected_read_length": selected.trimmed_len if selected else "NA",
                    "primer_offset": selected.primer_offset if selected else "NA",
                    "nominal_start_1based": selected.nominal_start_1based if selected else "NA",
                    "alignment_status": aln.status,
                    "ambiguity_reason": ambiguity_reason,
                    "allele_class": allele_class,
                    "allele_key": aln.allele_key,
                    "full_event_key": aln.full_event_key,
                    "cutsite_event_type": aln.cutsite_event_type,
                    "alignment_score": fmt_value(aln.score),
                    "ref_start_1based": fmt_value(aln.ref_start_1based),
                    "ref_end_1based": fmt_value(aln.ref_end_1based),
                    "cut_window_min_baseq": fmt_value(aln.cut_window_min_baseq),
                    "n_inserted": aln.n_inserted,
                    "n_deleted": aln.n_deleted,
                    "n_substituted": aln.n_substituted,
                    "n_inserted_cut": aln.n_inserted_cut,
                    "n_deleted_cut": aln.n_deleted_cut,
                    "n_substituted_cut": aln.n_substituted_cut,
                    "net_indel_length_cut": aln.net_indel_length_cut,
                    "frameshift_status": aln.frameshift_status,
                    "consequence": aln.consequence,
                    "predicted_lof": aln.predicted_lof,
                    "codon12_status": codon12_out,
                    "codon13_status": codon13_out,
                    "codon61_status": codon61_out,
                    "known_oncogenic_status": known_oncogenic_out,
                })

            existing = editing_by_name.get(analysis_name, {})
            raw_pairs = stats["raw_pairs"]
            classifiable = stats["cutsite_classifiable_pairs"]
            cutsite_covering = stats["cutsite_covering_pairs"]
            non_hdr_classifiable = sum(non_hdr_class_counts.values())
            qc_flags = []
            if classifiable < args.min_classifiable_pairs:
                qc_flags.append("fail_low_cutsite_coverage")
            if cutsite_covering and 100 * stats["alignment_failure_pairs"] / cutsite_covering > 5:
                qc_flags.append("warn_high_alignment_failure")
            if cutsite_covering and 100 * stats["ambiguous_low_confidence_pairs"] / cutsite_covering > 20:
                qc_flags.append("fail_high_ambiguous_rate")
            elif cutsite_covering and 100 * stats["ambiguous_low_confidence_pairs"] / cutsite_covering > 10:
                qc_flags.append("warn_high_ambiguous_rate")
            if stats["pair_id_mismatches"]:
                qc_flags.append("warn_pair_id_mismatches")
            qc_flag = ";".join(qc_flags) if qc_flags else "pass"

            summary_row = {
                "specimen_id": row["specimen_id"],
                "cohort": row["cohort"],
                "sample_name": row["sample_name"],
                "analysis_name": analysis_name,
                "raw_pairs": raw_pairs,
                "r1_records": stats["r1_records"],
                "r2_records": stats["r2_records"],
                "pair_id_mismatches": stats["pair_id_mismatches"],
                "amplicon_matched_pairs": stats["amplicon_matched_pairs"],
                "cutsite_covering_pairs": cutsite_covering,
                "cutsite_classifiable_pairs": classifiable,
                "alignment_failure_pairs": stats["alignment_failure_pairs"],
                "low_quality_pairs": stats["low_quality_pairs"],
                "conflicting_reverse_mate_pairs": stats["conflicting_reverse_mate_pairs"],
                "ambiguous_low_confidence_pairs": stats["ambiguous_low_confidence_pairs"],
                "wt_unmodified_cutsite_pairs": class_counts["wt_unmodified_cutsite"],
                "cutsite_indel_pairs": class_counts["cutsite_indel"],
                "substitution_only_cutsite_window_pairs": class_counts["substitution_only_cutsite_window"],
                "other_variant_outside_window_pairs": class_counts["other_variant_outside_window"],
                "hdr_kras_g12d_pair_calls": class_counts["hdr_kras_g12d"],
                "non_hdr_classifiable_pairs": non_hdr_classifiable,
                "non_hdr_wt_unmodified_cutsite_pairs": non_hdr_class_counts["wt_unmodified_cutsite"],
                "non_hdr_cutsite_indel_pairs": non_hdr_class_counts["cutsite_indel"],
                "non_hdr_substitution_only_cutsite_window_pairs": non_hdr_class_counts["substitution_only_cutsite_window"],
                "non_hdr_other_variant_outside_window_pairs": non_hdr_class_counts["other_variant_outside_window"],
                "paired_hdr_exact_pairs": hdr_status_counts["hdr_exact"],
                "paired_hdr_wt_exact_pairs": hdr_status_counts["wt_exact"],
                "paired_hdr_other_pattern_pairs": hdr_status_counts["other_pattern"],
                "paired_hdr_short_pairs": hdr_status_counts["short"],
                "paired_hdr_no_forward_primer_pairs": hdr_status_counts["no_forward_primer_read"],
                "existing_kras_hdr_exact_pct": existing.get("kras_hdr_exact_pct", "NA"),
                "existing_kras_hdr_exact_reads": existing.get("kras_hdr_exact_reads", "NA"),
                "existing_kras_hdr_denominator": existing.get("kras_informative_reads", "NA"),
                "existing_kras_qc_flag": existing.get("kras_qc_flag", "NA"),
                "cutsite_classifiable_pct_raw": fmt_pct(classifiable, raw_pairs),
                "alignment_failure_pct_cutsite_covering": fmt_pct(stats["alignment_failure_pairs"], cutsite_covering),
                "ambiguous_low_confidence_pct_raw": fmt_pct(stats["ambiguous_low_confidence_pairs"], raw_pairs),
                "ambiguous_low_confidence_pct_cutsite_covering": fmt_pct(stats["ambiguous_low_confidence_pairs"], cutsite_covering),
                "cutsite_indel_pct_classifiable": fmt_pct(class_counts["cutsite_indel"], classifiable),
                "non_hdr_cutsite_indel_pct_non_hdr_classifiable": fmt_pct(non_hdr_class_counts["cutsite_indel"], non_hdr_classifiable),
                "wt_unmodified_cutsite_pct_classifiable": fmt_pct(class_counts["wt_unmodified_cutsite"], classifiable),
                "qc_flag": qc_flag,
            }
            sample_rows.append(summary_row)

            for (cons, frame), count in sorted(consequence_counts.items()):
                consequence_counter[(analysis_name, cons, frame)] = count

            hdr_pct = existing.get("kras_hdr_exact_pct", "NA")
            hdr_den = existing.get("kras_informative_reads", "NA")
            hdr_reads = existing.get("kras_hdr_exact_reads", "NA")
            if hdr_pct != "NA" and hdr_den != "NA":
                figure_rows.append({
                    "specimen_id": row["specimen_id"], "cohort": row["cohort"],
                    "sample_name": row["sample_name"], "analysis_name": analysis_name,
                    "panel": "hdr_vs_non_hdr", "category": "HDR_KrasG12D_existing_prefix67",
                    "count": hdr_reads, "denominator": hdr_den, "denominator_name": "existing_kras_informative_reads",
                    "pct": hdr_pct,
                })
                figure_rows.append({
                    "specimen_id": row["specimen_id"], "cohort": row["cohort"],
                    "sample_name": row["sample_name"], "analysis_name": analysis_name,
                    "panel": "hdr_vs_non_hdr", "category": "non_HDR_existing_prefix67",
                    "count": "NA", "denominator": hdr_den, "denominator_name": "existing_kras_informative_reads",
                    "pct": f"{100.0 - float(hdr_pct):.6f}",
                })
            for category, count in [
                ("WT_cutsite_non_HDR", non_hdr_class_counts["wt_unmodified_cutsite"]),
                ("indel_cutsite_non_HDR", non_hdr_class_counts["cutsite_indel"]),
                ("substitution_only_cutsite_non_HDR", non_hdr_class_counts["substitution_only_cutsite_window"]),
                ("other_variant_outside_window_non_HDR", non_hdr_class_counts["other_variant_outside_window"]),
                ("ambiguous_low_confidence", stats["ambiguous_low_confidence_pairs"]),
            ]:
                denom = non_hdr_classifiable if category != "ambiguous_low_confidence" else raw_pairs
                figure_rows.append({
                    "specimen_id": row["specimen_id"], "cohort": row["cohort"],
                    "sample_name": row["sample_name"], "analysis_name": analysis_name,
                    "panel": "non_hdr_cutsite_breakdown", "category": category,
                    "count": count, "denominator": denom,
                    "denominator_name": "non_hdr_classifiable_pairs" if category != "ambiguous_low_confidence" else "raw_pairs",
                    "pct": fmt_pct(count, denom),
                })
            total_indel = class_counts["cutsite_indel"]
            for (cons, frame), count in sorted(consequence_counts.items()):
                figure_rows.append({
                    "specimen_id": row["specimen_id"], "cohort": row["cohort"],
                    "sample_name": row["sample_name"], "analysis_name": analysis_name,
                    "panel": "indel_consequence", "category": f"{cons}:{frame}",
                    "count": count, "denominator": total_indel,
                    "denominator_name": "cutsite_indel_pairs", "pct": fmt_pct(count, total_indel),
                })

            logging.info(
                "Finished %s: raw_pairs=%s classifiable=%s indels=%s qc=%s elapsed=%.1fs unique_alignments=%s",
                analysis_name, raw_pairs, classifiable, class_counts["cutsite_indel"], qc_flag,
                time.time() - sample_start, len(alignment_cache),
            )

    allele_rows: list[dict[str, object]] = []
    raw_pairs_by_analysis = {row["analysis_name"]: int(row["raw_pairs"]) for row in sample_rows}
    classifiable_by_analysis = {row["analysis_name"]: int(row["cutsite_classifiable_pairs"]) for row in sample_rows}
    for key, count in allele_counter.most_common():
        meta = dict(allele_meta[key])
        analysis_name = str(meta["analysis_name"])
        meta["count"] = count
        meta["pct_of_raw_pairs"] = fmt_pct(count, raw_pairs_by_analysis.get(analysis_name, 0))
        meta["pct_of_classifiable_pairs"] = fmt_pct(count, classifiable_by_analysis.get(analysis_name, 0))
        allele_rows.append(meta)

    consequence_rows: list[dict[str, object]] = []
    sample_lookup = {row["analysis_name"]: row for row in sample_rows}
    indel_by_analysis = {row["analysis_name"]: int(row["cutsite_indel_pairs"]) for row in sample_rows}
    for (analysis_name, cons, frame), count in sorted(consequence_counter.items()):
        sample = sample_lookup[analysis_name]
        consequence_rows.append({
            "specimen_id": sample["specimen_id"],
            "cohort": sample["cohort"],
            "sample_name": sample["sample_name"],
            "analysis_name": analysis_name,
            "consequence": cons,
            "frameshift_status": frame,
            "count": count,
            "pct_of_cutsite_indel_pairs": fmt_pct(count, indel_by_analysis.get(analysis_name, 0)),
            "pct_of_classifiable_pairs": fmt_pct(count, classifiable_by_analysis.get(analysis_name, 0)),
        })

    write_tsv(per_allele_dir / "kras_cutsite_allele_counts.tsv", allele_rows, per_allele_fields)
    write_tsv(summary_dir / "kras_cutsite_sample_summary.tsv", sample_rows, sample_fields)
    write_tsv(summary_dir / "kras_cutsite_consequence_summary.tsv", consequence_rows, consequence_fields)
    write_tsv(outdir / "figures" / "kras_cutsite_figure_ready.tsv", figure_rows, figure_fields)

    run_manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(sys.argv),
        "parameters": vars(args),
        "inputs": {
            "analysis_manifest": str(Path(args.analysis_manifest).resolve()),
            "amplicons": str(Path(args.amplicons).resolve()),
            "editing_summary": str(Path(args.editing_summary).resolve()),
        },
        "input_sha256": {
            "analysis_manifest": sha256_file(Path(args.analysis_manifest)),
            "amplicons": sha256_file(Path(args.amplicons)),
            "editing_summary": sha256_file(Path(args.editing_summary)),
        },
        "outputs": {
            "per_pair": str(per_pair_path),
            "per_allele": str(per_allele_dir / "kras_cutsite_allele_counts.tsv"),
            "sample_summary": str(summary_dir / "kras_cutsite_sample_summary.tsv"),
            "consequence_summary": str(summary_dir / "kras_cutsite_consequence_summary.tsv"),
            "figure_ready": str(outdir / "figures" / "kras_cutsite_figure_ready.tsv"),
        },
        "amplicon": {
            "locus": args.locus,
            "wt_amplicon_length": len(wt_amplicon),
            "cut_site_after_base_1based": args.cut_site,
            "cut_window_1based_inclusive": list(cut_window),
            "context_window_1based_inclusive": list(context_window),
            "guide_seq": amp["guide_seq"],
            "guide_orientation_in_amplicon": "reverse_complement",
            "coding_frame_start_1based": coding_frame_start_1based,
            "codon_map": {k: list(v) for k, v in codon_map.items()},
        },
        "method_notes": [
            "Reverse-primer-side reads are reverse-complemented into amplicon-forward coordinates.",
            "One cut-site-covering read is selected per pair, preferring longer reads then higher nominal cut-window quality.",
            "Alignment is anchored to the right amplicon end with a free reference prefix to avoid terminal truncation overcalling.",
            "Existing 67 bp Kras HDR barcode percentages are retained as the primary HDR metric.",
        ],
    }
    with (metadata_dir / "run_manifest.yaml").open("w") as handle:
        write_simple_yaml(run_manifest, handle)

    run_info = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": round(time.time() - start_time, 3),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python": sys.version,
        "executable": sys.executable,
        "cwd": os.getcwd(),
        "git_head": command_output(["git", "rev-parse", "HEAD"], Path.cwd()),
        "python_freeze": command_output([sys.executable, "-m", "pip", "freeze"]),
        "sample_count": len(sample_rows),
    }
    (metadata_dir / "run_info.json").write_text(json.dumps(run_info, indent=2) + "\n")
    logging.info("Wrote outputs under %s", outdir)


def write_simple_yaml(value: object, handle, indent: int = 0) -> None:
    prefix = " " * indent
    if isinstance(value, dict):
        for key, val in value.items():
            if isinstance(val, (dict, list)):
                handle.write(f"{prefix}{key}:\n")
                write_simple_yaml(val, handle, indent + 2)
            else:
                handle.write(f"{prefix}{key}: {yaml_scalar(val)}\n")
    elif isinstance(value, list):
        for item in value:
            if isinstance(item, (dict, list)):
                handle.write(f"{prefix}-\n")
                write_simple_yaml(item, handle, indent + 2)
            else:
                handle.write(f"{prefix}- {yaml_scalar(item)}\n")
    else:
        handle.write(f"{prefix}{yaml_scalar(value)}\n")


def yaml_scalar(value: object) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    text = str(value)
    if not text or any(ch in text for ch in ":#[]{}\n\t") or text.lower() in {"null", "true", "false", "na"}:
        return json.dumps(text)
    return text


if __name__ == "__main__":
    main()
