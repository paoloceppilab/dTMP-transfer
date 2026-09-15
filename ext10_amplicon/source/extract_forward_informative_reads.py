#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def find_primer_start(seq: str, primer_kmer: str, search_window: int) -> int:
    max_start = min(search_window, max(0, len(seq) - len(primer_kmer)))
    for start in range(max_start + 1):
        if seq[start : start + len(primer_kmer)] == primer_kmer:
            return start
    return -1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract forward-oriented informative amplicon reads from paired FASTQ files. "
            "Reads are retained only if the forward-primer k-mer is found near the read start; "
            "reverse-complemented reads are oriented to the forward direction before writing."
        )
    )
    parser.add_argument("--fastq-r1", required=True)
    parser.add_argument("--fastq-r2", required=True)
    parser.add_argument("--primer-fwd", required=True)
    parser.add_argument("--output-fastq", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--primer-kmer-len", type=int, default=10)
    parser.add_argument("--search-window", type=int, default=5)
    return parser.parse_args()


def iter_fastq(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip().upper()
            plus = handle.readline()
            qual = handle.readline().strip()
            yield header, seq, plus, qual


def process_record(
    header: str,
    seq: str,
    plus: str,
    qual: str,
    primer_kmer: str,
    search_window: int,
):
    start = find_primer_start(seq, primer_kmer, search_window)
    if start != -1:
        return header, seq[start:], plus, qual[start:], "forward"

    rc_seq = revcomp(seq)
    rc_qual = qual[::-1]
    rc_start = find_primer_start(rc_seq, primer_kmer, search_window)
    if rc_start != -1:
        return header, rc_seq[rc_start:], plus, rc_qual[rc_start:], "reverse_complemented"

    return None


def main() -> None:
    args = parse_args()
    primer_kmer = args.primer_fwd.strip().upper()[: args.primer_kmer_len]

    total_reads = 0
    kept_reads = 0
    forward_reads = 0
    reverse_complemented_reads = 0

    output_path = Path(args.output_fastq)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path = Path(args.summary_json)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(output_path, "wt") as out_handle:
        for fastq_path in [Path(args.fastq_r1), Path(args.fastq_r2)]:
            for header, seq, plus, qual in iter_fastq(fastq_path):
                total_reads += 1
                processed = process_record(
                    header,
                    seq,
                    plus,
                    qual,
                    primer_kmer,
                    args.search_window,
                )
                if processed is None:
                    continue

                kept_header, kept_seq, kept_plus, kept_qual, orientation = processed
                out_handle.write(kept_header)
                out_handle.write(kept_seq + "\n")
                out_handle.write(kept_plus)
                out_handle.write(kept_qual + "\n")
                kept_reads += 1

                if orientation == "forward":
                    forward_reads += 1
                else:
                    reverse_complemented_reads += 1

    summary = {
        "fastq_r1": args.fastq_r1,
        "fastq_r2": args.fastq_r2,
        "output_fastq": args.output_fastq,
        "primer_kmer": primer_kmer,
        "search_window": args.search_window,
        "total_reads_seen": total_reads,
        "kept_reads": kept_reads,
        "forward_reads": forward_reads,
        "reverse_complemented_reads": reverse_complemented_reads,
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
