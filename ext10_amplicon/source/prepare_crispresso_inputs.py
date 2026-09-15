#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare analysis manifest and CRISPRessoBatch settings files "
            "for the Alignment workflow."
        )
    )
    parser.add_argument("--sample-manifest", required=True)
    parser.add_argument("--amplicons", required=True)
    parser.add_argument("--remote-fastq-dir", required=True)
    parser.add_argument("--remote-processed-dir", required=True)
    parser.add_argument("--analysis-manifest-out", required=True)
    parser.add_argument("--kras-batch-out", required=True)
    parser.add_argument("--trp53-batch-out", required=True)
    parser.add_argument("--prefix-length", type=int, default=67)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()

    sample_rows = read_tsv(Path(args.sample_manifest))
    amplicon_rows = read_tsv(Path(args.amplicons))
    amplicons = {row["locus"]: row for row in amplicon_rows}

    analysis_rows: list[dict[str, str]] = []
    kras_rows: list[dict[str, str]] = []
    trp53_rows: list[dict[str, str]] = []

    remote_fastq_dir = args.remote_fastq_dir.rstrip("/")
    remote_processed_dir = args.remote_processed_dir.rstrip("/")
    prefix_length = args.prefix_length

    for row in sample_rows:
        sample_name = row["sample_name"]
        locus = row["locus"]
        specimen_id = row["specimen_id"]
        remote_r1 = f"{remote_fastq_dir}/{row['r1']}"
        remote_r2 = f"{remote_fastq_dir}/{row['r2']}"
        analysis_name = f"{specimen_id}_{locus}"
        processed_fastq = f"{remote_processed_dir}/{analysis_name}.forward.fastq.gz"
        summary_json = f"{remote_processed_dir}/{analysis_name}.summary.json"

        is_excluded = sample_name == "L742"
        analysis_rows.append(
            {
                "specimen_id": specimen_id,
                "cohort": row["cohort"],
                "locus": locus,
                "sample_name": sample_name,
                "r1": row["r1"],
                "r2": row["r2"],
                "remote_fastq_r1": remote_r1,
                "remote_fastq_r2": remote_r2,
                "analysis_read": "R1",
                "analysis_fastq": processed_fastq,
                "analysis_name": analysis_name,
                "processed_summary_json": summary_json,
                "fastq_status": "excluded" if is_excluded else "present",
                "analysis_include": "no" if is_excluded else "yes",
                "exclusion_reason": "failed_qc_missing_fastq" if is_excluded else "",
            }
        )

        if is_excluded:
            continue

        amp = amplicons[locus]
        truncated_wt = amp["wt_amplicon_seq"][:prefix_length]
        batch_row = {
            "name": analysis_name,
            "fastq_r1": processed_fastq,
            "amplicon_seq": truncated_wt,
            "guide_seq": amp["guide_seq"].upper() if locus == "Trp53" else "",
            "amplicon_name": f"{locus}_prefix{prefix_length}",
        }

        if locus == "Kras":
            batch_row["expected_hdr_amplicon_seq"] = amp["expected_hdr_amplicon_seq"][:prefix_length]
            kras_rows.append(batch_row)
        elif locus == "Trp53":
            trp53_rows.append(batch_row)
        else:
            raise ValueError(f"Unexpected locus: {locus}")

    analysis_fieldnames = [
        "specimen_id",
        "cohort",
        "locus",
        "sample_name",
        "r1",
        "r2",
        "remote_fastq_r1",
        "remote_fastq_r2",
        "analysis_read",
        "analysis_fastq",
        "analysis_name",
        "processed_summary_json",
        "fastq_status",
        "analysis_include",
        "exclusion_reason",
    ]
    batch_fieldnames = ["name", "fastq_r1", "amplicon_seq", "guide_seq", "amplicon_name"]

    write_tsv(Path(args.analysis_manifest_out), analysis_rows, analysis_fieldnames)
    write_tsv(
        Path(args.kras_batch_out),
        kras_rows,
        batch_fieldnames + ["expected_hdr_amplicon_seq"],
    )
    write_tsv(Path(args.trp53_batch_out), trp53_rows, batch_fieldnames)


if __name__ == "__main__":
    main()
