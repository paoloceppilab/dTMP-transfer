#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
from datetime import datetime, timezone
from pathlib import Path


LOGGER = logging.getLogger("summarize_clean_editing_metrics")
COHORT_ORDER = {"Normal": 0, "KP": 1, "KPTT": 2}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize cleaned Kras HDR and Trp53 indel metrics from "
            "CRISPResso2 outputs and processed forward reads."
        )
    )
    parser.add_argument("--analysis-manifest", required=True)
    parser.add_argument("--amplicons", required=True)
    parser.add_argument("--processed-dir", required=True)
    parser.add_argument("--trp53-batch-dir", required=True)
    parser.add_argument("--editing-summary-out", required=True)
    parser.add_argument("--editing-summary-long-out", required=True)
    parser.add_argument("--qc-metrics-out", required=True)
    parser.add_argument("--metrics-json-out", required=True)
    parser.add_argument(
        "--prefix-length",
        type=int,
        default=67,
        help="Prefix length used for the adapted CRISPResso workflow.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    return parser.parse_args()


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(levelname)s: %(message)s",
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def specimen_sort_key(specimen_id: str, cohort: str) -> tuple[int, int, str]:
    suffix = specimen_id.split("_")[-1]
    try:
        numeric_suffix = int(suffix)
    except ValueError:
        numeric_suffix = 10**9
    return (COHORT_ORDER.get(cohort, 99), numeric_suffix, specimen_id)


def format_pct(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.6f}"


def format_int(value: int | None) -> str:
    if value is None:
        return "NA"
    return str(value)


def is_pure_right_edge_truncation(
    aligned_sequence: str,
    reference_sequence: str,
    n_deleted: int,
    n_inserted: int,
) -> bool:
    if n_inserted != 0 or n_deleted == 0:
        return False
    if not aligned_sequence.endswith("-"):
        return False
    trimmed_aligned = aligned_sequence.rstrip("-")
    trimmed_reference = reference_sequence[: len(trimmed_aligned)]
    if "-" in trimmed_aligned or "-" in trimmed_reference:
        return False
    return True


def summarize_kras_fastq(
    fastq_path: Path,
    wt_prefix: str,
    hdr_prefix: str,
    informative_positions: list[int],
) -> dict[str, int | float]:
    total_reads = 0
    covering_reads = 0
    hdr_exact_reads = 0
    wt_exact_reads = 0
    other_pattern_reads = 0
    short_reads = 0

    wt_key = "".join(wt_prefix[pos] for pos in informative_positions)
    hdr_key = "".join(hdr_prefix[pos] for pos in informative_positions)
    min_required_length = max(informative_positions) + 1

    with gzip.open(fastq_path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip()
            handle.readline()
            handle.readline()
            total_reads += 1
            if len(sequence) < min_required_length:
                short_reads += 1
                continue
            covering_reads += 1
            key = "".join(sequence[pos] for pos in informative_positions)
            if key == hdr_key:
                hdr_exact_reads += 1
            elif key == wt_key:
                wt_exact_reads += 1
            else:
                other_pattern_reads += 1

    hdr_exact_pct = 100.0 * hdr_exact_reads / covering_reads if covering_reads else None
    hdr_exact_pct_of_classified = (
        100.0 * hdr_exact_reads / (hdr_exact_reads + wt_exact_reads)
        if (hdr_exact_reads + wt_exact_reads)
        else None
    )
    other_pattern_pct = (
        100.0 * other_pattern_reads / covering_reads if covering_reads else None
    )

    return {
        "processed_reads": total_reads,
        "informative_reads": covering_reads,
        "hdr_exact_reads": hdr_exact_reads,
        "wt_exact_reads": wt_exact_reads,
        "other_pattern_reads": other_pattern_reads,
        "short_reads": short_reads,
        "hdr_exact_pct": hdr_exact_pct,
        "hdr_exact_pct_of_classified": hdr_exact_pct_of_classified,
        "other_pattern_pct": other_pattern_pct,
    }


def summarize_trp53_alleles(table_path: Path) -> dict[str, int | float]:
    total_window_reads = 0
    truncation_reads = 0
    indel_reads = 0
    true_indel_reads = 0
    substitution_only_reads = 0

    with table_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            reads = int(row["#Reads"])
            n_deleted = int(row["n_deleted"])
            n_inserted = int(row["n_inserted"])
            n_mutated = int(row["n_mutated"])
            total_window_reads += reads

            is_truncation = is_pure_right_edge_truncation(
                aligned_sequence=row["Aligned_Sequence"],
                reference_sequence=row["Reference_Sequence"],
                n_deleted=n_deleted,
                n_inserted=n_inserted,
            )
            if is_truncation:
                truncation_reads += reads

            has_indel = n_deleted > 0 or n_inserted > 0
            if has_indel:
                indel_reads += reads
                if not is_truncation:
                    true_indel_reads += reads
            elif n_mutated > 0:
                substitution_only_reads += reads

    usable_reads = total_window_reads - truncation_reads
    indel_pct_usable = (
        100.0 * true_indel_reads / usable_reads if usable_reads else None
    )
    truncation_pct_total = (
        100.0 * truncation_reads / total_window_reads if total_window_reads else None
    )

    return {
        "window_reads": total_window_reads,
        "truncation_reads": truncation_reads,
        "usable_reads": usable_reads,
        "indel_reads_raw": indel_reads,
        "true_indel_reads": true_indel_reads,
        "substitution_only_reads": substitution_only_reads,
        "indel_pct_usable": indel_pct_usable,
        "truncation_pct_total": truncation_pct_total,
    }


def build_qc_flag(locus: str, metrics: dict[str, int | float] | None, missing: bool) -> str:
    if missing:
        return "missing_after_qc"
    if metrics is None:
        return "missing_metrics"
    if locus == "Kras":
        other_pattern_pct = metrics["other_pattern_pct"]
        if isinstance(other_pattern_pct, float) and other_pattern_pct > 20.0:
            return "high_noncanonical_fraction"
        return "pass"
    truncation_pct_total = metrics["truncation_pct_total"]
    if isinstance(truncation_pct_total, float) and truncation_pct_total > 30.0:
        return "edge_truncation_present"
    return "pass"


def main() -> None:
    args = parse_args()
    configure_logging(args.log_level)

    analysis_manifest = read_tsv(Path(args.analysis_manifest))
    amplicon_rows = read_tsv(Path(args.amplicons))
    amplicons = {row["locus"]: row for row in amplicon_rows}
    processed_dir = Path(args.processed_dir)
    trp53_batch_dir = Path(args.trp53_batch_dir)

    kras_wt_prefix = amplicons["Kras"]["wt_amplicon_seq"][: args.prefix_length]
    kras_hdr_prefix = amplicons["Kras"]["expected_hdr_amplicon_seq"][: args.prefix_length]
    kras_informative_positions = [
        idx
        for idx, (wt_base, hdr_base) in enumerate(zip(kras_wt_prefix, kras_hdr_prefix))
        if wt_base != hdr_base
    ]
    trp53_guide = amplicons["Trp53"]["guide_seq"].upper()

    long_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    specimen_summary: dict[str, dict[str, object]] = {}

    for row in analysis_manifest:
        specimen_id = row["specimen_id"]
        cohort = row["cohort"]
        locus = row["locus"]
        analysis_name = row["analysis_name"]
        missing = row["analysis_include"] != "yes"

        if specimen_id not in specimen_summary:
            specimen_summary[specimen_id] = {
                "specimen_id": specimen_id,
                "cohort": cohort,
                "kras_analysis_name": "NA",
                "kras_hdr_exact_pct": "NA",
                "kras_hdr_exact_pct_of_classified": "NA",
                "kras_hdr_exact_reads": "NA",
                "kras_informative_reads": "NA",
                "kras_wt_exact_reads": "NA",
                "kras_other_pattern_reads": "NA",
                "kras_short_reads": "NA",
                "kras_qc_flag": "NA",
                "trp53_analysis_name": "NA",
                "trp53_indel_pct": "NA",
                "trp53_true_indel_reads": "NA",
                "trp53_usable_reads": "NA",
                "trp53_truncation_reads": "NA",
                "trp53_window_reads": "NA",
                "trp53_substitution_only_reads": "NA",
                "trp53_qc_flag": "NA",
            }

        if missing:
            metrics: dict[str, int | float] | None = None
        elif locus == "Kras":
            fastq_path = processed_dir / f"{analysis_name}.forward.fastq.gz"
            metrics = summarize_kras_fastq(
                fastq_path=fastq_path,
                wt_prefix=kras_wt_prefix,
                hdr_prefix=kras_hdr_prefix,
                informative_positions=kras_informative_positions,
            )
        elif locus == "Trp53":
            table_path = (
                trp53_batch_dir
                / f"CRISPResso_on_{analysis_name}"
                / f"Trp53_prefix67.Alleles_frequency_table_around_sgRNA_{trp53_guide}.txt"
            )
            metrics = summarize_trp53_alleles(table_path=table_path)
        else:
            raise ValueError(f"Unexpected locus: {locus}")

        qc_flag = build_qc_flag(locus=locus, metrics=metrics, missing=missing)

        if locus == "Kras":
            primary_metric_name = "kras_hdr_exact_pct"
            primary_metric_pct = metrics["hdr_exact_pct"] if metrics else None
            specimen_summary[specimen_id].update(
                {
                    "kras_analysis_name": analysis_name,
                    "kras_hdr_exact_pct": format_pct(
                        metrics["hdr_exact_pct"] if metrics else None
                    ),
                    "kras_hdr_exact_pct_of_classified": format_pct(
                        metrics["hdr_exact_pct_of_classified"] if metrics else None
                    ),
                    "kras_hdr_exact_reads": format_int(
                        metrics["hdr_exact_reads"] if metrics else None
                    ),
                    "kras_informative_reads": format_int(
                        metrics["informative_reads"] if metrics else None
                    ),
                    "kras_wt_exact_reads": format_int(
                        metrics["wt_exact_reads"] if metrics else None
                    ),
                    "kras_other_pattern_reads": format_int(
                        metrics["other_pattern_reads"] if metrics else None
                    ),
                    "kras_short_reads": format_int(
                        metrics["short_reads"] if metrics else None
                    ),
                    "kras_qc_flag": qc_flag,
                }
            )
            qc_rows.append(
                {
                    "specimen_id": specimen_id,
                    "cohort": cohort,
                    "locus": locus,
                    "analysis_name": analysis_name,
                    "processed_reads": format_int(
                        metrics["processed_reads"] if metrics else None
                    ),
                    "usable_reads": format_int(
                        metrics["informative_reads"] if metrics else None
                    ),
                    "primary_metric_name": primary_metric_name,
                    "primary_metric_pct": format_pct(primary_metric_pct),
                    "technical_artifact_reads": format_int(
                        metrics["short_reads"] if metrics else None
                    ),
                    "technical_artifact_pct": format_pct(
                        (
                            100.0 * metrics["short_reads"] / metrics["processed_reads"]
                            if metrics and metrics["processed_reads"]
                            else None
                        )
                    ),
                    "qc_flag": qc_flag,
                }
            )
        else:
            primary_metric_name = "trp53_indel_pct"
            primary_metric_pct = metrics["indel_pct_usable"] if metrics else None
            specimen_summary[specimen_id].update(
                {
                    "trp53_analysis_name": analysis_name,
                    "trp53_indel_pct": format_pct(
                        metrics["indel_pct_usable"] if metrics else None
                    ),
                    "trp53_true_indel_reads": format_int(
                        metrics["true_indel_reads"] if metrics else None
                    ),
                    "trp53_usable_reads": format_int(
                        metrics["usable_reads"] if metrics else None
                    ),
                    "trp53_truncation_reads": format_int(
                        metrics["truncation_reads"] if metrics else None
                    ),
                    "trp53_window_reads": format_int(
                        metrics["window_reads"] if metrics else None
                    ),
                    "trp53_substitution_only_reads": format_int(
                        metrics["substitution_only_reads"] if metrics else None
                    ),
                    "trp53_qc_flag": qc_flag,
                }
            )
            qc_rows.append(
                {
                    "specimen_id": specimen_id,
                    "cohort": cohort,
                    "locus": locus,
                    "analysis_name": analysis_name,
                    "processed_reads": "NA",
                    "usable_reads": format_int(
                        metrics["usable_reads"] if metrics else None
                    ),
                    "primary_metric_name": primary_metric_name,
                    "primary_metric_pct": format_pct(primary_metric_pct),
                    "technical_artifact_reads": format_int(
                        metrics["truncation_reads"] if metrics else None
                    ),
                    "technical_artifact_pct": format_pct(
                        metrics["truncation_pct_total"] if metrics else None
                    ),
                    "qc_flag": qc_flag,
                }
            )

        long_rows.append(
            {
                "specimen_id": specimen_id,
                "cohort": cohort,
                "locus": locus,
                "analysis_name": analysis_name,
                "sample_name": row["sample_name"],
                "fastq_status": row["fastq_status"],
                "analysis_include": row["analysis_include"],
                "exclusion_reason": row["exclusion_reason"] or "NA",
                "primary_metric_name": primary_metric_name,
                "primary_metric_pct": format_pct(primary_metric_pct),
                "processed_reads": format_int(
                    metrics["processed_reads"] if metrics and "processed_reads" in metrics else None
                ),
                "usable_reads": format_int(
                    (
                        metrics["informative_reads"]
                        if metrics and locus == "Kras"
                        else metrics["usable_reads"]
                        if metrics and locus == "Trp53"
                        else None
                    )
                ),
                "hdr_exact_reads": format_int(
                    metrics["hdr_exact_reads"] if metrics and locus == "Kras" else None
                ),
                "wt_exact_reads": format_int(
                    metrics["wt_exact_reads"] if metrics and locus == "Kras" else None
                ),
                "other_pattern_reads": format_int(
                    metrics["other_pattern_reads"] if metrics and locus == "Kras" else None
                ),
                "short_reads": format_int(
                    metrics["short_reads"] if metrics and locus == "Kras" else None
                ),
                "window_reads": format_int(
                    metrics["window_reads"] if metrics and locus == "Trp53" else None
                ),
                "truncation_reads": format_int(
                    metrics["truncation_reads"] if metrics and locus == "Trp53" else None
                ),
                "true_indel_reads": format_int(
                    metrics["true_indel_reads"] if metrics and locus == "Trp53" else None
                ),
                "substitution_only_reads": format_int(
                    metrics["substitution_only_reads"] if metrics and locus == "Trp53" else None
                ),
                "qc_flag": qc_flag,
            }
        )

        LOGGER.info(
            "Summarized %s (%s): %s=%s",
            analysis_name,
            locus,
            primary_metric_name,
            format_pct(primary_metric_pct),
        )

    summary_rows = [
        specimen_summary[specimen_id]
        for specimen_id in sorted(
            specimen_summary,
            key=lambda item: specimen_sort_key(
                specimen_summary[item]["specimen_id"],
                specimen_summary[item]["cohort"],
            ),
        )
    ]
    long_rows.sort(key=lambda row: specimen_sort_key(row["specimen_id"], row["cohort"]))
    qc_rows.sort(key=lambda row: specimen_sort_key(row["specimen_id"], row["cohort"]))

    write_tsv(
        Path(args.editing_summary_out),
        summary_rows,
        [
            "specimen_id",
            "cohort",
            "kras_analysis_name",
            "kras_hdr_exact_pct",
            "kras_hdr_exact_pct_of_classified",
            "kras_hdr_exact_reads",
            "kras_informative_reads",
            "kras_wt_exact_reads",
            "kras_other_pattern_reads",
            "kras_short_reads",
            "kras_qc_flag",
            "trp53_analysis_name",
            "trp53_indel_pct",
            "trp53_true_indel_reads",
            "trp53_usable_reads",
            "trp53_truncation_reads",
            "trp53_window_reads",
            "trp53_substitution_only_reads",
            "trp53_qc_flag",
        ],
    )
    write_tsv(
        Path(args.editing_summary_long_out),
        long_rows,
        [
            "specimen_id",
            "cohort",
            "locus",
            "analysis_name",
            "sample_name",
            "fastq_status",
            "analysis_include",
            "exclusion_reason",
            "primary_metric_name",
            "primary_metric_pct",
            "processed_reads",
            "usable_reads",
            "hdr_exact_reads",
            "wt_exact_reads",
            "other_pattern_reads",
            "short_reads",
            "window_reads",
            "truncation_reads",
            "true_indel_reads",
            "substitution_only_reads",
            "qc_flag",
        ],
    )
    write_tsv(
        Path(args.qc_metrics_out),
        qc_rows,
        [
            "specimen_id",
            "cohort",
            "locus",
            "analysis_name",
            "processed_reads",
            "usable_reads",
            "primary_metric_name",
            "primary_metric_pct",
            "technical_artifact_reads",
            "technical_artifact_pct",
            "qc_flag",
        ],
    )

    metrics_json = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "prefix_length": args.prefix_length,
        "kras": {
            "metric_name": "kras_hdr_exact_pct",
            "informative_positions_1based": [pos + 1 for pos in kras_informative_positions],
            "wt_bases": "".join(kras_wt_prefix[pos] for pos in kras_informative_positions),
            "hdr_bases": "".join(kras_hdr_prefix[pos] for pos in kras_informative_positions),
            "description": (
                "Percent of processed forward-oriented reads spanning all informative "
                "positions that exactly match the 8-base HDR barcode pattern."
            ),
        },
        "trp53": {
            "metric_name": "trp53_indel_pct",
            "description": (
                "Percent of reads in the guide-centered allele table with indels after "
                "excluding pure right-edge truncation artifacts."
            ),
            "truncation_rule": (
                "No insertions, at least one deletion, and the aligned sequence consists "
                "only of a trailing right-edge deletion with no internal gaps."
            ),
        },
    }
    metrics_path = Path(args.metrics_json_out)
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    metrics_path.write_text(json.dumps(metrics_json, indent=2))


if __name__ == "__main__":
    main()
