#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import platform
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from plot_kras_cutsite_pm20_threebin_background83 import (
    COHORT_ORDER,
    classify_nonhdr_row,
    configure_style,
    is_full_window_classifiable,
    is_hdr_row,
    pct,
    sort_key,
)


FRAME_LABELS = ["frameshift_indel", "in_frame_indel"]
OPTIONAL_FRAME_LABELS = ["other_indel_frame_status"]
COLORS = {
    "frameshift_indel": "#D55E00",
    "in_frame_indel": "#E69F00",
    "other_indel_frame_status": "#8C8C8C",
}
PRETTY_LABELS = {
    "frameshift_indel": "frameshift indels",
    "in_frame_indel": "in-frame indels",
    "other_indel_frame_status": "other frame status",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot tumor-only Kras cleavage-site +/-20 bp indel reads split by frame status "
            "from existing per-pair cut-site calls. Uses the same classifiable non-HDR "
            "denominator and sub:83:T>C background handling as the three-bin plot."
        )
    )
    parser.add_argument("--per-pair")
    parser.add_argument("--outdir")
    parser.add_argument("--prefix", default="kras_cutsite_pm20_indel_frame_background83_tumor_only")
    parser.add_argument("--cut-site", type=int, default=91)
    parser.add_argument("--flank-bp", type=int, default=20)
    parser.add_argument("--min-window-baseq", type=int, default=20)
    parser.add_argument("--background-substitution", default="sub:83:T>C")
    parser.add_argument(
        "--tumor-only",
        action="store_true",
        help="Plot only non-Normal samples, while keeping all samples in the figure-ready table.",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run synthetic frame-status tests and exit without reading input tables.",
    )
    args = parser.parse_args()
    if args.self_test:
        return args
    missing = [name for name in ["per_pair", "outdir"] if getattr(args, name) is None]
    if missing:
        parser.error("Missing required arguments unless --self-test is used: " + ", ".join(f"--{m.replace('_', '-')}" for m in missing))
    return args


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def frame_label(row: dict[str, str]) -> str:
    status = row.get("frameshift_status", "")
    if status == "frameshift":
        return "frameshift_indel"
    if status == "in_frame":
        return "in_frame_indel"
    return "other_indel_frame_status"


def run_self_tests() -> None:
    tests = {
        "frameshift": (frame_label({"frameshift_status": "frameshift"}), "frameshift_indel"),
        "in_frame": (frame_label({"frameshift_status": "in_frame"}), "in_frame_indel"),
        "unknown": (frame_label({"frameshift_status": "not_evaluated"}), "other_indel_frame_status"),
    }
    failures = {name: {"observed": obs, "expected": exp} for name, (obs, exp) in tests.items() if obs != exp}
    if failures:
        raise AssertionError(json.dumps(failures, indent=2))
    print("Synthetic indel frame-status tests passed.", file=sys.stderr)


def main() -> None:
    args = parse_args()
    if args.self_test:
        run_self_tests()
        return

    configure_style()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cut_start = args.cut_site - args.flank_bp
    cut_end = args.cut_site + args.flank_bp

    counts: dict[tuple[str, str, str, str], int] = Counter()
    raw_by_sample: dict[tuple[str, str, str], int] = Counter()
    full_classifiable_by_sample: dict[tuple[str, str, str], int] = Counter()
    display_denominator_by_sample: dict[tuple[str, str, str], int] = Counter()
    indel_denominator_by_sample: dict[tuple[str, str, str], int] = Counter()
    qc_by_sample: dict[tuple[str, str, str], Counter[str]] = defaultdict(Counter)
    g12_status_by_sample_frame: dict[tuple[tuple[str, str, str], str], Counter[str]] = defaultdict(Counter)
    g13_status_by_sample_frame: dict[tuple[tuple[str, str, str], str], Counter[str]] = defaultdict(Counter)
    oncogenic_status_by_sample_frame: dict[tuple[tuple[str, str, str], str], Counter[str]] = defaultdict(Counter)

    opener = gzip.open if str(args.per_pair).endswith(".gz") else open
    with opener(args.per_pair, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["specimen_id"], row["cohort"], row["sample_name"])
            raw_by_sample[key] += 1
            if not is_full_window_classifiable(row, cut_start, cut_end, args.min_window_baseq):
                qc_by_sample[key]["not_full_window_classifiable_pairs"] += 1
                continue

            full_classifiable_by_sample[key] += 1
            if is_hdr_row(row):
                qc_by_sample[key]["hdr_excluded_pairs"] += 1
                continue

            display_denominator_by_sample[key] += 1
            display_call, flags = classify_nonhdr_row(row, cut_start, cut_end, args.background_substitution)
            if bool(flags["has_background_substitution"]):
                qc_by_sample[key]["background_substitution_nonhdr_classifiable_pairs"] += 1
            if display_call != "in_del":
                continue

            indel_denominator_by_sample[key] += 1
            label = frame_label(row)
            counts[(*key, label)] += 1
            qc_by_sample[key][f"{label}_pairs"] += 1
            g12_status_by_sample_frame[(key, label)][row["codon12_status"]] += 1
            g13_status_by_sample_frame[(key, label)][row["codon13_status"]] += 1
            oncogenic_status_by_sample_frame[(key, label)][row["known_oncogenic_status"]] += 1

    rows: list[dict[str, object]] = []
    for key in sorted(raw_by_sample, key=lambda x: sort_key({"specimen_id": x[0], "cohort": x[1]})):
        specimen, cohort, sample_name = key
        denominator = display_denominator_by_sample[key]
        indel_denominator = indel_denominator_by_sample[key]
        row: dict[str, object] = {
            "specimen_id": specimen,
            "cohort": cohort,
            "sample_name": sample_name,
            "plot_include": "yes" if (not args.tumor_only or cohort != "Normal") else "no",
            "raw_pairs": raw_by_sample[key],
            "not_full_window_classifiable_pairs": qc_by_sample[key]["not_full_window_classifiable_pairs"],
            "full_pm20_classifiable_pairs_including_HDR": full_classifiable_by_sample[key],
            "hdr_excluded_pairs": qc_by_sample[key]["hdr_excluded_pairs"],
            "displayed_nonhdr_classifiable_pairs": denominator,
            "indel_pairs": indel_denominator,
            "indel_pct_of_displayed_nonhdr_classifiable": pct(indel_denominator, denominator),
            "pm20_window_start": cut_start,
            "pm20_window_end": cut_end,
            "min_window_baseq": args.min_window_baseq,
            "background_substitution": args.background_substitution,
            "background_substitution_nonhdr_classifiable_pairs": qc_by_sample[key][
                "background_substitution_nonhdr_classifiable_pairs"
            ],
        }
        displayed_pct_sum = 0.0
        indel_pct_sum = 0.0
        for label in FRAME_LABELS + OPTIONAL_FRAME_LABELS:
            n = counts[(*key, label)]
            row[f"{label}_pairs"] = n
            row[f"{label}_pct_of_displayed_nonhdr_classifiable"] = pct(n, denominator)
            row[f"{label}_pct_of_indel_reads"] = pct(n, indel_denominator)
            displayed_pct_sum += float(row[f"{label}_pct_of_displayed_nonhdr_classifiable"])
            indel_pct_sum += float(row[f"{label}_pct_of_indel_reads"])
        row["frame_split_pct_of_displayed_nonhdr_classifiable_sum"] = displayed_pct_sum
        row["frame_split_pct_of_indel_reads_sum"] = indel_pct_sum if indel_denominator else 0.0
        row["hdr_denominator_check_pass"] = (
            full_classifiable_by_sample[key] - qc_by_sample[key]["hdr_excluded_pairs"]
        ) == denominator
        row["indel_count_check_pass"] = sum(int(row[f"{label}_pairs"]) for label in FRAME_LABELS + OPTIONAL_FRAME_LABELS) == indel_denominator
        rows.append(row)

    if not rows:
        raise ValueError("No samples found in per-pair table.")

    figure_ready = outdir / f"{args.prefix}_figure_ready.tsv"
    with figure_ready.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    plot_rows = [row for row in rows if row["plot_include"] == "yes"]
    if not plot_rows:
        raise ValueError("No samples selected for plotting.")
    labels_to_plot = list(FRAME_LABELS)
    if any(int(row["other_indel_frame_status_pairs"]) for row in plot_rows):
        labels_to_plot += OPTIONAL_FRAME_LABELS

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.24, top=0.78)
    x = np.arange(len(plot_rows))
    bottom = np.zeros(len(plot_rows))
    for label in labels_to_plot:
        vals = np.array([float(row[f"{label}_pct_of_displayed_nonhdr_classifiable"]) for row in plot_rows])
        ax.bar(
            x,
            vals,
            bottom=bottom,
            color=COLORS[label],
            edgecolor="white",
            linewidth=0.25,
            label=PRETTY_LABELS[label],
        )
        bottom += vals

    ymax = max(10.0, min(100.0, float(np.ceil((bottom.max() + 5.0) / 10.0) * 10.0)))
    ax.set_ylim(0, ymax)
    ax.set_ylabel("Classifiable non-HDR reads (%)")
    ax.set_xticks(x)
    ax.set_xticklabels([str(row["specimen_id"]) for row in plot_rows], rotation=60, ha="right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, ncol=len(labels_to_plot), loc="upper center", bbox_to_anchor=(0.5, 1.12))
    fig.suptitle(
        "KRAS cut-site +/-20 bp indel frame status",
        y=0.97,
        fontsize=10,
        fontweight="bold",
    )

    cohorts = [str(row["cohort"]) for row in plot_rows]
    start_idx = 0
    for idx in range(1, len(cohorts) + 1):
        if idx == len(cohorts) or cohorts[idx] != cohorts[start_idx]:
            end_idx = idx - 1
            if end_idx < len(cohorts) - 1:
                ax.axvline(end_idx + 0.5, color="#999999", linewidth=0.6)
            start_idx = idx

    outputs: dict[str, str] = {}
    for ext in ["pdf", "svg", "png"]:
        path = outdir / f"{args.prefix}.{ext}"
        fig.savefig(path, dpi=600 if ext == "png" else None, bbox_inches="tight")
        outputs[f"figure_{ext}"] = str(path)
    plt.close(fig)

    aggregate_denominator = sum(int(row["displayed_nonhdr_classifiable_pairs"]) for row in plot_rows)
    aggregate_indel_denominator = sum(int(row["indel_pairs"]) for row in plot_rows)
    aggregate_counts = {label: sum(int(row[f"{label}_pairs"]) for row in plot_rows) for label in FRAME_LABELS + OPTIONAL_FRAME_LABELS}
    aggregate_g12_status_by_frame: dict[str, Counter[str]] = defaultdict(Counter)
    aggregate_g13_status_by_frame: dict[str, Counter[str]] = defaultdict(Counter)
    aggregate_oncogenic_status_by_frame: dict[str, Counter[str]] = defaultdict(Counter)
    for row in plot_rows:
        key = (str(row["specimen_id"]), str(row["cohort"]), str(row["sample_name"]))
        for label in FRAME_LABELS + OPTIONAL_FRAME_LABELS:
            aggregate_g12_status_by_frame[label].update(g12_status_by_sample_frame[(key, label)])
            aggregate_g13_status_by_frame[label].update(g13_status_by_sample_frame[(key, label)])
            aggregate_oncogenic_status_by_frame[label].update(oncogenic_status_by_sample_frame[(key, label)])
    aggregate_percentages_displayed = {
        label: pct(count, aggregate_denominator) for label, count in aggregate_counts.items()
    }
    aggregate_percentages_indel = {
        label: pct(count, aggregate_indel_denominator) for label, count in aggregate_counts.items()
    }
    aggregate_qc = {
        "plot_sample_count": len(plot_rows),
        "displayed_nonhdr_classifiable_pairs": aggregate_denominator,
        "indel_pairs": aggregate_indel_denominator,
        "indel_pct_of_displayed_nonhdr_classifiable": pct(aggregate_indel_denominator, aggregate_denominator),
        "hdr_excluded_pairs": sum(int(row["hdr_excluded_pairs"]) for row in plot_rows),
        "other_indel_frame_status_pairs": aggregate_counts["other_indel_frame_status"],
        "g12_status_by_frame": {label: dict(aggregate_g12_status_by_frame[label]) for label in FRAME_LABELS + OPTIONAL_FRAME_LABELS},
        "g13_status_by_frame": {label: dict(aggregate_g13_status_by_frame[label]) for label in FRAME_LABELS + OPTIONAL_FRAME_LABELS},
        "known_oncogenic_status_by_frame": {
            label: dict(aggregate_oncogenic_status_by_frame[label]) for label in FRAME_LABELS + OPTIONAL_FRAME_LABELS
        },
    }
    qc_checks = {
        "hdr_denominators_consistent": all(bool(row["hdr_denominator_check_pass"]) for row in rows),
        "indel_counts_equal_frame_split": all(bool(row["indel_count_check_pass"]) for row in rows),
        "aggregate_counts_equal_indel_denominator": sum(aggregate_counts.values()) == aggregate_indel_denominator,
        "plotted_indel_pct_matches_frame_split_pct": abs(
            pct(aggregate_indel_denominator, aggregate_denominator) - sum(aggregate_percentages_displayed.values())
        )
        < 1e-9,
        "g12_status_counts_match_frame_counts": all(
            sum(aggregate_g12_status_by_frame[label].values()) == aggregate_counts[label]
            for label in FRAME_LABELS + OPTIONAL_FRAME_LABELS
        ),
        "g13_status_counts_match_frame_counts": all(
            sum(aggregate_g13_status_by_frame[label].values()) == aggregate_counts[label]
            for label in FRAME_LABELS + OPTIONAL_FRAME_LABELS
        ),
        "known_oncogenic_status_counts_match_frame_counts": all(
            sum(aggregate_oncogenic_status_by_frame[label].values()) == aggregate_counts[label]
            for label in FRAME_LABELS + OPTIONAL_FRAME_LABELS
        ),
    }
    if not all(qc_checks.values()):
        raise AssertionError(f"QC checks failed: {qc_checks}")

    run_info_path = outdir / f"{args.prefix}_plot_run_info.json"
    output_paths = {**outputs, "figure_ready": str(figure_ready), "run_info": str(run_info_path)}
    checksums = {
        "script": sha256_file(Path(__file__).resolve()),
        "helper_script": sha256_file(Path(__file__).resolve().with_name("plot_kras_cutsite_pm20_threebin_background83.py")),
        "inputs": {"per_pair": sha256_file(Path(args.per_pair))},
        "outputs": {name: sha256_file(Path(path)) for name, path in output_paths.items() if name != "run_info"},
    }
    run_info = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).resolve()),
        "argv": sys.argv,
        "window": {"start": cut_start, "end": cut_end, "cut_after": args.cut_site},
        "min_window_baseq": args.min_window_baseq,
        "background_substitution": args.background_substitution,
        "displayed_denominator": "classifiable non-HDR reads",
        "inputs": {"per_pair": args.per_pair},
        "outputs": output_paths,
        "aggregate_counts": aggregate_counts,
        "aggregate_percentages_of_displayed_nonhdr_classifiable": aggregate_percentages_displayed,
        "aggregate_percentages_of_indel_reads": aggregate_percentages_indel,
        "aggregate_qc": aggregate_qc,
        "qc_checks": qc_checks,
        "checksums_sha256": checksums,
        "method_note": (
            "This independent plot splits the three-bin plot's in/del group into frameshift and in-frame "
            "indels using the parser-provided frameshift_status field. The plotted denominator remains "
            "all classifiable non-HDR tumor reads, matching the parent plot's orange in/del segment. "
            "G12/G13 and known-oncogenic statuses are captured in run metadata for annotation planning."
        ),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "matplotlib_version": mpl.__version__,
        "numpy_version": np.__version__,
        "cohort_order": COHORT_ORDER,
    }
    run_info_path.write_text(json.dumps(run_info, indent=2) + "\n")

    print(
        json.dumps(
            {
                "outputs": output_paths,
                "aggregate_counts": aggregate_counts,
                "aggregate_percentages_of_displayed_nonhdr_classifiable": aggregate_percentages_displayed,
                "aggregate_percentages_of_indel_reads": aggregate_percentages_indel,
                "aggregate_qc": aggregate_qc,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
