#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import platform
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


COHORT_ORDER = {"Normal": 0, "KP": 1, "KPTT": 2}
DISPLAY_LABELS = ["WT_background", "in_del", "other_substitution"]
COLORS = {
    "WT_background": "#0072B2",
    "in_del": "#D55E00",
    "other_substitution": "#CC79A7",
}
PRETTY_LABELS = {
    "WT_background": "WT",
    "in_del": "in/del",
    "other_substitution": "other substitutions",
}
PAM_COLOR = "#D55E00"
GUIDE_COLOR = "#6E6E6E"
G12_COLOR = "#009E73"
SUB_RE = re.compile(r"^sub:(\d+):([ACGTN])>([ACGTN])$")
DEL_RE = re.compile(r"^del:(\d+)-(\d+):([ACGTN-]+)$")
INS_RE = re.compile(r"^ins_after:(\d+):([ACGTN]+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot three-bin Kras cleavage-site +/-20 bp allele classes from existing "
            "per-pair cut-site calls, treating a recurrent substitution as background "
            "and excluding HDR reads from the displayed denominator."
        )
    )
    parser.add_argument("--per-pair")
    parser.add_argument("--amplicons")
    parser.add_argument("--outdir")
    parser.add_argument("--prefix", default="kras_cutsite_pm20_threebin_background83")
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
        help="Run synthetic classifier tests and exit without reading input tables.",
    )
    args = parser.parse_args()
    if args.self_test:
        return args
    missing = [name for name in ["per_pair", "amplicons", "outdir"] if getattr(args, name) is None]
    if missing:
        parser.error("Missing required arguments unless --self-test is used: " + ", ".join(f"--{m.replace('_', '-')}" for m in missing))
    return args


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7.5,
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
        }
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_amplicon(path: Path) -> dict[str, str]:
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["locus"] == "Kras":
                return row
    raise ValueError(f"Kras row not found in {path}")


def parse_int_field(row: dict[str, str], field: str) -> int | None:
    value = row.get(field, "")
    if value in {"", "NA"}:
        return None
    return int(value)


def is_full_window_classifiable(
    row: dict[str, str],
    start: int,
    end: int,
    min_window_baseq: int,
) -> bool:
    if row.get("allele_class") == "ambiguous_low_confidence":
        return False
    if row.get("alignment_status") != "ok":
        return False
    ref_start = parse_int_field(row, "ref_start_1based")
    ref_end = parse_int_field(row, "ref_end_1based")
    minq = parse_int_field(row, "cut_window_min_baseq")
    if ref_start is None or ref_end is None or minq is None:
        return False
    return ref_start <= start and ref_end >= end and minq >= min_window_baseq


def is_hdr_row(row: dict[str, str]) -> bool:
    return row["paired_hdr_status"] == "hdr_exact" and row["allele_class"] == "hdr_kras_g12d"


def event_flags(
    event_key: str,
    start: int,
    end: int,
    background_substitution: str,
) -> dict[str, int | bool]:
    flags: dict[str, int | bool] = {
        "has_window_indel": False,
        "has_other_window_substitution": False,
        "has_background_substitution": False,
        "background_substitution_event_count": 0,
    }
    if event_key in {"", "NA", "WT_observed_region"}:
        return flags

    for event in event_key.split(";"):
        if event == background_substitution:
            flags["has_background_substitution"] = True
            flags["background_substitution_event_count"] = int(flags["background_substitution_event_count"]) + 1
            continue

        sub = SUB_RE.match(event)
        if sub:
            pos = int(sub.group(1))
            if start <= pos <= end:
                flags["has_other_window_substitution"] = True
            continue

        deletion = DEL_RE.match(event)
        if deletion:
            dstart = int(deletion.group(1))
            dend = int(deletion.group(2))
            if not (dend < start or dstart > end):
                flags["has_window_indel"] = True
            continue

        insertion = INS_RE.match(event)
        if insertion:
            anchor = int(insertion.group(1))
            if start <= anchor <= end:
                flags["has_window_indel"] = True
            continue

        raise ValueError(f"Unrecognized event token: {event!r}")

    return flags


def classify_nonhdr_row(
    row: dict[str, str],
    start: int,
    end: int,
    background_substitution: str,
) -> tuple[str, dict[str, int | bool]]:
    flags = event_flags(row.get("full_event_key", ""), start, end, background_substitution)
    if bool(flags["has_window_indel"]):
        return "in_del", flags
    if bool(flags["has_other_window_substitution"]):
        return "other_substitution", flags
    return "WT_background", flags


def pct(n: int, d: int) -> float:
    return 100.0 * n / d if d else 0.0


def sort_key(row: dict[str, object]) -> tuple[int, int, str]:
    specimen = str(row["specimen_id"])
    suffix = specimen.split("_")[-1]
    try:
        numeric = int(suffix)
    except ValueError:
        numeric = 10**9
    return (COHORT_ORDER.get(str(row["cohort"]), 99), numeric, specimen)


def panel_label(ax: plt.Axes, label: str, outside: bool = False) -> None:
    if outside:
        ax.text(
            -0.055,
            1.02,
            label,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=10,
            fontweight="bold",
            clip_on=False,
        )
        return
    ax.text(
        0.005,
        0.98,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        fontweight="bold",
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.0, "alpha": 0.9},
    )


def draw_window_sequence(ax: plt.Axes, wt_seq: str, start: int, end: int, cut_site: int) -> None:
    ax.text(1, 0.07, f"WT bases {start}-{end}", ha="left", va="center", fontsize=6.0)
    for pos in range(start, end + 1):
        base = wt_seq[pos - 1]
        color = PAM_COLOR if 86 <= pos <= 88 else "black"
        weight = "bold" if pos in {cut_site, cut_site + 1} else "normal"
        ax.text(
            pos,
            0.07,
            base,
            ha="center",
            va="center",
            fontsize=5.2,
            family="monospace",
            color=color,
            fontweight=weight,
        )
    ax.text(
        cut_site + 0.5,
        0.07,
        "|",
        ha="center",
        va="center",
        fontsize=5.8,
        family="monospace",
        fontweight="bold",
        color="black",
    )


def draw_scheme(ax: plt.Axes, wt_seq: str, cut_site: int, start: int, end: int) -> None:
    ax.set_xlim(1, len(wt_seq))
    ax.set_ylim(0, 1)
    ax.axis("off")
    y = 0.5
    ax.hlines(y, 1, len(wt_seq), color="black", linewidth=1.2)
    ax.axvspan(start, end, ymin=0.38, ymax=0.62, color="#56B4E9", alpha=0.25)
    ax.vlines(cut_site + 0.5, y - 0.33, y + 0.33, color="black", linestyle="--", linewidth=1.0)
    ax.text(start, 0.87, f"cleavage site +/-20 bp window (bases {start}-{end})", ha="left", va="center")
    ax.text(cut_site + 0.5, 0.16, f"cut after {cut_site}", ha="center", va="center", fontsize=7)
    ax.text(22, 0.74, "G12", ha="center", va="center", fontsize=7, color=G12_COLOR)
    ax.text(86.5, 0.74, "PAM", ha="center", va="center", fontsize=7, color=PAM_COLOR)
    ax.text(89.5, 0.67, "guide reverse-complement", ha="left", va="center", fontsize=7, color=GUIDE_COLOR)
    ax.axvspan(21, 23, ymin=0.62, ymax=0.70, color=G12_COLOR, alpha=0.65)
    ax.axvspan(86, 88, ymin=0.62, ymax=0.70, color=PAM_COLOR, alpha=0.65)
    ax.axvspan(89, 108, ymin=0.52, ymax=0.60, color=GUIDE_COLOR, alpha=0.45)
    ax.text(22, 0.30, "21-23", ha="center", va="center", fontsize=6)
    for pos in [1, start, 86, 88, 91, 92, end, len(wt_seq)]:
        ax.text(pos, 0.30, str(pos), ha="center", va="center", fontsize=6)
    draw_window_sequence(ax, wt_seq, start, end, cut_site)


def add_cohort_separators(ax: plt.Axes, rows: list[dict[str, object]]) -> None:
    cohorts = [str(row["cohort"]) for row in rows]
    start = 0
    for idx in range(1, len(cohorts) + 1):
        if idx == len(cohorts) or cohorts[idx] != cohorts[start]:
            end = idx - 1
            if end < len(cohorts) - 1:
                ax.axvline(end + 0.5, color="#999999", linewidth=0.6)
            start = idx


def run_self_tests() -> None:
    base_row = {
        "paired_hdr_status": "other_pattern",
        "allele_class": "substitution_only_cutsite_window",
        "full_event_key": "WT_observed_region",
    }

    def classify(event_key: str) -> str:
        row = {**base_row, "full_event_key": event_key}
        label, _ = classify_nonhdr_row(row, 71, 111, "sub:83:T>C")
        return label

    tests = {
        "background_only": (classify("sub:83:T>C"), "WT_background"),
        "background_plus_other_sub": (classify("sub:83:T>C;sub:93:A>T"), "other_substitution"),
        "indel_plus_subs": (classify("sub:83:T>C;del:89-89:T;sub:93:A>T"), "in_del"),
        "insertion_plus_subs": (classify("sub:83:T>C;ins_after:90:G;sub:93:A>T"), "in_del"),
        "non_background_pos83": (classify("sub:83:T>G"), "other_substitution"),
    }
    failures = {name: {"observed": obs, "expected": exp} for name, (obs, exp) in tests.items() if obs != exp}
    hdr_test = is_hdr_row({"paired_hdr_status": "hdr_exact", "allele_class": "hdr_kras_g12d"})
    if not hdr_test:
        failures["hdr_exclusion"] = {"observed": hdr_test, "expected": True}
    if failures:
        raise AssertionError(json.dumps(failures, indent=2))
    print("Synthetic classifier tests passed.", file=sys.stderr)


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

    amp = read_amplicon(Path(args.amplicons))
    wt_seq = amp["wt_amplicon_seq"].upper()

    counts: dict[tuple[str, str, str, str], int] = Counter()
    raw_by_sample: dict[tuple[str, str, str], int] = Counter()
    full_classifiable_by_sample: dict[tuple[str, str, str], int] = Counter()
    display_denominator_by_sample: dict[tuple[str, str, str], int] = Counter()
    qc_by_sample: dict[tuple[str, str, str], Counter[str]] = defaultdict(Counter)

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
            call, flags = classify_nonhdr_row(row, cut_start, cut_end, args.background_substitution)
            counts[(*key, call)] += 1

            if bool(flags["has_background_substitution"]):
                qc_by_sample[key]["background_substitution_nonhdr_classifiable_pairs"] += 1
                if bool(flags["has_window_indel"]):
                    qc_by_sample[key]["background_substitution_in_indel_pairs"] += 1
                elif bool(flags["has_other_window_substitution"]):
                    qc_by_sample[key]["background_substitution_plus_other_sub_no_indel_pairs"] += 1
                else:
                    qc_by_sample[key]["moved_other_substitution_to_WT_background_pairs"] += 1

    rows: list[dict[str, object]] = []
    for key in sorted(raw_by_sample, key=lambda x: sort_key({"specimen_id": x[0], "cohort": x[1]})):
        specimen, cohort, sample_name = key
        denominator = display_denominator_by_sample[key]
        hdr_excluded = qc_by_sample[key]["hdr_excluded_pairs"]
        full_classifiable = full_classifiable_by_sample[key]
        row: dict[str, object] = {
            "specimen_id": specimen,
            "cohort": cohort,
            "sample_name": sample_name,
            "plot_include": "yes" if (not args.tumor_only or cohort != "Normal") else "no",
            "raw_pairs": raw_by_sample[key],
            "not_full_window_classifiable_pairs": qc_by_sample[key]["not_full_window_classifiable_pairs"],
            "full_pm20_classifiable_pairs_including_HDR": full_classifiable,
            "hdr_excluded_pairs": hdr_excluded,
            "displayed_nonhdr_classifiable_pairs": denominator,
            "pm20_window_start": cut_start,
            "pm20_window_end": cut_end,
            "min_window_baseq": args.min_window_baseq,
            "background_substitution": args.background_substitution,
            "background_substitution_nonhdr_classifiable_pairs": qc_by_sample[key][
                "background_substitution_nonhdr_classifiable_pairs"
            ],
            "background_substitution_nonhdr_classifiable_pct": pct(
                qc_by_sample[key]["background_substitution_nonhdr_classifiable_pairs"], denominator
            ),
            "moved_other_substitution_to_WT_background_pairs": qc_by_sample[key][
                "moved_other_substitution_to_WT_background_pairs"
            ],
            "moved_other_substitution_to_WT_background_pct": pct(
                qc_by_sample[key]["moved_other_substitution_to_WT_background_pairs"], denominator
            ),
            "background_substitution_plus_other_sub_no_indel_pairs": qc_by_sample[key][
                "background_substitution_plus_other_sub_no_indel_pairs"
            ],
            "background_substitution_in_indel_pairs": qc_by_sample[key]["background_substitution_in_indel_pairs"],
        }
        pct_sum = 0.0
        for label in DISPLAY_LABELS:
            n = counts[(*key, label)]
            row[f"{label}_pairs"] = n
            row[f"{label}_pct"] = pct(n, denominator)
            pct_sum += float(row[f"{label}_pct"])
        row["three_bin_pct_sum"] = pct_sum
        row["hdr_denominator_check_pass"] = (full_classifiable - hdr_excluded) == denominator
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
    zero_denominator = [row["specimen_id"] for row in plot_rows if int(row["displayed_nonhdr_classifiable_pairs"]) == 0]
    if zero_denominator:
        raise ValueError(f"Displayed non-HDR denominator is zero for plotted samples: {zero_denominator}")

    fig = plt.figure(figsize=(7.2, 5.7))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 2.4], hspace=0.55)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    draw_scheme(ax_a, wt_seq, args.cut_site, cut_start, cut_end)
    panel_label(ax_a, "A")

    x = np.arange(len(plot_rows))
    bottom = np.zeros(len(plot_rows))
    for label in DISPLAY_LABELS:
        vals = np.array([float(row[f"{label}_pct"]) for row in plot_rows])
        ax_b.bar(
            x,
            vals,
            bottom=bottom,
            color=COLORS[label],
            edgecolor="white",
            linewidth=0.25,
            label=PRETTY_LABELS[label],
        )
        bottom += vals
    ax_b.set_ylim(0, 100)
    ax_b.set_ylabel("Classifiable non-HDR reads (%)")
    ax_b.set_xticks(x)
    ax_b.set_xticklabels([str(row["specimen_id"]) for row in plot_rows], rotation=60, ha="right")
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    ax_b.grid(axis="y", color="#E5E5E5", linewidth=0.5)
    ax_b.set_axisbelow(True)
    ax_b.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.16))
    add_cohort_separators(ax_b, plot_rows)
    panel_label(ax_b, "B", outside=True)

    fig.suptitle(
        "KRAS cut-site +/-20 bp background-adjusted allele classes",
        y=0.985,
        fontsize=10,
        fontweight="bold",
    )
    outputs: dict[str, str] = {}
    for ext in ["pdf", "svg", "png"]:
        path = outdir / f"{args.prefix}.{ext}"
        fig.savefig(path, dpi=600 if ext == "png" else None, bbox_inches="tight")
        outputs[f"figure_{ext}"] = str(path)
    plt.close(fig)

    aggregate_counts = {label: sum(int(row[f"{label}_pairs"]) for row in plot_rows) for label in DISPLAY_LABELS}
    aggregate_denominator = sum(int(row["displayed_nonhdr_classifiable_pairs"]) for row in plot_rows)
    aggregate_qc = {
        "plot_sample_count": len(plot_rows),
        "displayed_nonhdr_classifiable_pairs": aggregate_denominator,
        "full_pm20_classifiable_pairs_including_HDR": sum(
            int(row["full_pm20_classifiable_pairs_including_HDR"]) for row in plot_rows
        ),
        "hdr_excluded_pairs": sum(int(row["hdr_excluded_pairs"]) for row in plot_rows),
        "background_substitution_nonhdr_classifiable_pairs": sum(
            int(row["background_substitution_nonhdr_classifiable_pairs"]) for row in plot_rows
        ),
        "moved_other_substitution_to_WT_background_pairs": sum(
            int(row["moved_other_substitution_to_WT_background_pairs"]) for row in plot_rows
        ),
        "background_substitution_plus_other_sub_no_indel_pairs": sum(
            int(row["background_substitution_plus_other_sub_no_indel_pairs"]) for row in plot_rows
        ),
        "background_substitution_in_indel_pairs": sum(
            int(row["background_substitution_in_indel_pairs"]) for row in plot_rows
        ),
    }
    aggregate_percentages = {label: pct(count, aggregate_denominator) for label, count in aggregate_counts.items()}
    aggregate_qc.update(
        {
            "background_substitution_nonhdr_classifiable_pct": pct(
                int(aggregate_qc["background_substitution_nonhdr_classifiable_pairs"]), aggregate_denominator
            ),
            "moved_other_substitution_to_WT_background_pct": pct(
                int(aggregate_qc["moved_other_substitution_to_WT_background_pairs"]), aggregate_denominator
            ),
            "three_bin_pct_sum": sum(aggregate_percentages.values()),
        }
    )
    qc_checks = {
        "plot_rows_pct_sum_close_to_100": all(abs(float(row["three_bin_pct_sum"]) - 100.0) < 1e-9 for row in plot_rows),
        "hdr_denominators_consistent": all(bool(row["hdr_denominator_check_pass"]) for row in rows),
        "aggregate_counts_equal_displayed_denominator": sum(aggregate_counts.values()) == aggregate_denominator,
    }
    if not all(qc_checks.values()):
        raise AssertionError(f"QC checks failed: {qc_checks}")

    run_info_path = outdir / f"{args.prefix}_plot_run_info.json"
    output_paths = {**outputs, "figure_ready": str(figure_ready), "run_info": str(run_info_path)}
    checksums = {
        "script": sha256_file(Path(__file__).resolve()),
        "inputs": {
            "per_pair": sha256_file(Path(args.per_pair)),
            "amplicons": sha256_file(Path(args.amplicons)),
        },
        "outputs": {name: sha256_file(Path(path)) for name, path in output_paths.items() if name != "run_info"},
    }
    run_info = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).resolve()),
        "argv": sys.argv,
        "window": {"start": cut_start, "end": cut_end, "cut_after": args.cut_site},
        "min_window_baseq": args.min_window_baseq,
        "background_substitution": args.background_substitution,
        "displayed_bins": PRETTY_LABELS,
        "inputs": {"per_pair": args.per_pair, "amplicons": args.amplicons},
        "outputs": output_paths,
        "aggregate_counts": aggregate_counts,
        "aggregate_percentages": aggregate_percentages,
        "aggregate_qc": aggregate_qc,
        "qc_checks": qc_checks,
        "checksums_sha256": checksums,
        "method_note": (
            "Uses per-pair cut-site alignments from the existing +/-20 bp parser output; includes only "
            "non-ambiguous aligned reads spanning bases 71-111 with Q20 over the window. HDR reads are "
            "excluded from the displayed denominator. For non-HDR reads, exact sub:83:T>C is treated as "
            "background for display classification; any overlapping indel takes priority over substitutions."
        ),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "matplotlib_version": mpl.__version__,
        "numpy_version": np.__version__,
    }
    run_info_path.write_text(json.dumps(run_info, indent=2) + "\n")

    print(json.dumps({"outputs": output_paths, "aggregate_counts": aggregate_counts, "aggregate_qc": aggregate_qc}, indent=2))


if __name__ == "__main__":
    main()
