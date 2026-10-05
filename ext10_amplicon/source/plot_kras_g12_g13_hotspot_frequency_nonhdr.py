#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
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


COHORT_ORDER = {"Normal": 0, "KP": 1, "KPTT": 2}
CATEGORY_ORDER = [
    "WT",
    "G12D",
    "G12V",
    "G12C",
    "G12A",
    "G12R",
    "G12S",
    "G13D",
    "G13C",
    "G13R",
    "other",
]
SELECTED_VARIANT_CATEGORIES = [
    category for category in CATEGORY_ORDER if category not in {"WT", "other"}
]
COLORS = {
    "WT": "#D9D9D9",
    "G12D": "#0072B2",
    "G12V": "#56B4E9",
    "G12C": "#009E73",
    "G12A": "#F0E442",
    "G12R": "#E69F00",
    "G12S": "#D55E00",
    "G13D": "#CC79A7",
    "G13C": "#882255",
    "G13R": "#44AA99",
    "other": "#666666",
}
REQUESTED_G12 = {"D", "V", "C", "A", "R", "S"}
REQUESTED_G13 = {"D", "C", "R"}

EXPECTED_OLD_TUMOR_DENOMINATOR = 19_490_107
EXPECTED_OLD_G12D_READS = 15_269
EXPECTED_HDR_POSITIVE_REMOVED_READS = 2_594
EXPECTED_HDR_POSITIVE_G12D_READS = 2_589
EXPECTED_CORRECTED_DENOMINATOR = 19_487_513
EXPECTED_CORRECTED_G12D_READS = 12_680
EXPECTED_CORRECTED_COUNTS = {
    "WT": 18_730_422,
    "G12D": 12_680,
    "G12V": 50_182,
    "G12C": 34_035,
    "G12A": 7_091,
    "G12R": 9_014,
    "G12S": 18_370,
    "G13D": 13_814,
    "G13C": 8_684,
    "G13R": 6_153,
    "other": 597_068,
}
FIGURE_TITLE = (
    "Frequencies of point mutations affecting Kras codons 12 and 13 among "
    "non-HDR alleles"
)
DENOMINATOR_TYPE = "non-HDR tumor-only allele reads from codon13-window spectrum"
REQUIRED_COLUMNS = {
    "specimen_id",
    "cohort",
    "read_count",
    "strict_hdr_barcode",
    "window_class",
    "codon12_effect",
    "codon12_aa",
    "codon13_effect",
    "codon13_aa",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate corrected Supplementary Figure 5D KRAS G12/G13 point-mutation "
            "frequencies after excluding HDR-barcode-positive codon-window alleles."
        )
    )
    parser.add_argument("--alleles", required=False)
    parser.add_argument("--outdir", required=False)
    parser.add_argument(
        "--prefix",
        default="kras_g12_g13_hotspot_frequency_tumor_only_nonhdr",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run synthetic category tests and exit without reading input tables.",
    )
    args = parser.parse_args()
    if args.self_test:
        return args
    missing = [name for name in ["alleles", "outdir"] if getattr(args, name) is None]
    if missing:
        parser.error(
            "Missing required arguments unless --self-test is used: "
            + ", ".join(f"--{name}" for name in missing)
        )
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
            "legend.fontsize": 6.2,
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


def pct(numerator: int | float, denominator: int | float) -> float:
    return 100.0 * numerator / denominator if denominator else 0.0


def sort_key(row: dict[str, object]) -> tuple[int, int, str]:
    specimen = str(row["specimen_id"])
    suffix = specimen.split("_")[-1]
    try:
        numeric = int(suffix)
    except ValueError:
        numeric = 10**9
    return (COHORT_ORDER.get(str(row["cohort"]), 99), numeric, specimen)


def requested_g12_g13_hits(row: dict[str, str]) -> list[str]:
    hits: list[str] = []
    if (
        row["codon12_effect"] == "known_oncogenic_residue_change"
        and row["codon12_aa"] in REQUESTED_G12
    ):
        hits.append(f"G12{row['codon12_aa']}")
    if (
        row["codon13_effect"] == "known_oncogenic_residue_change"
        and row["codon13_aa"] in REQUESTED_G13
    ):
        hits.append(f"G13{row['codon13_aa']}")
    return hits


def classify_allele_row(row: dict[str, str]) -> tuple[str, list[str]]:
    hits = requested_g12_g13_hits(row)
    if (
        not hits
        and row["codon12_effect"] == "exact_wt_codon"
        and row["codon13_effect"] == "exact_wt_codon"
    ):
        return "WT", hits
    if len(hits) == 1:
        return hits[0], hits
    return "other", hits


def run_self_tests() -> None:
    base = {
        "codon12_effect": "exact_wt_codon",
        "codon12_aa": "G",
        "codon13_effect": "exact_wt_codon",
        "codon13_aa": "G",
    }
    tests = {
        "WT": (classify_allele_row(base)[0], "WT"),
        "G12D": (
            classify_allele_row(
                {
                    **base,
                    "codon12_effect": "known_oncogenic_residue_change",
                    "codon12_aa": "D",
                }
            )[0],
            "G12D",
        ),
        "G13R": (
            classify_allele_row(
                {
                    **base,
                    "codon13_effect": "known_oncogenic_residue_change",
                    "codon13_aa": "R",
                }
            )[0],
            "G13R",
        ),
        "compound_selected_to_other": (
            classify_allele_row(
                {
                    **base,
                    "codon12_effect": "known_oncogenic_residue_change",
                    "codon12_aa": "V",
                    "codon13_effect": "known_oncogenic_residue_change",
                    "codon13_aa": "C",
                }
            )[0],
            "other",
        ),
        "unrequested_residue_to_other": (
            classify_allele_row(
                {
                    **base,
                    "codon13_effect": "known_oncogenic_residue_change",
                    "codon13_aa": "V",
                }
            )[0],
            "other",
        ),
        "uncallable_to_other": (
            classify_allele_row(
                {
                    **base,
                    "codon12_effect": "uncallable_indel_or_ambiguous",
                    "codon12_aa": "NA",
                }
            )[0],
            "other",
        ),
    }
    failures = {
        name: {"observed": observed, "expected": expected}
        for name, (observed, expected) in tests.items()
        if observed != expected
    }
    if failures:
        raise AssertionError(json.dumps(failures, indent=2))
    print("Synthetic G12/G13 non-HDR category tests passed.", file=sys.stderr)


def validate_header(fieldnames: list[str] | None) -> None:
    if not fieldnames:
        raise ValueError("Input allele TSV has no header.")
    missing = sorted(REQUIRED_COLUMNS - set(fieldnames))
    if missing:
        raise ValueError(f"Input allele TSV is missing required columns: {missing}")


def add_cohort_separators(ax: plt.Axes, rows: list[dict[str, object]]) -> None:
    cohorts = [str(row["cohort"]) for row in rows]
    start = 0
    for idx in range(1, len(cohorts) + 1):
        if idx == len(cohorts) or cohorts[idx] != cohorts[start]:
            end = idx - 1
            if end < len(cohorts) - 1:
                ax.axvline(end + 0.5, color="#999999", linewidth=0.6)
            start = idx


def stacked_bar(
    ax: plt.Axes,
    rows: list[dict[str, object]],
    categories: list[str],
    *,
    ylabel: str,
    ylim: float,
) -> None:
    x = np.arange(len(rows))
    bottom = np.zeros(len(rows))
    for category in categories:
        values = np.array([float(row[f"{category}_pct"]) for row in rows])
        ax.bar(
            x,
            values,
            bottom=bottom,
            color=COLORS[category],
            edgecolor="white",
            linewidth=0.25,
            label=category,
        )
        bottom += values
    ax.set_ylim(0, ylim)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels([str(row["specimen_id"]) for row in rows], rotation=60, ha="right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.5)
    ax.set_axisbelow(True)
    add_cohort_separators(ax, rows)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows available for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_readme(path: Path, prefix: str, run_info_path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "# Corrected SF5D KRAS G12/G13 Non-HDR Panel",
                "",
                "This directory contains the corrected Supplementary Figure 5D point-mutation panel for Kras codons 12 and 13.",
                "",
                "## Correction",
                "",
                "- The previous bundled G12/G13 hotspot panel used all tumor-only codon13-window allele reads.",
                "- HDR-barcode-positive alleles are excluded here with `strict_hdr_barcode == no`.",
                "- The plotted denominator is non-HDR tumor-only allele reads.",
                "- G12D is retained only when it is non-HDR/unintended G12D after this filtering.",
                "",
                "## Source",
                "",
                "- Source allele table: `remote_kras_codon13_window_spectrum_copy/results/kras_codon13_window_spectrum/per_allele/kras_codon13_window_allele_counts.tsv`",
                "- Tumor filter: `cohort != Normal`",
                "- HDR exclusion: `strict_hdr_barcode == yes` rows are excluded from numerator and denominator.",
                "",
                "## Outputs",
                "",
                f"- `figures/{prefix}.png`",
                f"- `figures/{prefix}.pdf`",
                f"- `figures/{prefix}.svg`",
                f"- `data/{prefix}_figure_ready.tsv`",
                f"- `data/{prefix}_data.xlsx`",
                f"- `{run_info_path.relative_to(path.parent)}`",
                "",
                "## Expected aggregate QC",
                "",
                "- Old tumor-only denominator: 19,490,107 allele reads.",
                "- HDR-positive tumor alleles removed: 2,594 reads.",
                "- Corrected non-HDR tumor-only denominator: 19,487,513 allele reads.",
                "- Old bundled G12D: 15,269 reads / 0.078342%.",
                "- HDR-positive G12D removed: 2,589 reads.",
                "- Corrected non-HDR G12D: 12,680 reads / 0.065067%.",
                "",
            ]
        )
    )


def fail_if_expected_values_drift(
    *,
    old_denominator: int,
    old_counts: Counter[str],
    hdr_positive_denominator: int,
    hdr_positive_counts: Counter[str],
    corrected_denominator: int,
    corrected_counts: Counter[str],
) -> None:
    observed = {
        "old_tumor_denominator": old_denominator,
        "old_g12d_reads": old_counts["G12D"],
        "hdr_positive_removed_reads": hdr_positive_denominator,
        "hdr_positive_g12d_reads": hdr_positive_counts["G12D"],
        "corrected_denominator": corrected_denominator,
        "corrected_g12d_reads": corrected_counts["G12D"],
        "corrected_counts": {category: corrected_counts[category] for category in CATEGORY_ORDER},
    }
    expected = {
        "old_tumor_denominator": EXPECTED_OLD_TUMOR_DENOMINATOR,
        "old_g12d_reads": EXPECTED_OLD_G12D_READS,
        "hdr_positive_removed_reads": EXPECTED_HDR_POSITIVE_REMOVED_READS,
        "hdr_positive_g12d_reads": EXPECTED_HDR_POSITIVE_G12D_READS,
        "corrected_denominator": EXPECTED_CORRECTED_DENOMINATOR,
        "corrected_g12d_reads": EXPECTED_CORRECTED_G12D_READS,
        "corrected_counts": EXPECTED_CORRECTED_COUNTS,
    }
    if observed != expected:
        raise AssertionError(
            "Corrected SF5D aggregate values do not match expected locked values:\n"
            + json.dumps({"observed": observed, "expected": expected}, indent=2)
        )


def main() -> None:
    args = parse_args()
    if args.self_test:
        run_self_tests()
        return

    configure_style()
    allele_path = Path(args.alleles).resolve()
    package_dir = Path(args.outdir).resolve()
    figure_dir = package_dir / "figures"
    data_dir = package_dir / "data"
    metadata_dir = package_dir / "metadata"
    for directory in [figure_dir, data_dir, metadata_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    counts_by_sample: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    denominator_by_sample: Counter[tuple[str, str]] = Counter()
    old_counts_by_sample: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    old_denominator_by_sample: Counter[tuple[str, str]] = Counter()
    hdr_positive_counts_by_sample: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    hdr_positive_denominator_by_sample: Counter[tuple[str, str]] = Counter()
    compound_by_sample: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    window_class_by_sample: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    effect_pair_by_sample: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    included_strict_hdr_status: Counter[str] = Counter()
    source_accounting = Counter()

    with allele_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_header(reader.fieldnames)
        for row in reader:
            read_count = int(row["read_count"])
            category, hits = classify_allele_row(row)

            if row["cohort"] == "Normal":
                source_accounting["normal_rows_excluded"] += 1
                source_accounting["normal_reads_excluded"] += read_count
                continue

            key = (row["specimen_id"], row["cohort"])
            old_counts_by_sample[key][category] += read_count
            old_denominator_by_sample[key] += read_count
            source_accounting["tumor_rows_seen"] += 1
            source_accounting["tumor_reads_seen"] += read_count

            if row["strict_hdr_barcode"] == "yes":
                hdr_positive_counts_by_sample[key][category] += read_count
                hdr_positive_denominator_by_sample[key] += read_count
                source_accounting["tumor_strict_hdr_yes_rows_excluded"] += 1
                source_accounting["tumor_strict_hdr_yes_reads_excluded"] += read_count
                continue
            if row["strict_hdr_barcode"] != "no":
                raise ValueError(
                    "Unexpected strict_hdr_barcode value "
                    f"{row['strict_hdr_barcode']!r} for {row['specimen_id']}"
                )

            included_strict_hdr_status[row["strict_hdr_barcode"]] += read_count
            counts_by_sample[key][category] += read_count
            denominator_by_sample[key] += read_count
            source_accounting["tumor_nonhdr_rows_included"] += 1
            source_accounting["tumor_nonhdr_reads_included"] += read_count
            if len(hits) > 1:
                compound_by_sample[key][";".join(hits)] += read_count
            window_class_by_sample[key][row["window_class"]] += read_count
            effect_pair = f"{row['codon12_effect']}|{row['codon13_effect']}"
            effect_pair_by_sample[key][effect_pair] += read_count

    if not denominator_by_sample:
        raise ValueError("No tumor non-HDR allele rows found.")

    rows: list[dict[str, object]] = []
    for key in sorted(
        denominator_by_sample,
        key=lambda item: sort_key({"specimen_id": item[0], "cohort": item[1]}),
    ):
        specimen, cohort = key
        denominator = denominator_by_sample[key]
        old_denominator = old_denominator_by_sample[key]
        hdr_removed = hdr_positive_denominator_by_sample[key]
        row: dict[str, object] = {
            "specimen_id": specimen,
            "cohort": cohort,
            "plot_include": "yes",
            "denominator_type": DENOMINATOR_TYPE,
            "full_tumor_allele_reads_before_hdr_filter": old_denominator,
            "strict_hdr_barcode_yes_reads_excluded": hdr_removed,
            "denominator_nonhdr_tumor_allele_reads": denominator,
            "strict_hdr_barcode_yes_reads_contributed_to_denominator": 0,
            "compound_selected_g12_g13_reads": sum(compound_by_sample[key].values()),
        }
        pct_sum = 0.0
        for category in CATEGORY_ORDER:
            count = counts_by_sample[key][category]
            row[f"{category}_reads"] = count
            row[f"{category}_pct"] = pct(count, denominator)
            pct_sum += float(row[f"{category}_pct"])
        row["category_pct_sum"] = pct_sum
        row["category_count_check_pass"] = (
            sum(int(row[f"{category}_reads"]) for category in CATEGORY_ORDER)
            == denominator
        )
        rows.append(row)

    figure_ready = data_dir / f"{args.prefix}_figure_ready.tsv"
    write_tsv(figure_ready, rows)

    non_wt_categories = [category for category in CATEGORY_ORDER if category != "WT"]
    non_wt_max = max(
        sum(float(row[f"{category}_pct"]) for category in non_wt_categories)
        for row in rows
    )
    zoom_ylim = max(5.0, float(np.ceil((non_wt_max + 1.0) / 2.5) * 2.5))

    fig, (ax_full, ax_zoom) = plt.subplots(
        2,
        1,
        figsize=(8.2, 6.8),
        gridspec_kw={"height_ratios": [1.25, 1.45]},
    )
    fig.subplots_adjust(left=0.09, right=0.995, bottom=0.16, top=0.79, hspace=0.58)

    stacked_bar(
        ax_full,
        rows,
        CATEGORY_ORDER,
        ylabel="Non-HDR allele reads (%)",
        ylim=100.0,
    )
    ax_full.set_title("All categories", fontsize=8, pad=5)
    ax_full.legend(
        frameon=False,
        ncol=6,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.58),
        columnspacing=1.1,
        handlelength=1.5,
    )

    stacked_bar(
        ax_zoom,
        rows,
        non_wt_categories,
        ylabel="Non-HDR allele reads (%)",
        ylim=zoom_ylim,
    )
    ax_zoom.set_title("Non-WT zoom, same non-HDR denominator", fontsize=8, pad=5)
    fig.suptitle(FIGURE_TITLE, y=0.965, fontsize=10, fontweight="bold")

    output_paths: dict[str, str] = {}
    for ext in ["pdf", "svg", "png"]:
        figure_path = figure_dir / f"{args.prefix}.{ext}"
        fig.savefig(figure_path, dpi=600 if ext == "png" else None, bbox_inches="tight")
        output_paths[f"figure_{ext}"] = str(figure_path)
    plt.close(fig)

    aggregate_counts = {
        category: sum(int(row[f"{category}_reads"]) for row in rows)
        for category in CATEGORY_ORDER
    }
    aggregate_denominator = sum(
        int(row["denominator_nonhdr_tumor_allele_reads"]) for row in rows
    )
    aggregate_percentages = {
        category: pct(count, aggregate_denominator)
        for category, count in aggregate_counts.items()
    }
    old_aggregate_counts = Counter(
        {
            category: sum(
                old_counts_by_sample[key][category] for key in old_denominator_by_sample
            )
            for category in CATEGORY_ORDER
        }
    )
    hdr_positive_counts = Counter(
        {
            category: sum(
                hdr_positive_counts_by_sample[key][category]
                for key in hdr_positive_denominator_by_sample
            )
            for category in CATEGORY_ORDER
        }
    )
    old_aggregate_denominator = sum(old_denominator_by_sample.values())
    hdr_positive_denominator = sum(hdr_positive_denominator_by_sample.values())

    fail_if_expected_values_drift(
        old_denominator=old_aggregate_denominator,
        old_counts=old_aggregate_counts,
        hdr_positive_denominator=hdr_positive_denominator,
        hdr_positive_counts=hdr_positive_counts,
        corrected_denominator=aggregate_denominator,
        corrected_counts=Counter(aggregate_counts),
    )

    aggregate_compound: Counter[str] = Counter()
    aggregate_window_classes: Counter[str] = Counter()
    aggregate_effect_pairs: Counter[str] = Counter()
    for row in rows:
        key = (str(row["specimen_id"]), str(row["cohort"]))
        aggregate_compound.update(compound_by_sample[key])
        aggregate_window_classes.update(window_class_by_sample[key])
        aggregate_effect_pairs.update(effect_pair_by_sample[key])

    qc_checks = {
        "no_strict_hdr_barcode_yes_reads_contribute_to_corrected_sf5d": (
            included_strict_hdr_status.get("yes", 0) == 0
        ),
        "normal_samples_excluded_from_corrected_sf5d": all(
            row["cohort"] != "Normal" for row in rows
        ),
        "normal_458_and_normal_460_absent_from_corrected_sf5d": not any(
            row["specimen_id"] in {"Normal_458", "Normal_460"} for row in rows
        ),
        "per_sample_counts_sum_to_nonhdr_denominator": all(
            bool(row["category_count_check_pass"]) for row in rows
        ),
        "per_sample_pct_sum_close_to_100": all(
            abs(float(row["category_pct_sum"]) - 100.0) < 1e-9 for row in rows
        ),
        "aggregate_counts_sum_to_nonhdr_denominator": (
            sum(aggregate_counts.values()) == aggregate_denominator
        ),
        "aggregate_denominator_matches_expected": (
            aggregate_denominator == EXPECTED_CORRECTED_DENOMINATOR
        ),
        "old_denominator_minus_hdr_positive_equals_corrected_denominator": (
            old_aggregate_denominator - hdr_positive_denominator == aggregate_denominator
        ),
        "old_g12d_minus_hdr_positive_g12d_equals_corrected_g12d": (
            old_aggregate_counts["G12D"] - hdr_positive_counts["G12D"]
            == aggregate_counts["G12D"]
        ),
    }
    if not all(qc_checks.values()):
        raise AssertionError(f"QC checks failed: {json.dumps(qc_checks, indent=2)}")

    readme_path = package_dir / "README.md"
    run_info_path = metadata_dir / f"{args.prefix}_run_info.json"
    write_readme(readme_path, args.prefix, run_info_path)

    output_paths.update(
        {
            "figure_ready_tsv": str(figure_ready),
            "readme": str(readme_path),
            "run_info": str(run_info_path),
        }
    )
    checksums = {
        "script": sha256_file(Path(__file__).resolve()),
        "inputs": {"alleles": sha256_file(allele_path)},
        "outputs": {
            name: sha256_file(Path(path))
            for name, path in output_paths.items()
            if name != "run_info"
        },
    }
    run_info = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).resolve()),
        "argv": sys.argv,
        "analysis_goal": (
            "Correct Supplementary Figure 5D by removing HDR-barcode-positive Kras "
            "G12D/codon-window alleles from the codon 12/13 point-mutation spectrum."
        ),
        "figure_title": FIGURE_TITLE,
        "denominator_type": DENOMINATOR_TYPE,
        "categories": CATEGORY_ORDER,
        "filtering_rule": {
            "include_tumor_samples_only": "cohort != Normal",
            "exclude_hdr_barcode_positive": "strict_hdr_barcode == yes rows are excluded",
            "include_nonhdr_rows": "strict_hdr_barcode == no",
        },
        "classification_rules": {
            "WT": "exact WT codons at both codon 12 and codon 13",
            "mutation_category": (
                "exactly one selected G12/G13 residue change maps to that category"
            ),
            "other": (
                "compound selected mutations, unrequested residues, uncallable codons, "
                "indels/complex calls, or no selected mutation"
            ),
            "selected_g12_residues": sorted(REQUESTED_G12),
            "selected_g13_residues": sorted(REQUESTED_G13),
        },
        "inputs": {"alleles": str(allele_path)},
        "outputs": output_paths,
        "source_accounting": dict(source_accounting),
        "old_bundled_reference": {
            "tumor_denominator_allele_reads": old_aggregate_denominator,
            "aggregate_counts": dict(old_aggregate_counts),
            "aggregate_percentages": {
                category: pct(old_aggregate_counts[category], old_aggregate_denominator)
                for category in CATEGORY_ORDER
            },
        },
        "hdr_positive_removed": {
            "tumor_strict_hdr_barcode_yes_reads": hdr_positive_denominator,
            "aggregate_counts": dict(hdr_positive_counts),
            "aggregate_percent_of_old_denominator": {
                category: pct(hdr_positive_counts[category], old_aggregate_denominator)
                for category in CATEGORY_ORDER
            },
        },
        "corrected_nonhdr": {
            "aggregate_denominator_allele_reads": aggregate_denominator,
            "aggregate_counts": aggregate_counts,
            "aggregate_percentages": aggregate_percentages,
            "per_sample_denominators": {
                specimen: denominator_by_sample[(specimen, cohort)]
                for specimen, cohort in sorted(
                    denominator_by_sample,
                    key=lambda item: sort_key(
                        {"specimen_id": item[0], "cohort": item[1]}
                    ),
                )
            },
            "compound_selected_g12_g13_reads": sum(aggregate_compound.values()),
            "compound_selected_g12_g13_patterns": dict(aggregate_compound),
            "window_class_counts": dict(aggregate_window_classes),
            "top_codon12_codon13_effect_pairs": aggregate_effect_pairs.most_common(20),
            "zoom_ylim": zoom_ylim,
            "interpretation_note": (
                "Corrected G12D is non-HDR/unintended G12D only; HDR-derived Kras "
                "G12D remains represented by the separate HDR panel."
            ),
        },
        "expected_values_checked": {
            "old_bundled_g12d_reads": EXPECTED_OLD_G12D_READS,
            "old_bundled_g12d_pct": pct(
                EXPECTED_OLD_G12D_READS, EXPECTED_OLD_TUMOR_DENOMINATOR
            ),
            "hdr_positive_g12d_removed_reads": EXPECTED_HDR_POSITIVE_G12D_READS,
            "expected_corrected_g12d_reads": EXPECTED_CORRECTED_G12D_READS,
            "expected_corrected_denominator": EXPECTED_CORRECTED_DENOMINATOR,
        },
        "qc_checks": qc_checks,
        "checksums_sha256": checksums,
        "software": {
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "matplotlib_version": mpl.__version__,
            "numpy_version": np.__version__,
        },
    }
    run_info_path.write_text(json.dumps(run_info, indent=2) + "\n")

    print(
        json.dumps(
            {
                "outputs": output_paths,
                "corrected_nonhdr": run_info["corrected_nonhdr"],
                "hdr_positive_removed": run_info["hdr_positive_removed"],
                "qc_checks": qc_checks,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
