#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


COHORT_ORDER = {"Normal": 0, "KP": 1, "KPTT": 2}
GENE_COLORS = {
    "Kras": "#7FBF7B",
    "Trp53": "#F67E5F",
}
MISSING_COLOR = "#D9D9D9"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Figure 7B-like charts from editing summary data."
    )
    parser.add_argument("--editing-summary", required=True)
    parser.add_argument("--pdf-out", required=True)
    parser.add_argument("--png-out", required=True)
    parser.add_argument("--svg-out", required=True)
    parser.add_argument(
        "--layout",
        choices=["shared_axis_inset", "shared_axis_no_inset", "two_panel_no_zoom"],
        default="two_panel_no_zoom",
        help=(
            "Figure layout to render. Default matches the primary manuscript-style "
            "two-panel output."
        ),
    )
    return parser.parse_args()


def read_summary(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sort_key(row: dict[str, str]) -> tuple[int, int, str]:
    suffix = row["specimen_id"].split("_")[-1]
    try:
        numeric_suffix = int(suffix)
    except ValueError:
        numeric_suffix = 10**9
    return (COHORT_ORDER.get(row["cohort"], 99), numeric_suffix, row["specimen_id"])


def parse_optional_float(value: str) -> float | None:
    if value in {"", "NA", None}:
        return None
    return float(value)


def save_figure(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")


def add_cohort_brackets(
    axis: plt.Axes,
    cohorts: list[str],
    y_text: float = -0.47,
    y_line: float = -0.32,
) -> None:
    cohort_ranges: list[tuple[str, int, int]] = []
    start = 0
    for idx in range(1, len(cohorts) + 1):
        if idx == len(cohorts) or cohorts[idx] != cohorts[start]:
            cohort_ranges.append((cohorts[start], start, idx - 1))
            start = idx

    for _, _, end in cohort_ranges[:-1]:
        axis.axvline(end + 0.5, color="#999999", linewidth=0.8)

    for cohort, start, end in cohort_ranges:
        midpoint = (start + end) / 2
        axis.plot(
            [start - 0.35, end + 0.35],
            [y_line, y_line],
            transform=axis.get_xaxis_transform(),
            color="#777777",
            linewidth=1.0,
            clip_on=False,
        )
        axis.plot(
            [start - 0.35, start - 0.35],
            [y_line - 0.01, y_line + 0.01],
            transform=axis.get_xaxis_transform(),
            color="#777777",
            linewidth=1.0,
            clip_on=False,
        )
        axis.plot(
            [end + 0.35, end + 0.35],
            [y_line - 0.01, y_line + 0.01],
            transform=axis.get_xaxis_transform(),
            color="#777777",
            linewidth=1.0,
            clip_on=False,
        )
        axis.text(
            midpoint,
            y_text,
            cohort,
            ha="center",
            va="top",
            transform=axis.get_xaxis_transform(),
            fontsize=10,
        )


def plot_shared_axis(
    specimen_labels: list[str],
    cohorts: list[str],
    kras_values: list[float | None],
    trp53_values: list[float | None],
    show_inset: bool,
) -> plt.Figure:
    x_positions = list(range(len(specimen_labels)))
    fig, ax = plt.subplots(figsize=(15.5, 7.2))

    bar_width = 0.36
    missing_height = 0.8
    kras_plot_values = [value if value is not None else 0.0 for value in kras_values]
    trp53_plot_values = [value if value is not None else 0.0 for value in trp53_values]

    ax.bar(
        [x - bar_width / 2 for x in x_positions],
        kras_plot_values,
        width=bar_width,
        color=GENE_COLORS["Kras"],
        edgecolor="black",
        linewidth=0.4,
        label="Kras HDR",
    )
    ax.bar(
        [x + bar_width / 2 for x in x_positions],
        trp53_plot_values,
        width=bar_width,
        color=GENE_COLORS["Trp53"],
        edgecolor="black",
        linewidth=0.4,
        label="Trp53 Indel",
    )

    for xpos, value in zip(x_positions, kras_values):
        if value is None:
            ax.bar(
                xpos - bar_width / 2,
                missing_height,
                width=bar_width,
                color=MISSING_COLOR,
                edgecolor="black",
                linewidth=0.4,
                hatch="//",
            )
            ax.text(
                xpos - bar_width / 2,
                missing_height + 0.5,
                "NA",
                ha="center",
                va="bottom",
                fontsize=8,
                rotation=90,
            )

    for xpos, value in zip(x_positions, trp53_values):
        if value is None:
            ax.bar(
                xpos + bar_width / 2,
                missing_height,
                width=bar_width,
                color=MISSING_COLOR,
                edgecolor="black",
                linewidth=0.4,
                hatch="//",
            )
            ax.text(
                xpos + bar_width / 2,
                missing_height + 0.5,
                "NA",
                ha="center",
                va="bottom",
                fontsize=8,
                rotation=90,
            )

    ax.set_ylabel("Mutation % (Kras HDR and Trp53 indel)")
    ax.set_xlabel("Specimen", labelpad=18)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.6)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(specimen_labels, rotation=58, ha="right")
    ax.tick_params(axis="x", pad=4)
    ax.set_xlim(-0.8, len(specimen_labels) - 0.2)

    add_cohort_brackets(ax, cohorts)

    legend_handles = [
        Patch(facecolor=GENE_COLORS["Kras"], edgecolor="black", label="Kras HDR"),
        Patch(facecolor=GENE_COLORS["Trp53"], edgecolor="black", label="Trp53 Indel"),
        Patch(facecolor=MISSING_COLOR, edgecolor="black", hatch="//", label="Missing"),
    ]
    ax.legend(handles=legend_handles, frameon=False, ncol=3, loc="upper right")

    if show_inset:
        inset = ax.inset_axes([0.08, 0.57, 0.24, 0.28])
        inset.set_facecolor("white")
        inset.bar(
            x_positions,
            kras_plot_values,
            width=0.55,
            color=GENE_COLORS["Kras"],
            edgecolor="black",
            linewidth=0.35,
        )
        inset.set_xlim(-0.6, len(specimen_labels) - 0.4)
        valid_kras = [value for value in kras_values if value is not None]
        kras_ymax = max(valid_kras) * 1.2 if valid_kras else 0.01
        inset.set_ylim(0, max(kras_ymax, 0.01))
        inset.set_xticks([])
        inset.tick_params(axis="y", labelsize=8)
        inset.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.5)
        inset.spines["top"].set_visible(False)
        inset.spines["right"].set_visible(False)
        inset.set_ylabel("%", fontsize=8)
        inset.text(
            0.02,
            0.96,
            "Kras HDR (zoom)",
            transform=inset.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5},
        )

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.34)
    return fig


def plot_two_panel_no_zoom(
    specimen_labels: list[str],
    cohorts: list[str],
    kras_values: list[float | None],
    trp53_values: list[float | None],
) -> plt.Figure:
    x_positions = list(range(len(specimen_labels)))
    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        figsize=(15.5, 8.2),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 2]},
    )

    kras_plot_values = [value if value is not None else 0.0 for value in kras_values]
    trp53_plot_values = [value if value is not None else 0.0 for value in trp53_values]

    ax_top.bar(
        x_positions,
        kras_plot_values,
        color=GENE_COLORS["Kras"],
        edgecolor="black",
        linewidth=0.4,
        width=0.72,
    )
    ax_top.set_ylabel("Kras HDR (%)")
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.6)
    valid_kras = [value for value in kras_values if value is not None]
    kras_ymax = max(valid_kras) * 1.2 if valid_kras else 0.01
    ax_top.set_ylim(0, max(kras_ymax, 0.01))

    for xpos, value in zip(x_positions, kras_values):
        if value is None:
            ax_top.bar(
                xpos,
                max(kras_ymax * 0.2, 0.001),
                color=MISSING_COLOR,
                edgecolor="black",
                linewidth=0.4,
                hatch="//",
                width=0.72,
            )
            ax_top.text(
                xpos,
                max(kras_ymax * 0.22, 0.0012),
                "NA",
                ha="center",
                va="bottom",
                fontsize=8,
                rotation=90,
            )

    ax_bottom.bar(
        x_positions,
        trp53_plot_values,
        color=GENE_COLORS["Trp53"],
        edgecolor="black",
        linewidth=0.4,
        width=0.72,
    )
    ax_bottom.set_ylabel("Trp53 Indel (%)")
    ax_bottom.set_xlabel("Specimen", labelpad=18)
    ax_bottom.spines["top"].set_visible(False)
    ax_bottom.spines["right"].set_visible(False)
    ax_bottom.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.6)

    for xpos, value in zip(x_positions, trp53_values):
        if value is None:
            ax_bottom.bar(
                xpos,
                0.8,
                color=MISSING_COLOR,
                edgecolor="black",
                linewidth=0.4,
                hatch="//",
                width=0.72,
            )
            ax_bottom.text(
                xpos,
                1.0,
                "NA",
                ha="center",
                va="bottom",
                fontsize=8,
                rotation=90,
            )

    ax_bottom.set_xticks(x_positions)
    ax_bottom.set_xticklabels(specimen_labels, rotation=58, ha="right")
    ax_bottom.tick_params(axis="x", pad=4)
    ax_top.set_xlim(-0.8, len(specimen_labels) - 0.2)
    ax_bottom.set_xlim(-0.8, len(specimen_labels) - 0.2)

    add_cohort_brackets(ax_bottom, cohorts)

    legend_handles = [
        Patch(facecolor=GENE_COLORS["Kras"], edgecolor="black", label="Kras HDR"),
        Patch(facecolor=GENE_COLORS["Trp53"], edgecolor="black", label="Trp53 Indel"),
        Patch(facecolor=MISSING_COLOR, edgecolor="black", hatch="//", label="Missing"),
    ]
    ax_top.legend(handles=legend_handles, frameon=False, ncol=3, loc="upper right")

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.34)
    return fig


def main() -> None:
    args = parse_args()
    rows = sorted(read_summary(Path(args.editing_summary)), key=sort_key)

    specimen_labels = [row["specimen_id"] for row in rows]
    cohorts = [row["cohort"] for row in rows]
    kras_values = [parse_optional_float(row["kras_hdr_exact_pct"]) for row in rows]
    trp53_values = [parse_optional_float(row["trp53_indel_pct"]) for row in rows]

    if args.layout == "shared_axis_inset":
        fig = plot_shared_axis(
            specimen_labels,
            cohorts,
            kras_values,
            trp53_values,
            show_inset=True,
        )
    elif args.layout == "shared_axis_no_inset":
        fig = plot_shared_axis(
            specimen_labels,
            cohorts,
            kras_values,
            trp53_values,
            show_inset=False,
        )
    else:
        fig = plot_two_panel_no_zoom(specimen_labels, cohorts, kras_values, trp53_values)

    save_figure(fig, Path(args.pdf_out))
    save_figure(fig, Path(args.png_out))
    save_figure(fig, Path(args.svg_out))


if __name__ == "__main__":
    main()
