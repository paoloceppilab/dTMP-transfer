#!/usr/bin/env python3

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Derive a candidate full-length HDR amplicon from donor-positive FASTQ reads "
            "by taking the dominant longest donor-bearing read prefix and filling the "
            "remainder from the WT amplicon."
        )
    )
    parser.add_argument(
        "--fastq",
        nargs="+",
        required=True,
        help="Input FASTQ or FASTQ.GZ files.",
    )
    parser.add_argument(
        "--wt-amplicon",
        required=True,
        help="Full WT amplicon sequence in assay orientation.",
    )
    parser.add_argument(
        "--donor-motif",
        required=True,
        help="Known donor/HDR motif used to identify donor-positive reads.",
    )
    parser.add_argument(
        "--max-reads-per-file",
        type=int,
        default=500000,
        help="Maximum reads to inspect per file. Default: 500000",
    )
    parser.add_argument(
        "--min-donor-reads",
        type=int,
        default=20,
        help="Minimum donor-positive reads to collect before reporting. Default: 20",
    )
    parser.add_argument(
        "--output-json",
        required=True,
        help="Path to write the derivation summary JSON.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    wt = args.wt_amplicon.strip().upper()
    donor = args.donor_motif.strip().upper()
    donor_rc = revcomp(donor)

    counts: Counter[str] = Counter()
    donor_positive = 0
    per_file = []

    for fastq in [Path(p) for p in args.fastq]:
        file_total = 0
        file_donor = 0
        opener = gzip.open if fastq.suffix == ".gz" else open
        with opener(fastq, "rt") as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                seq = handle.readline().strip().upper()
                handle.readline()
                handle.readline()

                file_total += 1
                if file_total > args.max_reads_per_file:
                    break

                if donor in seq:
                    counts[seq] += 1
                    donor_positive += 1
                    file_donor += 1
                elif donor_rc in seq:
                    counts[revcomp(seq)] += 1
                    donor_positive += 1
                    file_donor += 1

        per_file.append(
            {
                "fastq": str(fastq),
                "reads_examined": min(file_total, args.max_reads_per_file),
                "donor_positive_reads": file_donor,
            }
        )

    if donor_positive < args.min_donor_reads:
        raise SystemExit(
            f"Only {donor_positive} donor-positive reads found; "
            f"need at least {args.min_donor_reads}."
        )

    max_len = max(len(seq) for seq in counts)
    top_longest = [(seq, count) for seq, count in counts.items() if len(seq) == max_len]
    top_longest.sort(key=lambda item: (-item[1], item[0]))
    longest_prefix, longest_prefix_count = top_longest[0]

    expected_hdr = longest_prefix + wt[len(longest_prefix) :]
    variants = [
        {"position_1based": i, "wt": w, "hdr": h}
        for i, (w, h) in enumerate(zip(wt, expected_hdr), start=1)
        if w != h
    ]

    payload = {
        "wt_amplicon": wt,
        "donor_motif": donor,
        "donor_positive_reads_total": donor_positive,
        "dominant_longest_prefix": longest_prefix,
        "dominant_longest_prefix_length": len(longest_prefix),
        "dominant_longest_prefix_count": longest_prefix_count,
        "expected_hdr_amplicon_seq": expected_hdr,
        "variant_sites_vs_wt": variants,
        "per_file_summary": per_file,
        "top_sequences": [
            {"sequence": seq, "length": len(seq), "count": count}
            for seq, count in counts.most_common(10)
        ],
        "assumption": (
            "Bases beyond the dominant donor-positive read prefix are copied from the WT "
            "amplicon because the available reads do not span the full 134 bp amplicon."
        ),
    }

    output_path = Path(args.output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
