#!/usr/bin/env python3
"""Aggregate block CSVs written by order8_m5_validation.cpp."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import fmean, stdev


def standard_error(values: list[float]) -> float:
    return stdev(values) / math.sqrt(len(values)) if len(values) > 1 else 0.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path, help="*_blocks.csv files")
    parser.add_argument("-o", "--output", required=True, type=Path)
    args = parser.parse_args()

    groups: dict[tuple[int, str], list[dict[str, str]]] = defaultdict(list)
    for path in args.inputs:
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                groups[(int(row["offset"]), row["obs"])].append(row)

    fieldnames = [
        "offset", "obs", "relx", "tree", "direct", "tree_se", "direct_se",
        "difference", "difference_se", "n_blocks", "paired_z",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for (offset, obs), rows in sorted(groups.items()):
            relx = [float(row["relx"]) for row in rows]
            theorem = [float(row["tree"]) for row in rows]
            direct = [float(row["direct"]) for row in rows]
            difference = [d - t for t, d in zip(theorem, direct)]
            diff_mean = fmean(difference)
            diff_se = standard_error(difference)
            writer.writerow(
                {
                    "offset": offset,
                    "obs": obs,
                    "relx": fmean(relx),
                    "tree": fmean(theorem),
                    "direct": fmean(direct),
                    "tree_se": standard_error(theorem),
                    "direct_se": standard_error(direct),
                    "difference": diff_mean,
                    "difference_se": diff_se,
                    "n_blocks": len(rows),
                    "paired_z": diff_mean / diff_se if diff_se else float("nan"),
                }
            )


if __name__ == "__main__":
    main()
