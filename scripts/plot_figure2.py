#!/usr/bin/env python3
"""Reproduce the Mach-5 order-eight collision-production comparison."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


STYLES = {
    "j2l4": {"label": r"$(j,\ell)=(2,4)$", "marker": "o", "color": "#1f77b4"},
    "j3l2": {"label": r"$(j,\ell)=(3,2)$", "marker": "s", "color": "#ff7f0e"},
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument("-o", "--output", type=Path, default=Path("figure2_mach5.png"))
    args = parser.parse_args()

    with args.data.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    theorem = [float(row["tree"]) for row in rows]
    direct = [float(row["direct"]) for row in rows]
    slope = sum(x * y for x, y in zip(theorem, direct)) / sum(x * x for x in theorem)
    nrms = 100.0 * math.sqrt(
        sum((x - y) ** 2 for x, y in zip(theorem, direct)) / sum(x * x for x in theorem)
    )

    fig, ax = plt.subplots(figsize=(9.9, 8.55), constrained_layout=True)
    for obs, style in STYLES.items():
        subset = [row for row in rows if row["obs"] == obs]
        x = [float(row["tree"]) for row in subset]
        y = [float(row["direct"]) for row in subset]
        xerr = [float(row["tree_se"]) for row in subset]
        yerr = [float(row["direct_se"]) for row in subset]
        ax.errorbar(
            x, y, xerr=xerr, yerr=yerr, linestyle="none", capsize=4,
            markersize=8, marker=style["marker"], color=style["color"],
            label=style["label"],
        )
        for row, px, py in zip(subset, x, y):
            offset = int(row["offset"])
            ax.annotate(f"{offset:+d}", (px, py), xytext=(6, 4),
                        textcoords="offset points", fontsize=11)

    limit = 1.4e5
    ax.plot([-limit, limit], [-limit, limit], "--", color="#2ca02c",
            linewidth=1.6, label="perfect match")
    ax.set_xscale("symlog", linthresh=10)
    ax.set_yscale("symlog", linthresh=10)
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_xlabel("Single-generator theorem prediction", fontsize=15)
    ax.set_ylabel("Direct hard-sphere weak-form sampling", fontsize=15)
    ax.set_title("Mach-5 shock: blind order-8 collision-production test", fontsize=17)
    ax.grid(True, which="major", alpha=0.25)
    ax.text(
        0.03, 0.975,
        f"slope through origin = {slope:.6f}\nnormalized RMS mismatch = {nrms:.3f}%",
        transform=ax.transAxes, va="top", fontsize=14,
    )

    handles, labels = ax.get_legend_handles_labels()
    order = [
        labels.index("perfect match"),
        labels.index(STYLES["j2l4"]["label"]),
        labels.index(STYLES["j3l2"]["label"]),
    ]
    ax.legend([handles[i] for i in order], [labels[i] for i in order],
              loc="lower right", fontsize=13)
    ax.tick_params(labelsize=12)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160)
    plt.close(fig)
    print(f"slope_through_origin={slope:.9f}")
    print(f"normalized_rms_mismatch_percent={nrms:.6f}")


if __name__ == "__main__":
    main()
