#!/usr/bin/env python3
"""Reproduce the Mach-5 order-eight collision-production comparison."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.lines import Line2D


STYLES = {
    "j2l4": {"label": r"$(j,\ell)=(2,4)$", "marker": "o", "color": "#1f77b4"},
    "j3l2": {"label": r"$(j,\ell)=(3,2)$", "marker": "s", "color": "#ff7f0e"},
}

# Deliberate offsets keep each station sign clear of its symbol.
LABEL_OFFSETS = {
    ("j2l4", -6): (12, 8),
    ("j2l4", -3): (10, 9),
    ("j2l4", 0): (10, -12),
    ("j2l4", 3): (10, 8),
    ("j2l4", 6): (-13, -15),
    ("j3l2", -6): (13, 9),
    ("j3l2", -3): (10, -13),
    ("j3l2", 0): (10, 8),
    ("j3l2", 3): (10, 8),
    ("j3l2", 6): (11, 12),
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

        for px, py in zip(x, y):
            ax.plot([px, px], [px, py], color=style["color"], linewidth=1.15,
                    alpha=0.68, zorder=2)

        # Open symbols and horizontal bars: theorem estimates on y=x.
        ax.errorbar(
            x, x, xerr=xerr, linestyle="none", capsize=4,
            markersize=10.5, marker=style["marker"], color=style["color"],
            markerfacecolor="white", markeredgecolor=style["color"],
            markeredgewidth=1.8, zorder=3,
        )

        # Filled symbols and vertical bars: direct weak-form estimates.
        ax.errorbar(
            x, y, yerr=yerr, linestyle="none", capsize=4,
            markersize=6.8, marker=style["marker"], color=style["color"],
            markerfacecolor=style["color"], markeredgecolor=style["color"],
            zorder=4,
        )

        for row, px, py in zip(subset, x, y):
            station = int(row["offset"])
            dx, dy = LABEL_OFFSETS[(obs, station)]
            annotation = ax.annotate(
                f"{station:+d}", (px, py), xytext=(dx, dy),
                textcoords="offset points", fontsize=11, zorder=6,
            )
            annotation.set_path_effects([
                path_effects.Stroke(linewidth=3.2, foreground="white"),
                path_effects.Normal(),
            ])

    limit = 5.0e4
    ax.plot([-limit, limit], [-limit, limit], "--", color="#2ca02c",
            linewidth=1.6, zorder=1)
    ax.set_xscale("symlog", linthresh=10)
    ax.set_yscale("symlog", linthresh=10)
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_xlabel("Single-generator theorem prediction", fontsize=15)
    ax.set_ylabel("Direct hard-sphere weak-form sampling", fontsize=15)
    ax.set_title("Mach-5 shock: blind order-8 collision-production test", fontsize=17)
    ax.grid(True, which="major", alpha=0.25)
    ax.tick_params(labelsize=12)
    ax.text(
        0.03, 0.975,
        "fit to filled symbols, constrained through origin\n"
        f"slope = {slope:.6f} (not shown)\n"
        f"normalized RMS mismatch = {nrms:.3f}%",
        transform=ax.transAxes, va="top", fontsize=13.5,
    )

    role_handles = [
        Line2D([], [], linestyle="none", marker="D", markersize=8.5,
               markerfacecolor="white", markeredgecolor="black", markeredgewidth=1.5,
               label="theorem value on $y=x$"),
        Line2D([], [], linestyle="none", marker="D", markersize=6.2,
               markerfacecolor="black", markeredgecolor="black",
               label="direct weak-form sample"),
    ]
    role_legend = ax.legend(handles=role_handles, loc="upper right", fontsize=11.5,
                            framealpha=0.93, handletextpad=0.7)
    ax.add_artist(role_legend)

    observable_handles = [
        Line2D([], [], linestyle="--", color="#2ca02c", linewidth=1.6,
               label=r"perfect match: $y=x$ (slope 1)"),
        Line2D([], [], linestyle="none", marker="o", markersize=8,
               markerfacecolor=STYLES["j2l4"]["color"],
               markeredgecolor=STYLES["j2l4"]["color"], label=STYLES["j2l4"]["label"]),
        Line2D([], [], linestyle="none", marker="s", markersize=8,
               markerfacecolor=STYLES["j3l2"]["color"],
               markeredgecolor=STYLES["j3l2"]["color"], label=STYLES["j3l2"]["label"]),
    ]
    ax.legend(handles=observable_handles, loc="lower right", fontsize=12.5,
              framealpha=0.93)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160)
    plt.close(fig)
    print(f"slope_through_origin={slope:.9f}")
    print(f"normalized_rms_mismatch_percent={nrms:.6f}")


if __name__ == "__main__":
    main()
