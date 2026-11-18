#!/usr/bin/env python3
"""Export plain scientific figures from saved pilot/audit tables (Matplotlib).

Plots are descriptive technical diagnostics for one control prefix. They do
not compare disease groups. No data are generated or filled in by this script.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "preprocessing"))


def make(args):
    from hisat2_pipeline import checksum, save

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if args.output_dir.exists():
        raise ValueError("output folder exists; choose a new folder")
    summary = json.loads(args.summary.read_text())
    with args.sensitivity.open() as handle:
        sensitivity = list(csv.DictReader(handle, delimiter="\t"))
    output = args.output_dir.resolve()
    output.mkdir(parents=True)
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )
    record = {
        "scope": "descriptive diagnostics for one real control prefix",
        "inputs": [
            checksum(args.summary),
            checksum(args.sensitivity),
            checksum(Path(__file__)),
        ],
        "matplotlib_version": matplotlib.__version__,
        "outputs": [],
    }

    def export(fig, name):
        fig.tight_layout()
        for extension in ["svg", "pdf", "png"]:
            path = output / f"{name}.{extension}"
            fig.savefig(
                path, dpi=220, bbox_inches="tight", metadata={"Creator": "Splice Girl"}
            )
            record["outputs"].append(checksum(path))
        plt.close(fig)

    counts = summary["observations_by_class"]
    kinds = [
        "assigned",
        "within_gene_exon_pairs",
        "cross_gene",
        "cross_gene_and_exon_pairs",
        "unannotated",
    ]
    labels = [
        "Unique exon pair",
        "One gene; multiple exon pairs",
        "Multiple genes",
        "Multiple genes and exon pairs",
        "Unannotated",
    ]
    fig, ax = plt.subplots(figsize=(8.8, 4.0))
    values = [counts[kind] for kind in kinds]
    bars = ax.barh(
        labels, values, color=["#333333", "#999999", "#bbbbbb", "#777777", "#dddddd"]
    )
    ax.invert_yaxis()
    ax.set_xlim(0, max(values) * 1.3)
    for bar, value in zip(bars, values):
        ax.text(
            value + max(values) * 0.015,
            bar.get_y() + bar.get_height() / 2,
            f"{value:,}",
            va="center",
            fontsize=9,
        )
    ax.set_xlabel("Read-junction observations after upstream filters")
    ax.set_title(
        "Why junction observations were retained or excluded\nSRR36710387: first 1,000,000 cleaned pairs",
        loc="left",
        fontsize=11,
    )
    export(fig, "junction_assignment")
    strand = summary["spliced_alignment_strand_counts"]
    fractions = []
    for mate in ["read1", "read2"]:
        sense, antisense = (
            strand.get(mate + "_sense", 0),
            strand.get(mate + "_antisense", 0),
        )
        fractions.append((sense / (sense + antisense), antisense / (sense + antisense)))
    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    ax.barh(
        ["Read 1", "Read 2"],
        [row[0] * 100 for row in fractions],
        color="#333333",
        label="Sense to transcript",
    )
    ax.barh(
        ["Read 1", "Read 2"],
        [row[1] * 100 for row in fractions],
        left=[row[0] * 100 for row in fractions],
        color="#bbbbbb",
        label="Antisense to transcript",
    )
    ax.set_xlim(0, 100)
    ax.set_xlabel("Primary unique spliced alignments (%)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False)
    ax.set_title(
        "Mate orientation in the technical pilot\nConsistent with RF; confirm the library protocol",
        loc="left",
        fontsize=11,
    )
    export(fig, "mate_orientation")
    fig, axes = plt.subplots(2, 1, figsize=(6.6, 4.8), sharex=True)
    support = [int(row["minimum_support"]) for row in sensitivity]
    axes[0].plot(
        support,
        [int(row["total_support"]) for row in sensitivity],
        "o-",
        color="#222222",
    )
    axes[0].set_ylabel("Retained junction support")
    axes[0].set_ylim(0, max(int(row["total_support"]) for row in sensitivity) * 1.2)
    axes[0].set_title(
        "CDH11: exploratory support-threshold check\nOne control prefix; no disease comparison",
        loc="left",
        fontsize=11,
    )
    axes[1].plot(
        support,
        [
            float(row["cycle_fraction"]) if row["cycle_fraction"] else float("nan")
            for row in sensitivity
        ],
        "o-",
        color="#222222",
        clip_on=False,
    )
    axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("Cycle fraction")
    axes[1].set_xlabel("Minimum read support per junction")
    axes[1].set_xticks(support)
    axes[1].text(
        0.02, 0.15, "Effectively zero; cycle rank is zero", transform=axes[1].transAxes
    )
    export(fig, "cdh11_support_sensitivity")
    save(output / "figure_manifest.json", record)
    print(f"Saved scientific figures to {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--sensitivity", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    try:
        make(parser.parse_args())
    except (OSError, ValueError, ImportError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
