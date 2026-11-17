#!/usr/bin/env python3
"""Recompute descriptive Hodge effects across an explicitly supplied support grid.

This checks count/coverage sensitivity only. It does not change alignment,
rerun anchor filtering, select a preferred threshold, or perform inference.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "preprocessing"))
from hisat2_pipeline import checksum, save


def thresholds(text):
    values = sorted(set(float(part) for part in text.split(",")))
    if not values or any(not math.isfinite(value) or value < 0 for value in values):
        raise ValueError(
            "thresholds must be comma-separated finite nonnegative numbers"
        )
    return values


def retained_rows(rows, minimum_edge, minimum_gene):
    retained = [
        row
        for row in rows
        if float(row["count"]) >= minimum_edge and float(row["count"]) > 0
    ]
    totals = defaultdict(float)
    for row in retained:
        totals[(row["sample_id"], row["gene_id"])] += float(row["count"])
    return [
        row
        for row in retained
        if totals[(row["sample_id"], row["gene_id"])] >= minimum_gene
    ]


def run(args):
    edges, coverage = (
        thresholds(args.minimum_junction_counts),
        thresholds(args.minimum_gene_support),
    )
    output = args.output_dir.resolve()
    if output.exists():
        raise ValueError(f"output folder already exists: {output}")
    binary = args.binary.resolve()
    # Check the complete original input before filtering can hide bad rows.
    subprocess.run([str(binary), "validate", str(args.junctions)], check=True)
    with args.junctions.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames
        rows = list(reader)
    genes = sorted({row["gene_id"] for row in rows})
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()
    manifest = {
        "status": "running",
        "mode": "descriptive support sensitivity",
        "inputs": [checksum(args.junctions), checksum(args.manifest), checksum(binary)],
        "scenarios": [],
        "outputs": [],
    }
    save(output / "run_manifest.json", manifest)
    try:
        command = [
            str(binary),
            "compare",
            str(args.junctions.resolve()),
            str(args.manifest.resolve()),
            "--output",
            str(output / "unfiltered"),
        ]
        with (output / "unfiltered.log").open("w") as log:
            result = subprocess.run(
                command, stdout=log, stderr=subprocess.STDOUT, check=False
            )
        if (
            result.returncode not in (0, 1)
            or not (output / "unfiltered/gene_summaries.tsv").is_file()
        ):
            raise ValueError(
                "unfiltered input/manifest analysis failed; see unfiltered.log"
            )
        manifest["unfiltered_command"] = command
        manifest["unfiltered_manifest"] = checksum(
            output / "unfiltered/run_manifest.json"
        )
        with (output / "sensitivity.tsv").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "scenario",
                    "minimum_junction_count",
                    "minimum_gene_support",
                    "gene_id",
                    "status",
                    "case_n",
                    "control_n",
                    "mean_difference_case_minus_control",
                    "retained_rows",
                ]
            )
            for index, (edge, total) in enumerate(
                itertools.product(edges, coverage), 1
            ):
                name = f"scenario_{index:03d}"
                directory = output / name
                directory.mkdir()
                kept = retained_rows(rows, edge, total)
                table = directory / "junctions.tsv"
                with table.open("x", encoding="utf-8", newline="") as target:
                    table_writer = csv.DictWriter(
                        target, fieldnames=header, delimiter="\t", lineterminator="\n"
                    )
                    table_writer.writeheader()
                    table_writer.writerows(kept)
                record = {
                    "scenario": name,
                    "minimum_junction_count": edge,
                    "minimum_gene_support": total,
                    "retained_rows": len(kept),
                    "input": checksum(table),
                }
                manifest["scenarios"].append(record)
                summaries = {}
                if kept:
                    command = [
                        str(binary),
                        "compare",
                        str(table),
                        str(args.manifest.resolve()),
                        "--output",
                        str(directory / "results"),
                    ]
                    record["command"] = command
                    with (directory / "analysis.log").open("w") as log:
                        result = subprocess.run(
                            command, stdout=log, stderr=subprocess.STDOUT, check=False
                        )
                    # compare saves a diagnostic bundle when no group pair survives.
                    summary_path = directory / "results/gene_summaries.tsv"
                    if not summary_path.is_file() or result.returncode not in (0, 1):
                        raise ValueError(f"{name}: analysis failed; see analysis.log")
                    with summary_path.open() as source:
                        summaries = {
                            row["gene_id"]: row
                            for row in csv.DictReader(source, delimiter="\t")
                        }
                    if result.returncode and any(
                        row["status"] == "complete" for row in summaries.values()
                    ):
                        raise ValueError(
                            f"{name}: unexpected analysis failure; see analysis.log"
                        )
                    record["analysis_manifest"] = checksum(
                        directory / "results/run_manifest.json"
                    )
                for gene in genes:
                    summary = summaries.get(gene, {})
                    writer.writerow(
                        [
                            name,
                            edge,
                            total,
                            gene,
                            summary.get("status", "no_retained_junctions"),
                            summary.get("case_n", "0"),
                            summary.get("control_n", "0"),
                            summary.get("mean_difference_case_minus_control", ""),
                            len(kept),
                        ]
                    )
                record["status"] = "complete" if kept else "no_retained_junctions"
                save(output / "run_manifest.json", manifest)
        manifest["outputs"] = [checksum(output / "sensitivity.tsv")]
        manifest["status"] = "complete"
    except Exception as error:
        manifest.update(status="failed", error=str(error))
        raise
    finally:
        save(output / "run_manifest.json", manifest)
    print(f"Descriptive sensitivity table: {output / 'sensitivity.tsv'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("junctions", "manifest", "output-dir"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--minimum-junction-counts", required=True)
    parser.add_argument("--minimum-gene-support", required=True)
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target/release/splice-girl"
    )
    try:
        run(parser.parse_args())
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
