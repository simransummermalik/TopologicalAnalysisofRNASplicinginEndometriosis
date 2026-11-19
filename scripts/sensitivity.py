#!/usr/bin/env python3
"""Recompute Hodge effects across explicit sample or shared-graph grids.

The default sample mode is descriptive. Shared mode can optionally run Rust
statistics with explicit inference settings. Neither mode changes alignment,
reruns anchor filtering or selects a preferred threshold.
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
from shared_graph import (
    build_shared,
    export_shared,
    load_cohort,
    positive_integer,
    read_table,
)


def thresholds(text):
    if not isinstance(text, str) or not text.strip():
        raise ValueError(
            "thresholds must be comma-separated finite nonnegative numbers"
        )
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
    if getattr(args, "mode", "sample") == "shared":
        return run_shared(args)
    if getattr(args, "run_stats", False):
        raise ValueError("--run-stats is supported only with --mode shared")
    edges, coverage = (
        thresholds(args.minimum_junction_counts),
        thresholds(args.minimum_gene_support),
    )
    output = args.output_dir.resolve()
    if output.exists():
        raise ValueError(f"output folder already exists: {output}")
    binary = args.binary.resolve()
    # Check the complete original input before filtering can hide bad rows.
    validation = subprocess.run(
        [str(binary), "validate", str(args.junctions)],
        capture_output=True,
        text=True,
        check=False,
    )
    if validation.returncode:
        raise ValueError(f"input validation failed: {validation.stderr.strip()}")
    with args.junctions.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames
        rows = list(reader)
    genes = sorted({row["gene_id"] for row in rows})
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()
    (output / "validation.log").write_text(validation.stdout + validation.stderr)
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


def shared_scenarios(args):
    """Primary plus unique one-factor changes; never a Cartesian product."""
    primary = (
        args.primary_edge_count,
        args.primary_sample_count,
        args.primary_gene_support,
    )
    if any(value is None for value in primary):
        raise ValueError(
            "shared mode requires --primary-edge-count, --primary-sample-count "
            "and --primary-gene-support"
        )
    grids = [
        thresholds(args.minimum_junction_counts),
        thresholds(args.minimum_samples_per_edge),
        thresholds(args.minimum_gene_support),
    ]
    if any(value < 1 or value != int(value) for grid in grids for value in grid):
        raise ValueError("shared-mode grids must contain positive integers")
    scenarios = [("primary", primary)]
    seen = {primary}
    for column, grid in enumerate(grids):
        for value in grid:
            setting = list(primary)
            setting[column] = int(value)
            setting = tuple(setting)
            if setting not in seen:
                seen.add(setting)
                scenarios.append((f"scenario_{len(scenarios):03d}", setting))
    return scenarios


def run_shared(args):
    scenarios = shared_scenarios(args)
    cohort = load_cohort(
        args.junctions, args.manifest, args.genes, args.complete_counts
    )
    if any(setting[1] > len(cohort.samples) for _, setting in scenarios):
        raise ValueError("minimum samples per edge exceeds the manifest cohort size")
    if args.run_stats:
        if not all(
            (args.genes, args.family, args.alternative, args.design, args.min_per_group)
        ):
            raise ValueError(
                "--run-stats requires --genes, --family, --alternative, "
                "--design independent and --min-per-group"
            )
        if args.min_per_group < 2:
            raise ValueError("--min-per-group must be at least 2")
        _, samples = read_table(args.manifest, ["sample_id", "condition"])
        if any(
            sum(row["condition"].lower() == group for row in samples)
            < args.min_per_group
            for group in ("case", "control")
        ):
            raise ValueError("manifest does not meet --min-per-group")
    output, binary = args.output_dir.resolve(), args.binary.resolve()
    if output.exists():
        raise ValueError(f"output folder already exists: {output}")
    inputs = [
        checksum(path)
        for path in (args.junctions, args.manifest, args.genes, binary)
        if path
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()
    manifest = dict(
        status="running",
        mode="shared graph one-factor sensitivity",
        complete_counts_declared=True,
        approval="explicit caller settings; not team approval",
        inputs=inputs,
        cohort=cohort.samples,
        scenarios=[],
        outputs=[],
        implementation=[
            checksum(Path(__file__)),
            checksum(ROOT / "scripts/shared_graph.py"),
        ],
        anchor_sensitivity="not evaluated; requires separate read-counter exports",
        statistics=(
            dict(
                family=args.family,
                alternative=args.alternative,
                design=args.design,
                min_per_group=args.min_per_group,
                role="sensitivity; no automatic primary selection or multiplicity pooling",
            )
            if args.run_stats
            else None
        ),
    )
    save(output / "run_manifest.json", manifest)
    try:
        with (output / "sensitivity.tsv").open(
            "x", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "scenario",
                    "minimum_junction_count",
                    "minimum_samples_per_edge",
                    "minimum_gene_support",
                    "gene_id",
                    "status",
                    "shared_edges",
                    "cycle_rank",
                    "covered_samples",
                    "required_samples",
                    "case_n",
                    "control_n",
                    "mean_difference_case_minus_control",
                ]
            )
            for name, (edge, samples, coverage) in scenarios:
                result = build_shared(cohort, edge, samples, coverage)
                directory = output / name
                export_shared(directory, cohort, result, args.manifest)
                record = dict(
                    scenario=name,
                    status="running",
                    minimum_junction_count=edge,
                    minimum_samples_per_edge=samples,
                    minimum_gene_support=coverage,
                    eligible_genes=len(result["eligible"]),
                    retained_rows=len(result["junctions"]),
                    inputs=[
                        checksum(directory / file)
                        for file in (
                            "junctions.tsv",
                            "samples.tsv",
                            "eligible_genes.tsv",
                        )
                    ],
                )
                manifest["scenarios"].append(record)
                save(output / "run_manifest.json", manifest)
                summaries = {}
                if result["junctions"]:
                    command = [
                        str(binary),
                        "compare",
                        str(directory / "junctions.tsv"),
                        str(directory / "samples.tsv"),
                        "--output",
                        str(directory / "results"),
                    ]
                    record["command"] = command
                    with (directory / "analysis.log").open(
                        "x", encoding="utf-8"
                    ) as log:
                        process = subprocess.run(
                            command, stdout=log, stderr=subprocess.STDOUT, check=False
                        )
                    record["returncode"] = process.returncode
                    summary_path = directory / "results/gene_summaries.tsv"
                    if process.returncode not in (0, 1) or not summary_path.is_file():
                        raise ValueError(f"{name}: analysis failed; see analysis.log")
                    with summary_path.open(encoding="utf-8") as source:
                        summaries = {
                            row["gene_id"]: row
                            for row in csv.DictReader(source, delimiter="\t")
                        }
                    if process.returncode and any(
                        row["status"] == "complete" for row in summaries.values()
                    ):
                        raise ValueError(
                            f"{name}: unexpected analysis failure; see analysis.log"
                        )
                    record["analysis_manifest"] = checksum(
                        directory / "results/run_manifest.json"
                    )
                for gene in result["genes"]:
                    summary = summaries.get(gene["gene_id"], {})
                    if (
                        gene["status"] == "eligible"
                        and summary.get("status") != "complete"
                    ):
                        raise ValueError(
                            f"{name}: eligible gene missing a complete Rust comparison"
                        )
                    writer.writerow(
                        [
                            name,
                            edge,
                            samples,
                            coverage,
                            gene["gene_id"],
                            gene["status"],
                            gene["shared_edges"],
                            gene["cycle_rank"],
                            gene["covered_samples"],
                            gene["required_samples"],
                            summary.get("case_n", "0"),
                            summary.get("control_n", "0"),
                            (
                                summary.get("mean_difference_case_minus_control", "")
                                if gene["status"] == "eligible"
                                else ""
                            ),
                        ]
                    )
                if args.run_stats and result["eligible"]:
                    command = [
                        str(binary),
                        "stats",
                        str(directory / "junctions.tsv"),
                        str(directory / "samples.tsv"),
                        "--genes",
                        str(directory / "eligible_genes.tsv"),
                        "--family",
                        args.family,
                        "--alternative",
                        args.alternative,
                        "--design",
                        args.design,
                        "--min-per-group",
                        str(args.min_per_group),
                        "--output",
                        str(directory / "statistics"),
                    ]
                    record["statistics_command"] = command
                    with (directory / "statistics.log").open(
                        "x", encoding="utf-8"
                    ) as log:
                        process = subprocess.run(
                            command, stdout=log, stderr=subprocess.STDOUT, check=False
                        )
                    record["statistics_returncode"] = process.returncode
                    if process.returncode:
                        raise ValueError(
                            f"{name}: statistics failed; see statistics.log"
                        )
                    record["statistics_manifest"] = checksum(
                        directory / "statistics/run_manifest.json"
                    )
                    record["statistics_status"] = "complete"
                else:
                    record["statistics_status"] = (
                        "no_eligible_genes" if args.run_stats else "not_requested"
                    )
                record["status"] = "complete"
                handle.flush()
                record["outputs"] = [
                    checksum(path)
                    for path in sorted(directory.iterdir())
                    if path.is_file()
                ]
                save(output / "run_manifest.json", manifest)
        manifest.update(
            status="complete", outputs=[checksum(output / "sensitivity.tsv")]
        )
    except BaseException as error:
        manifest.update(
            status="interrupted" if isinstance(error, KeyboardInterrupt) else "failed",
            error=str(error),
        )
        if manifest["scenarios"] and manifest["scenarios"][-1]["status"] == "running":
            manifest["scenarios"][-1].update(
                status=manifest["status"], error=str(error)
            )
        raise
    finally:
        save(output / "run_manifest.json", manifest)
    print(f"Shared-graph sensitivity table: {output / 'sensitivity.tsv'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("junctions", "manifest", "output-dir"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--minimum-junction-counts", required=True)
    parser.add_argument("--minimum-gene-support", required=True)
    parser.add_argument("--mode", choices=["sample", "shared"], default="sample")
    parser.add_argument("--genes", type=Path)
    parser.add_argument("--complete-counts", action="store_true")
    parser.add_argument("--minimum-samples-per-edge")
    parser.add_argument("--primary-edge-count", type=positive_integer)
    parser.add_argument("--primary-sample-count", type=positive_integer)
    parser.add_argument("--primary-gene-support", type=positive_integer)
    parser.add_argument("--run-stats", action="store_true")
    parser.add_argument("--family")
    parser.add_argument(
        "--alternative", choices=["two-sided", "case-greater", "case-less"]
    )
    parser.add_argument("--design", choices=["independent"])
    parser.add_argument("--min-per-group", type=positive_integer)
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
