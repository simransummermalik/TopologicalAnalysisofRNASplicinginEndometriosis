#!/usr/bin/env python3
"""Measure Rust analysis/export on observed candidate genes in a real pilot.

This is a technical benchmark, not a case/control analysis. Missing candidates
are counted and never replaced with zero-support graphs. Each timed process
has its own resource monitor so peak memory is not inherited from earlier runs.
"""

import argparse
import csv
import json
import platform
import statistics
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "preprocessing"))

WORKER = """
import json, resource, subprocess, sys, time
started = time.perf_counter()
result = subprocess.run(sys.argv[1:], stdout=subprocess.DEVNULL,
                        stderr=subprocess.PIPE, text=True, timeout=120)
elapsed = time.perf_counter() - started
peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
print(json.dumps(dict(seconds=elapsed, peak_rss_bytes=peak if sys.platform == 'darwin'
                     else peak * 1024, returncode=result.returncode, error=result.stderr)))
"""


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(args):
    from hisat2_pipeline import checksum, save

    if args.runs < 1 or args.output_dir.exists():
        raise ValueError("runs must be positive and output folder must be new")
    family = {row["gene_id"] for row in read(args.genes)}
    selected = [row for row in read(args.junctions) if row["gene_id"] in family]
    if not selected:
        raise ValueError("no candidate junction rows in this input")
    groups = defaultdict(list)
    for row in selected:
        groups[(row["sample_id"], row["gene_id"])].append(row)
    output = args.output_dir.resolve()
    output.mkdir(parents=True)
    table = output / "junctions.tsv"
    with table.open("x", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "gene_id", "from", "to", "count"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(selected)
    measurements = []
    for _ in range(args.runs):
        with tempfile.TemporaryDirectory(
            prefix="splice-girl-real-benchmark-"
        ) as folder:
            destination = Path(folder) / "results"
            command = [
                str(args.binary.resolve()),
                "analyze",
                str(table),
                "--output",
                str(destination),
            ]
            completed = subprocess.run(
                [sys.executable, "-c", WORKER, *command],
                check=True,
                capture_output=True,
                text=True,
                timeout=150,
            )
            measurement = json.loads(completed.stdout)
            if measurement["returncode"]:
                raise ValueError(measurement["error"])
            metrics = read(destination / "sample_metrics.tsv")
            if len(metrics) != len(groups) or any(
                row["checks_status"] != "PASS" for row in metrics
            ):
                raise ValueError(
                    "benchmark outputs have missing groups or failed mathematical checks"
                )
            measurement.pop("error")
            measurement["math_checks"] = "PASS"
            measurements.append(measurement)
    observed = {gene for _, gene in groups}
    record = {
        "scope": args.scope,
        "operation": "analyze plus TSV/manifest export; no figures or inference",
        "environment": platform.platform(),
        "inputs": [
            checksum(p)
            for p in [args.junctions, args.genes, args.binary, Path(__file__)]
        ],
        "filtered_input": checksum(table),
        "candidate_genes": len(family),
        "observed_genes": len(observed),
        "unobserved_genes": len(family - observed),
        "junction_rows": len(selected),
        "sample_gene_groups": len(groups),
        "largest_graph_nodes": max(
            len({node for r in rows for node in (r["from"], r["to"])})
            for rows in groups.values()
        ),
        "largest_graph_edges": max(map(len, groups.values())),
        "runs": measurements,
        "median_seconds": statistics.median(r["seconds"] for r in measurements),
        "maximum_peak_rss_bytes": max(r["peak_rss_bytes"] for r in measurements),
        "limits": "One deterministic control prefix; observed genes only; startup/filesystem included; not a full-cohort benchmark.",
    }
    save(output / "benchmark.json", record)
    print(
        json.dumps(
            {
                key: record[key]
                for key in [
                    "observed_genes",
                    "junction_rows",
                    "median_seconds",
                    "maximum_peak_rss_bytes",
                ]
            },
            indent=2,
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["junctions", "genes", "output-dir"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument(
        "--scope",
        required=True,
        help="describe the real source and subset, not just its filename",
    )
    parser.add_argument("--runs", type=int, default=5)
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
