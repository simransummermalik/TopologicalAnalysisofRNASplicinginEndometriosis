#!/usr/bin/env python3
"""Benchmark a deterministic candidate-sized synthetic Splice Girl input.

The generated table is deliberately synthetic: it exercises the software at a
larger, candidate-like scale without pretending to be real RNA-seq evidence.
The generated files live in a temporary directory and are removed when the
benchmark finishes.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def build_input(path: Path, genes: int, samples: int) -> None:
    rows = []
    for sample_index in range(samples):
        sample = f"sample_{sample_index + 1:02d}"
        for gene_index in range(genes):
            gene = f"GENE_{gene_index + 1:04d}"
            base = 100 + (gene_index % 17)
            rows.extend(
                [
                    (sample, gene, "E1", "E2", base),
                    (sample, gene, "E2", "E3", base - 11),
                    (sample, gene, "E1", "E3", 5 + sample_index % 3),
                ]
            )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "gene_id", "from", "to", "count"])
        writer.writerows(rows)


def build_manifest(path: Path, samples: int) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample_id",
                "condition",
                "subject_id",
                "tissue",
                "accession",
                "junction_path",
            ]
        )
        for sample_index in range(samples):
            sample = f"sample_{sample_index + 1:02d}"
            condition = "control" if sample_index < samples // 2 else "case"
            writer.writerow([sample, condition, sample, "synthetic", "NA", "NA"])


def command(binary: Path, junctions: Path, manifest: Path) -> list[str]:
    return [str(binary), "compare", str(junctions), str(manifest)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", type=int, default=200)
    parser.add_argument("--samples", type=int, default=8)
    parser.add_argument("--runs", type=int, default=10)
    parser.add_argument("--binary", type=Path, default=REPO / "target" / "release" / "splice-girl")
    args = parser.parse_args()
    if args.genes < 1 or args.samples < 2 or args.samples % 2 or args.runs < 1:
        parser.error("genes and runs must be positive; samples must be an even number of at least 2")
    if not args.binary.is_file():
        raise SystemExit(f"binary does not exist: {args.binary}; run cargo build --release first")

    with tempfile.TemporaryDirectory(prefix="splice-girl-benchmark-") as directory:
        folder = Path(directory)
        junctions = folder / "junctions.tsv"
        manifest = folder / "manifest.tsv"
        build_input(junctions, args.genes, args.samples)
        build_manifest(manifest, args.samples)
        run = command(args.binary, junctions, manifest)
        timings = []
        for _ in range(args.runs):
            started = time.perf_counter()
            completed = subprocess.run(run, capture_output=True, text=True, check=False)
            timings.append((time.perf_counter() - started) * 1000)
            if completed.returncode != 0:
                raise SystemExit(completed.stderr or completed.stdout)

    print(f"input: {args.genes} genes x {args.samples} samples x 3 junctions")
    print(f"rows: {args.genes * args.samples * 3}")
    print(f"runs: {args.runs}")
    print(f"median_ms: {statistics.median(timings):.3f}")
    print(f"mean_ms: {statistics.mean(timings):.3f}")
    print(f"min_ms: {min(timings):.3f}")
    print(f"max_ms: {max(timings):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
