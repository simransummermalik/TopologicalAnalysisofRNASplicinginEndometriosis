#!/usr/bin/env python3
"""Run the project's required local software checks; nothing is committed or pushed."""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def call(command):
    print("\n$ " + " ".join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), cwd=ROOT, check=True)


def smoke():
    """Actual external-tool test on synthetic reads, never on the human cohort."""
    sys.path.insert(0, str(ROOT / "preprocessing"))
    from hisat2_pipeline import require_tools

    require_tools(["cleaver", "hisat2", "hisat2-build"])
    with tempfile.TemporaryDirectory(prefix="splice-girl-smoke-") as directory:
        root = Path(directory)
        reference = root / "reference"
        fixtures = ROOT / "preprocessing/tests/fixtures"
        call(
            [
                sys.executable,
                "preprocessing/hisat2_pipeline.py",
                "prepare-reference",
                "--fasta",
                fixtures / "tiny_genome.fa",
                "--gtf",
                fixtures / "tiny_annotation.gtf",
                "--output-dir",
                reference,
            ]
        )
        # Two distinct copies demonstrate sample-preserving batch execution.
        for name in ("case", "control"):
            (root / f"{name}.fastq").write_bytes(
                (fixtures / "tiny_reads_R1.fastq").read_bytes()
            )
        sheet = root / "input.tsv"
        sheet.write_text(
            "sample_id\tcondition\tread1\tread2\nCASE\tcase\tcase.fastq\t\nCONTROL\tcontrol\tcontrol.fastq\t\n"
        )
        output = root / "batch"
        call(
            [
                sys.executable,
                "preprocessing/batch_preprocessing.py",
                "run",
                "--samples",
                sheet,
                "--index",
                reference / "genome",
                "--gene-map",
                reference / "junction_exon_map.tsv",
                "--splice-sites",
                reference / "splice_sites.tsv",
                "--strandedness",
                "unstranded",
                "--min-anchor",
                "8",
                "--min-unique",
                "1",
                "--output-dir",
                output,
            ]
        )
        call(
            [
                "cargo",
                "run",
                "--quiet",
                "--",
                "report",
                output / "junctions.tsv",
                output / "samples.tsv",
                "--output",
                root / "report",
            ]
        )
        import csv
        import json

        with (output / "junctions.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert {row["sample_id"] for row in rows} == {"CASE", "CONTROL"}
        assert len(rows) == 2 and all(row["count"] == "1" for row in rows), rows
        assert (
            json.loads((output / "batch_manifest.json").read_text())["status"]
            == "complete"
        )
        assert (root / "report/report.html").is_file()
    print("External-tool synthetic smoke test passed.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--with-tools",
        action="store_true",
        help="also run Cleaver/HISAT2 on a tiny synthetic fixture",
    )
    args = parser.parse_args()
    commands = [
        ["cargo", "fmt", "--all", "--", "--check"],
        ["cargo", "test", "--all-targets"],
        ["cargo", "clippy", "--all-targets", "--", "-D", "warnings"],
        ["cargo", "build", "--release"],
        [sys.executable, "-m", "unittest", "discover", "-s", "preprocessing/tests"],
        [sys.executable, "-m", "unittest", "discover", "-s", "gene_sets/tests"],
        [sys.executable, "scripts/benchmark_candidate.py", "--runs", "3"],
    ]
    try:
        for command in commands:
            call(command)
        if args.with_tools:
            smoke()
    except (OSError, subprocess.SubprocessError, ValueError, AssertionError) as error:
        print(f"Project check failed: {error}", file=sys.stderr)
        return 1
    print("\nAll requested project checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
