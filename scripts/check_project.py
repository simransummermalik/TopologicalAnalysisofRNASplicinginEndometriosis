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


def smoke(pairs=1):
    """Actual external-tool test on synthetic reads, never on the human cohort."""
    sys.path.insert(0, str(ROOT / "preprocessing"))
    from hisat2_pipeline import require_tools

    if pairs < 1:
        raise ValueError("smoke-test pair count must be positive")
    require_tools(["cleaver", "hisat2", "hisat2-build", "samtools"])
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
            # A reverse-complement mate spans the same synthetic junction.
            # Both mates count as supporting read alignments, not one fragment.
            read = (fixtures / "tiny_reads_R1.fastq").read_text().splitlines()
            reverse = read[1].translate(str.maketrans("ACGT", "TGCA"))[::-1]
            with (
                (root / f"{name}.fastq").open("x") as mate1,
                (root / f"{name}_R2.fastq").open("x") as mate2,
            ):
                for number in range(pairs):
                    mate1.write(f"@smoke_{number}/1\n{read[1]}\n+\n{read[3]}\n")
                    mate2.write(f"@smoke_{number}/2\n{reverse}\n+\n{read[3][::-1]}\n")
        sheet = root / "input.tsv"
        sheet.write_text(
            "sample_id\tcondition\tread1\tread2\nCASE\tcase\tcase.fastq\tcase_R2.fastq\nCONTROL\tcontrol\tcontrol.fastq\tcontrol_R2.fastq\n"
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
        assert len(rows) == 2 and all(int(row["count"]) == 2 * pairs for row in rows), (
            rows
        )
        assert (
            json.loads((output / "batch_manifest.json").read_text())["status"]
            == "complete"
        )
        assert (root / "report/report.html").is_file()
        # Exercise actual pipes, BAM decoding and Cleaver indexing. Counts and
        # every QC counter must agree with the retained-SAM workflow above.
        streamed = root / "streamed"
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
                "--storage-mode",
                "stream-bam",
                "--sort-memory-mb",
                "1",
                "--output-dir",
                streamed,
            ]
        )
        assert (streamed / "junctions.tsv").read_bytes() == (
            output / "junctions.tsv"
        ).read_bytes()
        for name in ("CASE", "CONTROL"):
            assert not (streamed / name / "aligned.sam").exists()
            assert not (streamed / name / "aligned.partial.bam").exists()
            assert (streamed / name / "aligned.bam.bai").is_file()
            legacy_qc = json.loads((output / name / "junction_qc.json").read_text())
            stream_qc = json.loads((streamed / name / "junction_qc.json").read_text())
            assert stream_qc == legacy_qc, (stream_qc, legacy_qc)
            record = json.loads((streamed / name / "run_manifest.json").read_text())
            assert (
                record["status"] == "complete"
                and not record["settings"]["sam_retained"]
            )
            assert all(stage["status"] == "complete" for stage in record["stages"])
            if pairs >= 12000:
                import re

                sort_stage = next(
                    stage for stage in record["stages"] if stage["name"] == "align_sort"
                )
                sort_log = (streamed / name / sort_stage["logs"][-1]).read_text()
                assert re.search(r"merging from [1-9][0-9]* files", sort_log), sort_log
                print(
                    f"Actual sort spill confirmed for {name} ({pairs} synthetic pairs)."
                )
        call(
            [
                "cargo",
                "run",
                "--quiet",
                "--",
                "report",
                streamed / "junctions.tsv",
                streamed / "samples.tsv",
                "--output",
                root / "streamed_report",
            ]
        )
        assert (root / "streamed_report/report.html").is_file()
    print(
        "External-tool synthetic smoke test passed in retain-sam and stream-bam modes; counts and QC agree."
    )


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
        [
            sys.executable,
            "-m",
            "unittest",
            "discover",
            "-s",
            "scripts",
            "-p",
            "test_*.py",
        ],
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
