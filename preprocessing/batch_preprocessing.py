#!/usr/bin/env python3
"""Plan or process a sample sheet sequentially using the Cleaver/HISAT2 adapter.

Required TSV columns: sample_id, condition, read1, read2. Paths are relative to
the sample sheet; read2 may be blank for single-end input. No downloads occur.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shlex
import subprocess
import sys
from pathlib import Path

from hisat2_pipeline import checksum, configuration, require_tools, run, save


def load_samples(path: Path):
    samples = []
    seen = set()
    used_reads = set()
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        if len(set(header)) != len(header) or not {"sample_id", "condition", "read1", "read2"}.issubset(header):
            raise ValueError("sample sheet needs unique sample_id, condition, read1, read2 columns")
        for number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"sample sheet line {number}: wrong number of fields")
            sample = row["sample_id"]
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", sample) or sample in seen:
                raise ValueError(f"line {number}: invalid or duplicate sample_id {sample!r}")
            if row["condition"] not in {"case", "control"}:
                raise ValueError(f"line {number}: condition must be case or control")
            if not row["read1"].strip():
                raise ValueError(f"line {number}: read1 is required")
            item = {"sample_id": sample, "condition": row["condition"]}
            for name in ("read1", "read2"):
                value = row[name]
                resolved = (path.resolve().parent / value).resolve() if value else None
                if resolved:
                    if not resolved.is_file():
                        raise ValueError(f"line {number}: missing {name}: {resolved}")
                    if resolved in used_reads:
                        raise ValueError(f"line {number}: read file reused across samples/mates: {resolved}")
                    used_reads.add(resolved)
                item[name] = resolved
            seen.add(sample)
            samples.append(item)
    if not samples:
        raise ValueError("sample sheet has no samples")
    return samples


def batch(args):
    samples = load_samples(args.samples)
    output = args.output_dir.resolve()
    if output.exists():
        raise ValueError(f"output folder already exists: {output}; choose a new folder")
    jobs = []
    for sample in samples:
        options = vars(args).copy()
        options.update(sample)
        options.update(output_dir=output / sample["sample_id"], command="plan" if args.command == "plan" else "run")
        job = argparse.Namespace(**options)
        _, commands = configuration(job)  # Validate every sample before any work.
        jobs.append(job)
        if args.command == "plan":
            print(f"# {job.sample_id} ({sample['condition']})")
            for command in commands:
                print(shlex.join(command))
    if args.command == "plan":
        return
    require_tools(["cleaver", "hisat2"])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()
    record = {
        "pipeline": "Splice Girl sequential Cleaver/HISAT2 batch",
        "status": "running",
        "sample_sheet": checksum(args.samples),
        "samples": [{"sample_id": s["sample_id"], "condition": s["condition"], "status": "not_started"} for s in samples],
        "outputs": [],
    }
    manifest = output / "batch_manifest.json"
    save(manifest, record)
    try:
        for job, item in zip(jobs, record["samples"]):
            item["status"] = "running"
            save(manifest, record)
            try:
                run(job)
            except Exception:
                item["status"] = "failed"
                raise
            item["status"] = "complete"
            item["run_manifest"] = checksum(job.output_dir / "run_manifest.json")
            save(manifest, record)
        # Publish a combined table only when every job completed. Identity is
        # preserved; absent genes/edges are not filled with fabricated zeros.
        temporary = output / "junctions.tsv.tmp"
        with temporary.open("x", encoding="utf-8", newline="") as dest:
            writer = csv.writer(dest, delimiter="\t", lineterminator="\n")
            header = ["sample_id", "gene_id", "from", "to", "count"]
            writer.writerow(header)
            for job in jobs:
                with (job.output_dir / "junctions.tsv").open(encoding="utf-8", newline="") as source:
                    reader = csv.reader(source, delimiter="\t")
                    if next(reader, None) != header:
                        raise ValueError(f"unexpected junction schema for {job.sample_id}")
                    for row in reader:
                        if len(row) != 5 or row[0] != job.sample_id:
                            raise ValueError(f"junction identity mismatch for {job.sample_id}")
                        writer.writerow(row)
        temporary.replace(output / "junctions.tsv")
        with (output / "samples.tsv").open("x", encoding="utf-8") as handle:
            handle.write("sample_id\tcondition\n")
            for sample in samples:
                handle.write(f"{sample['sample_id']}\t{sample['condition']}\n")
        record["outputs"] = [checksum(output / name) for name in ("junctions.tsv", "samples.tsv")]
        record["status"] = "complete"
    except Exception as error:
        record.update(status="failed", error=str(error))
        raise
    finally:
        save(manifest, record)
    print(f"Ready for Splice Girl: {output / 'junctions.tsv'} and {output / 'samples.tsv'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["plan", "run"])
    for name in ("samples", "index", "gene-map", "splice-sites", "output-dir"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--strandedness", choices=["unstranded", "F", "R", "FR", "RF"], required=True)
    parser.add_argument("--min-anchor", type=int, required=True)
    parser.add_argument("--min-unique", type=int, required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--skip-cleaning", action="store_true")
    parser.add_argument("--exclude-duplicates", action="store_true")
    try:
        batch(parser.parse_args())
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
