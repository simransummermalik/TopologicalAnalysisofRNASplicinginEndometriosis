#!/usr/bin/env python3
"""Prepare references or plan/run Cleaver -> HISAT2 -> annotated junctions.

This is the production-method adapter; run_preprocessing.py retains the legacy
STAR prototype. Thresholds and RNA strandedness must be explicitly supplied.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from build_junction_gene_map import build_gene_map

HERE = Path(__file__).resolve().parent


def checksum(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {
        "path": str(path.resolve()),
        "sha256": digest.hexdigest(),
        "bytes": path.stat().st_size,
    }


def save(path: Path, record):
    temporary = path.with_suffix(".json.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
    temporary.replace(path)


def require_tools(names):
    missing = [name for name in names if not shutil.which(name)]
    if missing:
        raise ValueError(
            "missing tools: " + ", ".join(missing) + "; see preprocessing/README.md"
        )


def versions(names):
    result = {}
    for name in names:
        call = subprocess.run(
            [name, "--version"], capture_output=True, text=True, check=False, timeout=15
        )
        lines = (call.stdout + call.stderr).strip().splitlines()
        result[name] = {
            "path": shutil.which(name),
            "version_output": "\n".join(lines[:12]),
        }
    return result


def index_files(prefix: Path):
    for suffix in ("ht2", "ht2l"):
        files = [Path(f"{prefix}.{index}.{suffix}") for index in range(1, 9)]
        if all(path.is_file() for path in files):
            return files
    raise ValueError(
        f"incomplete HISAT2 index prefix {prefix}; expected all eight .ht2 or .ht2l files"
    )


def splice_hints(gene_map: Path, output: Path):
    sites = set()
    with gene_map.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            # HISAT2 uses 0-based coordinates of the two flanking exon bases.
            sites.add(
                (
                    row["chromosome"],
                    int(row["intron_start"]) - 2,
                    int(row["intron_end"]),
                    {"1": "+", "2": "-"}[row["strand"]],
                )
            )
    with output.open("x", encoding="utf-8") as handle:
        for site in sorted(sites):
            handle.write("\t".join(map(str, site)) + "\n")


def reference(args):
    for path in (args.fasta, args.gtf):
        if not path.is_file():
            raise ValueError(f"missing reference input: {path}")
    if args.threads < 1:
        raise ValueError("threads must be positive")
    output = args.output_dir.resolve()
    fasta = output / "reference.fa"
    prefix = output / "genome"
    command = ["hisat2-build", "-p", str(args.threads), str(fasta), str(prefix)]
    if args.command == "plan-reference":
        print(
            f"Create new folder {output}; expand/copy FASTA, build exon map and 0-based splice hints, then:"
        )
        print(shlex.join(command))
        return
    require_tools(["hisat2-build"])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()
    record = {
        "status": "running",
        "method": "HISAT2 reference preparation",
        "inputs": [],
        "commands": [command],
        "tools": versions(["hisat2-build"]),
    }
    manifest = output / "reference_manifest.json"
    save(manifest, record)
    try:
        record["inputs"] = [checksum(args.fasta), checksum(args.gtf)]
        opener = gzip.open if args.fasta.suffix == ".gz" else open
        with opener(args.fasta, "rb") as source, fasta.open("xb") as dest:
            shutil.copyfileobj(source, dest)
        map_path = output / "junction_exon_map.tsv"
        record["annotation_summary"] = build_gene_map(args.gtf, map_path)
        splice_hints(map_path, output / "splice_sites.tsv")
        with (output / "index.log").open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        record["outputs"] = [
            checksum(path)
            for path in [*index_files(prefix), map_path, output / "splice_sites.tsv"]
        ]
        record["status"] = "complete"
    except Exception as error:
        record.update(status="failed", error=str(error))
        raise
    finally:
        save(manifest, record)
    print(f"Prepared HISAT2 reference: {prefix}")


def configuration(args):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", args.sample_id):
        raise ValueError(
            "sample ID must start with a letter/number and use only letters, numbers, dot, hyphen, underscore"
        )
    for path in [args.read1, args.gene_map, args.splice_sites] + (
        [args.read2] if args.read2 else []
    ):
        if not path.is_file():
            raise ValueError(f"input file does not exist: {path}")
    if args.read2 and args.read1.resolve() == args.read2.resolve():
        raise ValueError("paired reads must be two different files")
    if min(args.threads, args.min_unique, args.min_anchor) < 1:
        raise ValueError("threads, min-unique, and min-anchor must be positive")
    index_files(args.index)
    if args.read2 and args.strandedness not in {"unstranded", "FR", "RF"}:
        raise ValueError("paired reads require unstranded, FR, or RF strandedness")
    if not args.read2 and args.strandedness not in {"unstranded", "F", "R"}:
        raise ValueError("single-end reads require unstranded, F, or R strandedness")
    output = args.output_dir.resolve()
    cleaned1 = (
        args.read1.resolve() if args.skip_cleaning else output / "trimmed.R1.fastq.gz"
    )
    cleaned2 = (
        (args.read2.resolve() if args.skip_cleaning else output / "trimmed.R2.fastq.gz")
        if args.read2
        else None
    )
    commands = []
    if not args.skip_cleaning:
        clean = [
            "cleaver",
            "fastp",
            "-i",
            str(args.read1.resolve()),
            "-o",
            str(cleaned1),
            "-j",
            str(output / "fastp.json"),
        ]
        if args.read2:
            clean.extend(["-I", str(args.read2.resolve()), "-O", str(cleaned2), "-2"])
        commands.append(clean)
    align = [
        "hisat2",
        "-x",
        str(args.index.resolve()),
        "-p",
        str(args.threads),
        "--known-splicesite-infile",
        str(args.splice_sites.resolve()),
        "--summary-file",
        str(output / "alignment_summary.txt"),
        "-S",
        str(output / "aligned.sam"),
    ]
    if args.strandedness != "unstranded":
        align.extend(["--rna-strandness", args.strandedness])
    align.extend(
        ["-1", str(cleaned1), "-2", str(cleaned2)]
        if cleaned2
        else ["-U", str(cleaned1)]
    )
    commands.append(align)
    convert = [
        sys.executable,
        str(HERE / "sam_junctions.py"),
        "--sam",
        str(output / "aligned.sam"),
        "--gene-map",
        str(args.gene_map.resolve()),
        "--sample-id",
        args.sample_id,
        "--output",
        str(output / "junctions.tsv"),
        "--qc-json",
        str(output / "junction_qc.json"),
        "--min-anchor",
        str(args.min_anchor),
        "--min-unique",
        str(args.min_unique),
    ]
    if args.exclude_duplicates:
        convert.append("--exclude-duplicates")
    commands.extend(
        [
            convert,
            [
                "cleaver",
                "samtools",
                "sort",
                "-o",
                str(output / "aligned.bam"),
                str(output / "aligned.sam"),
            ],
            ["cleaver", "samtools", "index", str(output / "aligned.bam")],
        ]
    )
    return output, commands


def run(args):
    output, commands = configuration(args)
    if args.command == "plan":
        for command in commands:
            print(shlex.join(command))
        return
    require_tools(["cleaver", "hisat2"])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()
    record = {
        "pipeline": "Splice Girl Cleaver/HISAT2",
        "status": "running",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "commands": commands,
        "settings": {
            "sample_id": args.sample_id,
            "strandedness": args.strandedness,
            "minimum_unique_reads": args.min_unique,
            "minimum_anchor": args.min_anchor,
            "exclude_duplicates": args.exclude_duplicates,
            "cleaning_skipped": args.skip_cleaning,
            "count_unit": "primary uniquely mapped read alignments with NH:i:1 and XS strand",
        },
        "tools": versions(["cleaver", "hisat2"]),
    }
    manifest = output / "run_manifest.json"
    save(manifest, record)
    try:
        paths = [
            args.read1,
            args.gene_map,
            args.splice_sites,
            *index_files(args.index),
        ] + ([args.read2] if args.read2 else [])
        record["inputs"] = [checksum(path) for path in paths]
        for index, command in enumerate(commands, 1):
            print(f"[{index}/{len(commands)}] {shlex.join(command)}", flush=True)
            with (output / f"step_{index}.log").open("w") as log:
                subprocess.run(
                    command, stdout=log, stderr=subprocess.STDOUT, check=True
                )
        required = [
            output / "junctions.tsv",
            output / "junction_qc.json",
            output / "aligned.bam",
            output / "aligned.bam.bai",
        ]
        for path in required:
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"missing or empty expected output: {path}")
        record["outputs"] = [checksum(path) for path in required]
        record["status"] = "complete"
    except Exception as error:
        record.update(status="failed", error=str(error))
        raise
    finally:
        save(manifest, record)
    print(f"Splice Girl junction input: {output / 'junctions.tsv'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("check")
    for action in ("prepare-reference", "plan-reference"):
        ref = sub.add_parser(action)
        ref.add_argument("--fasta", type=Path, required=True)
        ref.add_argument("--gtf", type=Path, required=True)
        ref.add_argument("--output-dir", type=Path, required=True)
        ref.add_argument("--threads", type=int, default=2)
    for action in ("plan", "run"):
        item = sub.add_parser(action)
        item.add_argument("--sample-id", required=True)
        for name in ("read1", "index", "gene-map", "splice-sites", "output-dir"):
            item.add_argument(f"--{name}", type=Path, required=True)
        item.add_argument("--read2", type=Path)
        item.add_argument(
            "--strandedness",
            choices=["unstranded", "F", "R", "FR", "RF"],
            required=True,
        )
        item.add_argument("--min-anchor", type=int, required=True)
        item.add_argument("--min-unique", type=int, required=True)
        item.add_argument("--threads", type=int, default=2)
        item.add_argument("--skip-cleaning", action="store_true")
        item.add_argument("--exclude-duplicates", action="store_true")
    args = parser.parse_args()
    try:
        if args.command == "check":
            require_tools(["cleaver", "hisat2", "hisat2-build"])
            print(json.dumps(versions(["cleaver", "hisat2", "hisat2-build"]), indent=2))
        elif args.command in {"prepare-reference", "plan-reference"}:
            reference(args)
        else:
            run(args)
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
