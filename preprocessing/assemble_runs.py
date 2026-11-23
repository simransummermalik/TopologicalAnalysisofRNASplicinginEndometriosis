#!/usr/bin/env python3
"""Combine separately completed Cleaver/HISAT2 sample runs into one cohort bundle.

Use this when samples were processed as separate jobs (for example one Slurm
job per sample) instead of one sequential batch. It writes the same
batch_manifest.json, junctions.tsv and samples.tsv that the batch command
writes, so the configured Rust analysis can verify them. Nothing is realigned
or recounted, and the run folders are not modified or copied.

Required TSV columns: sample_id, condition, run_dir, strand_check. run_dir is
a completed run folder (run_manifest.json, junctions.tsv). strand_check may be
blank when the run manifest already records a strand check (runs from
2026-11-23 on); otherwise it is the JSON written by strand_orientation.py for
that sample's BAM. Paths are relative to the sheet.

Every sample must have completed every stage with the given strandedness,
anchor, support and duplicate settings, and all samples must share one gene
map, splice-site file and index. A sample held by the strand check is left
out of the combined table and listed under excluded_samples.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

from batch_preprocessing import write_combined
from hisat2_pipeline import checksum, save
from strand_orientation import EXPECTED_READ1, summarize

PIPELINE = "Splice Girl assembled Cleaver/HISAT2 runs"
RUN_PIPELINE = "Splice Girl Cleaver/HISAT2"
COUNT_UNIT = "primary uniquely mapped read alignments with NH:i:1 and XS strand"
COLUMNS = ("sample_id", "condition", "run_dir", "strand_check")


def load_sheet(path: Path) -> list[dict]:
    rows, seen, dirs = [], set(), set()
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        if len(set(header)) != len(header) or not set(COLUMNS).issubset(header):
            raise ValueError(f"run sheet needs unique columns: {', '.join(COLUMNS)}")
        for number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"run sheet line {number}: wrong number of fields")
            sample = row["sample_id"]
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", sample) or sample in seen:
                raise ValueError(f"line {number}: invalid or duplicate sample_id {sample!r}")
            if row["condition"] not in {"case", "control"}:
                raise ValueError(f"line {number}: condition must be case or control")
            base = path.resolve().parent
            run_dir = (base / row["run_dir"]).resolve()
            if not row["run_dir"].strip() or not run_dir.is_dir():
                raise ValueError(f"line {number}: run_dir is not a folder: {run_dir}")
            if run_dir in dirs:
                raise ValueError(f"line {number}: run_dir reused: {run_dir}")
            strand = row["strand_check"].strip()
            strand_path = (base / strand).resolve() if strand else None
            if strand_path and not strand_path.is_file():
                raise ValueError(f"line {number}: missing strand_check: {strand_path}")
            seen.add(sample)
            dirs.add(run_dir)
            rows.append(
                {
                    "sample_id": sample,
                    "condition": row["condition"],
                    "run_dir": run_dir,
                    "strand_check": strand_path,
                }
            )
    if not rows:
        raise ValueError("run sheet has no samples")
    return rows


def command_value(commands: list[list[str]], flag: str) -> str:
    values = [
        command[index + 1]
        for command in commands
        for index, argument in enumerate(command[:-1])
        if argument == flag
    ]
    if len(values) != 1:
        raise ValueError(f"expected one recorded {flag} argument, found {len(values)}")
    return values[0]


def reference_identity(run: dict) -> dict:
    """Hashes of the shared reference files one run recorded."""
    inputs = {record["path"]: record["sha256"] for record in run.get("inputs", [])}
    commands = run.get("commands", [])
    found = {}
    for name, flag in (("gene_map", "--gene-map"), ("splice_sites", "--known-splicesite-infile")):
        path = command_value(commands, flag)
        if path not in inputs:
            raise ValueError(f"no recorded checksum for {name} {path}")
        found[name] = inputs[path]
    index = sorted(
        sha for path, sha in inputs.items() if path.endswith((".ht2", ".ht2l"))
    )
    if len(index) != 8:
        raise ValueError(f"expected 8 recorded index files, found {len(index)}")
    found["index"] = index
    return found


def check_run(sample: dict, settings: dict) -> dict:
    """Validate one finished run; return its manifest."""
    sample_id = sample["sample_id"]
    manifest_path = sample["run_dir"] / "run_manifest.json"
    if not manifest_path.is_file():
        raise ValueError(f"{sample_id}: missing {manifest_path}")
    run = json.loads(manifest_path.read_text(encoding="utf-8"))
    stages = run.get("stages", [])
    if (
        run.get("pipeline") != RUN_PIPELINE
        or run.get("status") != "complete"
        or not stages
        or any(stage.get("status") != "complete" for stage in stages)
    ):
        raise ValueError(f"{sample_id}: run is not a completed {RUN_PIPELINE} run")
    recorded = run.get("settings", {})
    for key, expected in settings.items():
        if recorded.get(key) != expected:
            raise ValueError(
                f"{sample_id}: {key} is {recorded.get(key)!r}, expected {expected!r}"
            )
    counts = sample["run_dir"] / "junctions.tsv"
    outputs = [r for r in run.get("outputs", []) if Path(r["path"]).name == "junctions.tsv"]
    if len(outputs) != 1 or not counts.is_file():
        raise ValueError(f"{sample_id}: missing recorded junctions.tsv")
    if checksum(counts)["sha256"] != outputs[0]["sha256"]:
        raise ValueError(f"{sample_id}: junctions.tsv changed since the run finished")
    return run


def strand_result(sample: dict, run: dict, strandedness: str) -> dict:
    """The sample's strand call, recomputed from counts, with its source."""
    sample_id = sample["sample_id"]
    if run.get("strand_check"):
        counts, source = run["strand_check"], {"source": "run_manifest"}
    elif sample["strand_check"]:
        counts = json.loads(sample["strand_check"].read_text(encoding="utf-8"))
        if counts.get("sample_id") not in (None, sample_id):
            raise ValueError(
                f"{sample_id}: strand_check file is for {counts.get('sample_id')}"
            )
        source = {"source": "strand_orientation.py", "file": checksum(sample["strand_check"])}
    else:
        raise ValueError(
            f"{sample_id}: no strand check; run strand_orientation.py on its BAM "
            "and add the JSON to the strand_check column"
        )
    if counts.get("expected") not in (None, strandedness):
        raise ValueError(
            f"{sample_id}: strand check expected {counts['expected']}, not {strandedness}"
        )
    # Recompute the call from the counts rather than trusting a stored status.
    return {**summarize(Counter(counts), strandedness), **source}


def assemble(args) -> dict:
    if args.strandedness not in EXPECTED_READ1:
        raise ValueError("assemble requires RF or FR strandedness")
    sheet = load_sheet(args.samples)
    output = args.output_dir.resolve()
    if output.exists():
        raise ValueError(f"output folder already exists: {output}; choose a new folder")
    settings = {
        "strandedness": args.strandedness,
        "minimum_anchor": args.min_anchor,
        "minimum_unique_reads": args.min_unique,
        "exclude_duplicates": args.exclude_duplicates,
        "count_unit": COUNT_UNIT,
    }
    # Check every sample before writing anything.
    checked, reference = [], None
    for sample in sheet:
        run = check_run(sample, {**settings, "sample_id": sample["sample_id"]})
        identity = reference_identity(run)
        if reference is None:
            reference = identity
        elif identity != reference:
            different = [key for key in identity if identity[key] != reference[key]]
            raise ValueError(
                f"{sample['sample_id']}: {', '.join(different)} differ from "
                f"{sheet[0]['sample_id']}"
            )
        checked.append((sample, strand_result(sample, run, args.strandedness)))

    items, kept = [], []
    for sample, strand in checked:
        item = {
            "sample_id": sample["sample_id"],
            "condition": sample["condition"],
            "status": "complete",
            "run_manifest": checksum(sample["run_dir"] / "run_manifest.json"),
            "strand_check": strand,
        }
        if strand["status"] != "pass":
            item["status"] = "excluded"
            item["exclusion_reason"] = f"strand_check_{strand['status']}"
        else:
            kept.append(
                (sample["sample_id"], sample["condition"], sample["run_dir"] / "junctions.tsv")
            )
        items.append(item)
    if not kept:
        raise ValueError("every sample was excluded by the strand check")

    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()
    record = {
        "pipeline": PIPELINE,
        "status": "running",
        "sample_sheet": checksum(args.samples),
        "settings": settings,
        "reference_identity": reference,
        "samples": [item for item in items if item["status"] == "complete"],
        "excluded_samples": [item for item in items if item["status"] == "excluded"],
        "outputs": [],
    }
    manifest = output / "batch_manifest.json"
    save(manifest, record)
    try:
        write_combined(output, kept)
        record["outputs"] = [
            checksum(output / name) for name in ("junctions.tsv", "samples.tsv")
        ]
        record["status"] = "complete"
    except BaseException as error:
        record.update(status="failed", error=str(error) or type(error).__name__)
        raise
    finally:
        save(manifest, record)
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--samples", type=Path, required=True, help="run sheet TSV")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--strandedness", choices=sorted(EXPECTED_READ1), required=True)
    parser.add_argument("--min-anchor", type=int, required=True)
    parser.add_argument("--min-unique", type=int, required=True)
    parser.add_argument("--exclude-duplicates", action="store_true")
    args = parser.parse_args(argv)
    try:
        record = assemble(args)
    except (OSError, ValueError, KeyError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    output = args.output_dir
    print(f"Assembled {len(record['samples'])} samples: {output / 'junctions.tsv'} and {output / 'samples.tsv'}")
    for item in record["excluded_samples"]:
        check = item["strand_check"]
        agreement = check.get("read1_agreement")
        shown = "no read-1 records" if agreement is None else f"read 1 agreement {agreement:.2%}"
        print(f"WARNING: {item['sample_id']} excluded by the strand check ({shown}).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
