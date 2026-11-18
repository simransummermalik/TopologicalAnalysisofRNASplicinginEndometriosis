#!/usr/bin/env python3
"""Read-only cohort storage/input check. Never downloads or launches an alignment.

Estimate the current SAM-retaining workflow using the actual pilot file sizes
and saved FastQC counts. Estimates are planning values, not disk guarantees.
"""

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

from hisat2_pipeline import checksum, save

ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    with path.open(encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pair_counts(paths):
    counts = {}
    for path in paths:
        for row in rows(path):
            sample, mate = row["Sample"].rsplit("_", 1)
            if mate not in {"forward", "reverse"}:
                raise ValueError(f"unrecognized mate in {row['Sample']}")
            # This saved MultiQC export expresses total_sequences in millions.
            value = round(float(row["fastqc-total_sequences"]) * 1_000_000)
            if value <= 0 or mate in counts.setdefault(sample, {}):
                raise ValueError(f"invalid/duplicate QC count for {sample}/{mate}")
            counts[sample][mate] = value
    result = {}
    for sample, mates in counts.items():
        if set(mates) != {"forward", "reverse"} or mates["forward"] != mates["reverse"]:
            raise ValueError(f"unequal or missing mate counts for {sample}")
        result[sample] = mates["forward"]
    return result


def run(args):
    metadata = rows(args.metadata)
    expected = {row["sample_id"]: row["condition"] for row in metadata}
    if len(expected) != len(metadata) or not expected:
        raise ValueError("metadata needs unique samples")
    counts = pair_counts(args.qc)
    if counts.keys() != expected.keys():
        raise ValueError("QC sample IDs do not exactly match the cohort metadata")
    sheet = rows(args.samples) if args.samples else []
    supplied = {row["sample_id"]: row for row in sheet}
    if sheet and (len(supplied) != len(sheet) or supplied.keys() != expected.keys()):
        raise ValueError("sample sheet must contain exactly the cohort sample IDs")
    available = []
    used = set()
    for sample, condition in expected.items():
        item = supplied.get(sample)
        if item and item["condition"] != condition:
            raise ValueError(f"condition disagrees with metadata: {sample}")
        present = []
        for mate in ["read1", "read2"]:
            value = item.get(mate, "") if item else ""
            path = (args.samples.resolve().parent / value).resolve() if value else None
            if path and path in used:
                raise ValueError(f"read path reused: {path}")
            if path:
                used.add(path)
            present.append(bool(path and path.is_file() and path.stat().st_size))
        available.append(
            {
                "sample_id": sample,
                "condition": condition,
                "raw_pairs_from_saved_qc": counts[sample],
                "read1_present": present[0],
                "read2_present": present[1],
            }
        )
    pilot = json.loads(args.pilot_manifest.read_text())
    pairs = pilot["written_pairs"]
    clean_bytes = sum(item["bytes"] for item in pilot["outputs"])
    total_pairs = sum(counts.values())
    factor = total_pairs / pairs
    # Retain raw + cleaned reads, SAM and BAM; reserve 25% on top. Cleaning
    # compression is only a proxy for raw storage and is explicitly labeled.
    estimate = int(
        factor
        * (
            2 * clean_bytes
            + args.pilot_sam.stat().st_size
            + args.pilot_bam.stat().st_size
        )
        * 1.25
    )
    free = None
    if args.storage:
        if not args.storage.is_dir():
            raise ValueError("storage must be an existing mounted/accessible folder")
        free = shutil.disk_usage(args.storage).free
    record = {
        "scope": "full eight-sample input/storage preflight only; no processing launched",
        "samples": available,
        "input_paths_checked": bool(sheet),
        "raw_read_pairs": total_pairs,
        "estimated_working_bytes_with_25_percent_margin": estimate,
        "storage_path": str(args.storage.resolve()) if args.storage else None,
        "storage_free_bytes": free,
        "all_paired_inputs_present": bool(sheet)
        and all(row["read1_present"] and row["read2_present"] for row in available),
        "storage_estimate_fits": free >= estimate if free is not None else False,
        "inputs": [
            checksum(p)
            for p in [args.metadata, *args.qc, args.pilot_manifest, Path(__file__)]
            + ([args.samples] if args.samples else [])
        ],
        "pilot_file_sizes": {
            "pairs": pairs,
            "clean_gzip_bytes": clean_bytes,
            "sam_bytes": args.pilot_sam.stat().st_size,
            "bam_bytes": args.pilot_bam.stat().st_size,
        },
        "limits": [
            "Raw and cleaned storage use the pilot cleaned-gzip rate as a proxy; allow extra for references, indexes, tool temporary files and other users.",
            "Sample presence does not verify FASTQ integrity, identity, checksums or full-sample scope.",
            "Final strand, mapping-ambiguity, support and coverage choices still need the scientific review recorded in priorities.md.",
            "No path means no assumed drive/server; an SSH address must be checked on that server.",
        ],
    }
    if args.output:
        if args.output.exists():
            raise ValueError("output file exists; choose a new filename")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        save(args.output, record)
    print(json.dumps(record, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metadata", type=Path, default=ROOT / "data/metadata/samples.tsv"
    )
    parser.add_argument(
        "--qc",
        type=Path,
        nargs="+",
        default=[
            ROOT
            / "data/galaxy_data/october_20/Galaxy138-[controls_multiqc_stats].tabular",
            ROOT
            / "data/galaxy_data/october_20/Galaxy140-[endometriosis_multiqc_stats].tabular",
        ],
    )
    parser.add_argument(
        "--samples",
        type=Path,
        help="batch TSV: sample_id, condition, read1, read2; paths relative to this sheet",
    )
    parser.add_argument(
        "--storage", type=Path, help="existing local/mounted output folder"
    )
    base = ROOT / "preprocessing/runs/SRR36710387_pilot_nov17"
    parser.add_argument(
        "--pilot-manifest", type=Path, default=base / "subset/subset_manifest.json"
    )
    parser.add_argument(
        "--pilot-sam", type=Path, default=base / "alignment/aligned.sam"
    )
    parser.add_argument(
        "--pilot-bam", type=Path, default=base / "alignment/aligned.bam"
    )
    parser.add_argument("--output", type=Path)
    try:
        run(parser.parse_args())
    except (OSError, ValueError, KeyError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
