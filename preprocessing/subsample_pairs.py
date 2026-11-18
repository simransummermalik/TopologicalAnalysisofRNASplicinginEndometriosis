#!/usr/bin/env python3
"""Extract the first N synchronized FASTQ pairs for a bounded technical pilot.

This is a deterministic prefix, not a random or representative study sample.
The original files remain unchanged. Input/output hashes and scope are saved.
"""

from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
from contextlib import ExitStack
from pathlib import Path

from hisat2_pipeline import checksum, save


def fastq_record(handle):
    header = handle.readline()
    if not header:
        return None
    sequence, separator, quality = [handle.readline() for _ in range(3)]
    if not header.startswith(b"@") or not separator.startswith(b"+"):
        raise ValueError("malformed four-line FASTQ record")
    sequence, quality = sequence.rstrip(b"\r\n"), quality.rstrip(b"\r\n")
    if not sequence or len(sequence) != len(quality):
        raise ValueError("FASTQ sequence/quality length mismatch or truncated record")
    name = header.split()[0][1:]
    if name.endswith((b"/1", b"/2")):
        name = name[:-2]
    return (
        name,
        header.rstrip(b"\r\n")
        + b"\n"
        + sequence
        + b"\n"
        + separator.rstrip(b"\r\n")
        + b"\n"
        + quality
        + b"\n",
        len(sequence),
    )


def extract(read1, read2, output, pairs, expected1=None, expected2=None):
    if pairs < 1:
        raise ValueError("pairs must be positive")
    if read1.resolve() == read2.resolve():
        raise ValueError("paired inputs must be different files")
    if not read1.is_file() or not read2.is_file():
        raise ValueError("both FASTQ inputs must exist")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()  # Existing results are never overwritten.
    manifest = output / "subset_manifest.json"
    record = {
        "status": "running",
        "method": "first N synchronized cleaned read pairs",
        "requested_pairs": pairs,
        "written_pairs": 0,
        "scope": "technical pilot only; deterministic prefix may have ordering bias; not the full sample or a biological replicate",
        "inputs": [],
        "outputs": [],
    }
    save(manifest, record)
    try:
        for path, expected in [(read1, expected1), (read2, expected2)]:
            item = checksum(path)
            if expected is not None and item["sha256"] != expected:
                raise ValueError(f"source checksum mismatch: {path}")
            item["expected_sha256"] = expected
            record["inputs"].append(item)
        with ExitStack() as stack:
            readers = [
                stack.enter_context(
                    gzip.open(path, "rb") if path.suffix == ".gz" else path.open("rb")
                )
                for path in (read1, read2)
            ]
            targets = [output / "R1.fastq.gz", output / "R2.fastq.gz"]
            writers = [
                stack.enter_context(
                    gzip.GzipFile(filename=str(path), mode="wb", mtime=0)
                )
                for path in targets
            ]
            lengths = [None, None]
            for index in range(pairs):
                left, right = [fastq_record(reader) for reader in readers]
                if left is None or right is None:
                    raise ValueError(
                        f"input ended before requested pair {index + 1}; partial outputs are not complete"
                    )
                if left[0] != right[0]:
                    raise ValueError(
                        f"read names disagree at pair {index + 1}; mates are not synchronized"
                    )
                for mate, value in enumerate((left, right)):
                    writers[mate].write(value[1])
                    lengths[0] = (
                        value[2] if lengths[0] is None else min(lengths[0], value[2])
                    )
                    lengths[1] = (
                        value[2] if lengths[1] is None else max(lengths[1], value[2])
                    )
                record["written_pairs"] = index + 1
            record["read_length_range"] = lengths
        record["outputs"] = [checksum(path) for path in targets]
        record["status"] = "complete"
    except Exception as error:
        record.update(status="failed", error=str(error))
        raise
    finally:
        save(manifest, record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--pairs", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-sha256-r1")
    parser.add_argument("--expected-sha256-r2")
    args = parser.parse_args()
    try:
        record = extract(
            args.read1,
            args.read2,
            args.output_dir,
            args.pairs,
            args.expected_sha256_r1,
            args.expected_sha256_r2,
        )
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    print(
        f"Saved {record['written_pairs']} technical-pilot read pairs to {args.output_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
