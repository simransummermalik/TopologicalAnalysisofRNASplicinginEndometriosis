#!/usr/bin/env python3
"""Post-alignment strandedness check for spliced, uniquely mapped reads.

HISAT2's XS tag gives the transcript strand of a spliced alignment from the
splice-site motif, independently of which strand the read aligned to. A read
that aligned to the same strand as XS is "sense"; the opposite is
"antisense". In an RF library read 1 is antisense and read 2 is sense; in FR
it is the reverse.

The 2026-11-18 decision holds a sample for review when fewer than 95% of its
eligible read-1 spliced alignments have the expected orientation (antisense
for RF). Eligible records are primary, mapped, not QC-failed or
supplementary, spliced (CIGAR contains N), NH:i:1 and with XS:A:+ or XS:A:-.

The junction counter records these counts in junction_qc.json for every new
run. For an existing BAM:
  samtools view aligned.bam | strand_orientation.py --sam - --expected RF
Prints a JSON summary; exits 0 on pass, 3 when the sample is held for review.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable

MINIMUM = 0.95
EXPECTED_READ1 = {"RF": "antisense", "FR": "sense"}


def add_record(counts: Counter, flag: int, transcript_strand: str) -> None:
    """Count one eligible record; transcript_strand is '+' or '-' from XS."""
    read_strand = "-" if flag & 16 else "+"
    relation = "sense" if read_strand == transcript_strand else "antisense"
    if flag & 1 and flag & 64:
        mate = "read1"
    elif flag & 1 and flag & 128:
        mate = "read2"
    else:
        mate = "unpaired"
    counts[f"{mate}_{relation}"] += 1


def summarize(
    counts: Counter, expected: str | None = None, minimum: float = MINIMUM
) -> dict:
    """Fractions per mate and, for RF/FR, a pass or hold call on read 1."""
    summary: dict = {
        key: counts.get(key, 0)
        for key in (
            "read1_sense",
            "read1_antisense",
            "read2_sense",
            "read2_antisense",
            "unpaired_sense",
            "unpaired_antisense",
        )
    }
    for mate in ("read1", "read2", "unpaired"):
        total = summary[f"{mate}_sense"] + summary[f"{mate}_antisense"]
        summary[f"{mate}_antisense_fraction"] = (
            summary[f"{mate}_antisense"] / total if total else None
        )
    summary["unit"] = (
        "primary, uniquely mapped (NH:i:1) spliced alignments with XS strand"
    )
    if expected is None:
        return summary
    if expected not in EXPECTED_READ1:
        raise ValueError("expected strandedness must be RF or FR")
    antisense = summary["read1_antisense_fraction"]
    agreement = None
    if antisense is not None:
        agreement = antisense if expected == "RF" else 1 - antisense
    summary.update(
        expected=expected,
        minimum_read1_agreement=minimum,
        read1_agreement=agreement,
        status=(
            "no_read1_records"
            if agreement is None
            else "pass"
            if agreement >= minimum
            else "hold_for_review"
        ),
    )
    return summary


def eligible_strand(fields: list[str]) -> tuple[int, str] | None:
    """Return (flag, XS strand) for an eligible SAM record, else None."""
    flag = int(fields[1])
    if flag & (4 | 256 | 512 | 2048) or "N" not in fields[5]:
        return None
    tags = {tag[:2]: tag[5:] for tag in fields[11:] if tag[2:3] == ":"}
    if tags.get("NH") != "1" or tags.get("XS") not in {"+", "-"}:
        return None
    return flag, tags["XS"]


def count_sam(lines: Iterable[str]) -> Counter:
    counts = Counter()
    for line_number, raw in enumerate(lines, 1):
        if not raw.strip() or raw.startswith("@"):
            continue
        fields = raw.rstrip("\r\n").split("\t")
        if len(fields) < 11:
            raise ValueError(f"SAM line {line_number}: expected at least 11 fields")
        found = eligible_strand(fields)
        if found:
            add_record(counts, *found)
    return counts


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sam", type=Path, required=True, help="SAM file or - for stdin")
    parser.add_argument("--expected", choices=sorted(EXPECTED_READ1), required=True)
    parser.add_argument("--minimum", type=float, default=MINIMUM)
    parser.add_argument("--sample-id", default=None)
    args = parser.parse_args(argv)
    try:
        if args.sam == Path("-"):
            counts = count_sam(sys.stdin)
        else:
            with args.sam.open(encoding="utf-8") as handle:
                counts = count_sam(handle)
        summary = summarize(counts, args.expected, args.minimum)
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    if args.sample_id:
        summary["sample_id"] = args.sample_id
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["status"] == "pass" else 3


if __name__ == "__main__":
    sys.exit(main())
