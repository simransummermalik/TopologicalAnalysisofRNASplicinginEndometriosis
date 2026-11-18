#!/usr/bin/env python3
"""Infer paired-library orientation from read k-mers, without an aligner.

Each read is compared with the sense and antisense strands of a few reference
mRNAs (25-mers sampled every 8 bases). In an RF (dUTP-style) library, read 1
matches the antisense strand and read 2 the sense strand; FR is the opposite;
an unstranded library gives about half of each. K-mers found on both strands
are discarded. This is a quick check of the library, independent of the
alignment settings; it does not replace per-sample QC after alignment.

Usage:
  strand_kmer_check.py --mrna reference_mrnas.fa \
    --sample SRR1 reads_1.fastq reads_2.fastq [--sample SRR2 ...]
Prints one tab-separated row per sample and mate.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

K = 25
STEP = 8
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name, parts = None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith(">"):
            if name:
                sequences[name] = "".join(parts)
            name, parts = line[1:].split()[0], []
        elif line:
            parts.append(line.upper())
    if name:
        sequences[name] = "".join(parts)
    if not sequences:
        raise ValueError(f"{path}: no FASTA records")
    return sequences


def strand_kmers(sequences) -> tuple[set[str], set[str]]:
    sense, antisense = set(), set()
    for sequence in sequences:
        for strand, target in ((sequence, sense), (reverse_complement(sequence), antisense)):
            target.update(strand[i : i + K] for i in range(len(strand) - K + 1))
    shared = sense & antisense
    return sense - shared, antisense - shared


def classify_reads(lines, sense: set[str], antisense: set[str]) -> dict[str, int]:
    counts = {"reads": 0, "sense": 0, "antisense": 0, "both": 0}
    for index, line in enumerate(lines):
        if index % 4 != 1:
            continue
        read = line.strip().upper()
        counts["reads"] += 1
        probes = [read[i : i + K] for i in range(0, len(read) - K + 1, STEP)]
        in_sense = any(probe in sense for probe in probes)
        in_antisense = any(probe in antisense for probe in probes)
        if in_sense and in_antisense:
            counts["both"] += 1
        elif in_sense:
            counts["sense"] += 1
        elif in_antisense:
            counts["antisense"] += 1
    return counts


def antisense_fraction(counts: dict[str, int]) -> float | None:
    matched = counts["sense"] + counts["antisense"]
    return counts["antisense"] / matched if matched else None


def call_library(read1: dict[str, int], read2: dict[str, int], minimum: float = 0.9) -> str:
    """RF, FR or unresolved from both mates' antisense fractions."""
    first, second = antisense_fraction(read1), antisense_fraction(read2)
    if first is None or second is None:
        return "unresolved"
    if first >= minimum and second <= 1 - minimum:
        return "RF"
    if first <= 1 - minimum and second >= minimum:
        return "FR"
    return "unresolved"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--mrna", type=Path, required=True)
    parser.add_argument(
        "--sample", nargs=3, action="append", required=True,
        metavar=("SAMPLE_ID", "READ1", "READ2"),
    )
    args = parser.parse_args(argv)
    sense, antisense = strand_kmers(read_fasta(args.mrna).values())
    print("sample_id\tmate\treads\tsense\tantisense\tboth\tantisense_fraction\tlibrary_call")
    for sample_id, read1, read2 in args.sample:
        mates = []
        for path in (read1, read2):
            with open(path, encoding="utf-8") as handle:
                mates.append(classify_reads(handle, sense, antisense))
        call = call_library(*mates)
        for mate, counts in zip(("read1", "read2"), mates):
            fraction = antisense_fraction(counts)
            print(
                f"{sample_id}\t{mate}\t{counts['reads']}\t{counts['sense']}\t"
                f"{counts['antisense']}\t{counts['both']}\t"
                f"{'' if fraction is None else f'{fraction:.4f}'}\t{call}"
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
