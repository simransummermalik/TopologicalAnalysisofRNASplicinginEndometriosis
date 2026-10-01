#!/usr/bin/env python3
"""Build an annotated STAR-junction-to-gene map from a GTF file.

The script uses adjacent exons within each transcript. It intentionally keeps
the output small and explicit so the later STAR conversion never has to guess
which gene owns a splice junction.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        key, separator, value = item.partition(" ")
        if not separator:
            continue
        attributes[key] = value.strip().strip('"')
    return attributes


def build_gene_map(gtf_path: Path, output_path: Path) -> dict[str, int]:
    transcripts: dict[tuple[str, str, str, str], set[tuple[int, int]]] = (
        defaultdict(set)
    )
    exon_rows = 0

    with gtf_path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{gtf_path}:{line_number}: expected 9 GTF columns, "
                    f"found {len(fields)}"
                )
            chromosome, _, feature, start_text, end_text, _, strand, _, attrs = (
                fields
            )
            if feature != "exon":
                continue
            if strand not in {"+", "-"}:
                continue

            attributes = parse_attributes(attrs)
            gene_id = attributes.get("gene_id")
            transcript_id = attributes.get("transcript_id")
            if not gene_id or not transcript_id:
                raise ValueError(
                    f"{gtf_path}:{line_number}: exon is missing gene_id or "
                    "transcript_id"
                )

            start = int(start_text)
            end = int(end_text)
            if start < 1 or end < start:
                raise ValueError(
                    f"{gtf_path}:{line_number}: invalid exon coordinates"
                )

            transcripts[(chromosome, strand, gene_id, transcript_id)].add(
                (start, end)
            )
            exon_rows += 1

    junction_genes: dict[tuple[str, int, int, str], set[str]] = defaultdict(set)
    skipped_overlaps = 0
    for (chromosome, strand, gene_id, _), exon_set in transcripts.items():
        exons = sorted(exon_set)
        for left, right in zip(exons, exons[1:]):
            intron_start = left[1] + 1
            intron_end = right[0] - 1
            if intron_start > intron_end:
                skipped_overlaps += 1
                continue
            star_strand = "1" if strand == "+" else "2"
            junction_genes[
                (chromosome, intron_start, intron_end, star_strand)
            ].add(gene_id)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["chromosome", "intron_start", "intron_end", "strand", "gene_id"]
        )
        for key in sorted(junction_genes):
            for gene_id in sorted(junction_genes[key]):
                writer.writerow([*key, gene_id])

    return {
        "exon_rows": exon_rows,
        "transcripts": len(transcripts),
        "junctions": len(junction_genes),
        "rows_written": sum(len(genes) for genes in junction_genes.values()),
        "skipped_overlapping_exons": skipped_overlaps,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a STAR junction-to-gene map from adjacent GTF exons."
    )
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.gtf.is_file():
        raise SystemExit(f"GTF file does not exist: {args.gtf}")
    try:
        stats = build_gene_map(args.gtf, args.output)
    except (OSError, ValueError) as error:
        raise SystemExit(str(error)) from error

    print(f"Wrote {stats['rows_written']} junction-to-gene rows to {args.output}")
    print(
        "Source summary: "
        f"{stats['exon_rows']} exons, {stats['transcripts']} transcripts, "
        f"{stats['skipped_overlapping_exons']} overlapping exon pairs skipped"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
