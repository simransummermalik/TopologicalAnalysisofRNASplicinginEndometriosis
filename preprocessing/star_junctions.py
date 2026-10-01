#!/usr/bin/env python3
"""Convert STAR SJ.out.tab rows into the Splice Girl input contract."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

JunctionKey = tuple[str, int, int, str]


def load_gene_map(path: Path) -> dict[JunctionKey, set[str]]:
    mapping: dict[JunctionKey, set[str]] = defaultdict(set)
    required = {"chromosome", "intron_start", "intron_end", "strand", "gene_id"}

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != required:
            raise ValueError(
                f"{path}: expected tab-separated header "
                "chromosome, intron_start, intron_end, strand, gene_id"
            )
        for line_number, row in enumerate(reader, start=2):
            try:
                chromosome = row["chromosome"].strip()
                start = int(row["intron_start"])
                end = int(row["intron_end"])
                strand = row["strand"].strip()
                gene_id = row["gene_id"].strip()
            except (TypeError, ValueError) as error:
                raise ValueError(f"{path}:{line_number}: invalid map row") from error
            if not chromosome or not gene_id or strand not in {"1", "2"}:
                raise ValueError(f"{path}:{line_number}: invalid map values")
            if start < 1 or end < start:
                raise ValueError(f"{path}:{line_number}: invalid coordinates")
            mapping[(chromosome, start, end, strand)].add(gene_id)

    return mapping


def directed_nodes(key: JunctionKey) -> tuple[str, str]:
    chromosome, intron_start, intron_end, strand = key
    if strand == "1":
        return (
            f"{chromosome}:{intron_start - 1}:+",
            f"{chromosome}:{intron_end + 1}:+",
        )
    if strand == "2":
        return (
            f"{chromosome}:{intron_end + 1}:-",
            f"{chromosome}:{intron_start - 1}:-",
        )
    raise ValueError(f"unsupported STAR strand code: {strand}")


def convert_star_junctions(
    sj_path: Path,
    gene_map_path: Path,
    output_path: Path,
    sample_id: str,
    minimum_unique_reads: int = 1,
    include_multimapping: bool = False,
    qc_path: Path | None = None,
) -> dict[str, int | str | bool]:
    gene_map = load_gene_map(gene_map_path)
    counts: dict[tuple[str, str, str], int] = defaultdict(int)
    stats: dict[str, int | str | bool] = {
        "sample_id": sample_id,
        "star_rows": 0,
        "written_rows": 0,
        "below_minimum_unique_reads": 0,
        "unannotated_junctions": 0,
        "ambiguous_gene_junctions": 0,
        "include_multimapping": include_multimapping,
        "minimum_unique_reads": minimum_unique_reads,
    }

    with sj_path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{sj_path}:{line_number}: STAR SJ.out.tab requires 9 columns"
                )
            chromosome, start_text, end_text, strand = fields[:4]
            try:
                key = (chromosome, int(start_text), int(end_text), strand)
                unique_reads = int(fields[6])
                multimapping_reads = int(fields[7])
            except ValueError as error:
                raise ValueError(
                    f"{sj_path}:{line_number}: invalid numeric STAR field"
                ) from error

            stats["star_rows"] = int(stats["star_rows"]) + 1
            if unique_reads < minimum_unique_reads:
                stats["below_minimum_unique_reads"] = (
                    int(stats["below_minimum_unique_reads"]) + 1
                )
                continue

            genes = gene_map.get(key, set())
            if not genes:
                stats["unannotated_junctions"] = (
                    int(stats["unannotated_junctions"]) + 1
                )
                continue
            if len(genes) != 1:
                stats["ambiguous_gene_junctions"] = (
                    int(stats["ambiguous_gene_junctions"]) + 1
                )
                continue

            gene_id = next(iter(genes))
            from_node, to_node = directed_nodes(key)
            count = unique_reads
            if include_multimapping:
                count += multimapping_reads
            counts[(gene_id, from_node, to_node)] += count

    if not counts:
        raise ValueError(
            "no annotated junctions passed the filters; inspect the GTF build, "
            "STAR reference, strand values, and minimum read threshold"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "gene_id", "from", "to", "count"])
        for (gene_id, from_node, to_node), count in sorted(counts.items()):
            writer.writerow([sample_id, gene_id, from_node, to_node, count])

    stats["written_rows"] = len(counts)
    if qc_path is not None:
        qc_path.parent.mkdir(parents=True, exist_ok=True)
        with qc_path.open("w", encoding="utf-8") as handle:
            json.dump(stats, handle, indent=2, sort_keys=True)
            handle.write("\n")
    return stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert STAR SJ.out.tab to a Splice Girl junction TSV."
    )
    parser.add_argument("--sj", required=True, type=Path)
    parser.add_argument("--gene-map", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--qc-json", type=Path)
    parser.add_argument("--min-unique", type=int, default=1)
    parser.add_argument("--include-multimapping", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.min_unique < 0:
        raise SystemExit("--min-unique must be zero or greater")
    for label, path in (("STAR junction", args.sj), ("gene map", args.gene_map)):
        if not path.is_file():
            raise SystemExit(f"{label} file does not exist: {path}")
    try:
        stats = convert_star_junctions(
            sj_path=args.sj,
            gene_map_path=args.gene_map,
            output_path=args.output,
            sample_id=args.sample_id,
            minimum_unique_reads=args.min_unique,
            include_multimapping=args.include_multimapping,
            qc_path=args.qc_json,
        )
    except (OSError, ValueError) as error:
        raise SystemExit(str(error)) from error
    print(
        f"Wrote {stats['written_rows']} annotated junction rows to {args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
