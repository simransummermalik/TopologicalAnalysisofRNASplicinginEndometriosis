#!/usr/bin/env python3
"""Count unique, primary HISAT2 splice-junction reads under frozen exon rules.

SAM POS and skipped intervals use 1-based closed coordinates. XS supplies the
transcript strand; the read's reverse-complement flag is not a substitute.
Counts are per read alignment, not per fragment, and retain PCR duplicates
unless --exclude-duplicates is explicitly requested.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable

from star_junctions import load_gene_map
from strand_orientation import add_record, summarize

CIGAR = re.compile(r"([1-9][0-9]*)([MIDNSHP=X])")


def skipped_intervals(
    position: int, cigar: str, min_anchor: int
) -> list[tuple[int, int]]:
    tokens = CIGAR.findall(cigar)
    if position < 1 or not tokens or "".join(n + op for n, op in tokens) != cigar:
        raise ValueError(f"invalid mapped POS/CIGAR: {position} {cigar}")
    intervals = []
    reference = position
    for index, (length, op) in enumerate(tokens):
        size = int(length)
        if op == "N":
            left = tokens[index - 1] if index else ("0", "")
            right = tokens[index + 1] if index + 1 < len(tokens) else ("0", "")
            if (
                left[1] in {"M", "=", "X"}
                and right[1] in {"M", "=", "X"}
                and min(int(left[0]), int(right[0])) >= min_anchor
            ):
                intervals.append((reference, reference + size - 1))
        if op in {"M", "D", "N", "=", "X"}:
            reference += size
    return intervals


def count_junctions(
    lines: Iterable[str],
    gene_map_path: Path,
    sample_id: str,
    min_anchor: int,
    min_unique: int,
    exclude_duplicates: bool = False,
):
    if not sample_id or any(char in sample_id for char in "\t\r\n"):
        raise ValueError("sample_id must be nonempty and contain no tab or newline")
    if min_anchor < 1 or min_unique < 1:
        raise ValueError("min-anchor and min-unique must both be positive")
    mapping = load_gene_map(gene_map_path)
    stats = Counter()
    ambiguous_target_multiplicity = Counter()
    counts = Counter()
    orientation = Counter()
    for line_number, raw in enumerate(lines, 1):
        if not raw.strip() or raw.startswith("@"):
            continue
        fields = raw.rstrip("\r\n").split("\t")
        stats["sam_records"] += 1
        if len(fields) < 11:
            raise ValueError(f"SAM line {line_number}: expected at least 11 fields")
        try:
            flag, position = int(fields[1]), int(fields[3])
        except ValueError as error:
            raise ValueError(f"SAM line {line_number}: invalid FLAG/POS") from error
        if flag < 0:
            raise ValueError(f"SAM line {line_number}: negative FLAG")
        if flag & (4 | 256 | 512 | 2048):
            stats["unmapped_secondary_qcfail_or_supplementary"] += 1
            continue
        if exclude_duplicates and flag & 1024:
            stats["duplicate_records_excluded"] += 1
            continue
        if "N" not in fields[5]:
            stats["unspliced_records"] += 1
            continue
        tags = {}
        for tag in fields[11:]:
            parts = tag.split(":", 2)
            if len(parts) == 3 and parts[0] in {"NH", "XS"}:
                if parts[0] in tags:
                    raise ValueError(
                        f"SAM line {line_number}: duplicate {parts[0]} tag"
                    )
                tags[parts[0]] = (parts[1], parts[2])
        if tags.get("NH") != ("i", "1"):
            stats["nonunique_or_missing_NH_records"] += 1
            continue
        strand = {("A", "+"): "1", ("A", "-"): "2"}.get(tags.get("XS"))
        if strand is None:
            stats["unknown_transcript_strand_records"] += 1
            continue
        add_record(orientation, flag, tags["XS"][1])
        try:
            intervals = skipped_intervals(position, fields[5], min_anchor)
        except ValueError as error:
            raise ValueError(f"SAM line {line_number}: {error}") from error
        if not intervals:
            stats["insufficient_anchor_records"] += 1
        for start, end in intervals:
            targets = mapping.get((fields[2], start, end, strand), set())
            if not targets:
                stats["unannotated_junction_observations"] += 1
            elif len(targets) != 1:
                stats["ambiguous_junction_observations"] += 1
                ambiguous_target_multiplicity[len(targets)] += 1
                if len({target[0] for target in targets}) > 1:
                    stats["ambiguous_across_genes"] += 1
                else:
                    stats["ambiguous_within_one_gene"] += 1
            elif next(iter(targets))[1] == next(iter(targets))[2]:
                # Both ends fall in one node; a self-loop has no graph meaning.
                stats["same_node_junction_observations"] += 1
            else:
                counts[next(iter(targets))] += 1
                stats["assigned_junction_observations"] += 1
    filtered = {key: count for key, count in counts.items() if count >= min_unique}
    stats["below_minimum_unique_edges"] = len(counts) - len(filtered)
    stats["written_rows"] = len(filtered)
    assignment_counts = {
        "assigned": stats.get("assigned_junction_observations", 0),
        "ambiguous": stats.get("ambiguous_junction_observations", 0),
        "unannotated": stats.get("unannotated_junction_observations", 0),
        "same_node": stats.get("same_node_junction_observations", 0),
    }
    assignment_total = sum(assignment_counts.values())
    stats["junction_assignment_summary"] = {
        "denominator": assignment_total,
        "denominator_unit": (
            "anchor-qualified junction observations after primary, unique-alignment "
            "and known-strand filters"
        ),
        **assignment_counts,
        "ambiguous_across_genes": stats.get("ambiguous_across_genes", 0),
        "ambiguous_within_one_gene": stats.get("ambiguous_within_one_gene", 0),
        "fraction_assigned": (
            assignment_counts["assigned"] / assignment_total
            if assignment_total
            else None
        ),
        "fraction_ambiguous": (
            assignment_counts["ambiguous"] / assignment_total
            if assignment_total
            else None
        ),
        "fraction_unannotated": (
            assignment_counts["unannotated"] / assignment_total
            if assignment_total
            else None
        ),
        "fraction_same_node": (
            assignment_counts["same_node"] / assignment_total
            if assignment_total
            else None
        ),
        "ambiguous_target_multiplicity": {
            str(target_count): observations
            for target_count, observations in sorted(
                ambiguous_target_multiplicity.items()
            )
        },
    }
    stats["strand_orientation"] = summarize(orientation)
    return filtered, {
        **dict(stats),
        "sample_id": sample_id,
        "minimum_anchor": min_anchor,
        "minimum_unique_reads": min_unique,
        "count_unit": "primary read alignments",
        "exclude_duplicates": exclude_duplicates,
    }


def convert_sam(
    sam_path: Path,
    gene_map: Path,
    sample_id: str,
    output: Path,
    qc: Path,
    min_anchor: int,
    min_unique: int,
    exclude_duplicates=False,
):
    with sam_path.open(encoding="utf-8") as handle:
        return convert_lines(
            handle,
            gene_map,
            sample_id,
            output,
            qc,
            min_anchor,
            min_unique,
            exclude_duplicates,
        )


def convert_lines(
    lines,
    gene_map,
    sample_id,
    output,
    qc,
    min_anchor,
    min_unique,
    exclude_duplicates=False,
):
    """Apply the same converter to a file or a supervised SAM stream."""
    if output.exists() or qc.exists():
        raise ValueError("refusing to overwrite junction TSV or QC output")
    counts, stats = count_junctions(
        lines, gene_map, sample_id, min_anchor, min_unique, exclude_duplicates
    )
    qc.parent.mkdir(parents=True, exist_ok=True)
    with qc.open("x", encoding="utf-8") as handle:
        json.dump(stats, handle, indent=2, sort_keys=True)
    if not counts:
        raise ValueError(
            f"no uniquely mapped annotated junctions passed; inspect {qc}, SAM NH/XS tags, annotation, anchors, and support thresholds"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "gene_id", "from", "to", "count"])
        writer.writerows(
            (sample_id, *key, count) for key, count in sorted(counts.items())
        )
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("sam", "gene-map", "output", "qc-json"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--min-anchor", type=int, required=True)
    parser.add_argument("--min-unique", type=int, required=True)
    parser.add_argument("--exclude-duplicates", action="store_true")
    args = parser.parse_args()
    try:
        converter = convert_lines if args.sam == Path("-") else convert_sam
        stats = converter(
            sys.stdin if args.sam == Path("-") else args.sam,
            args.gene_map,
            args.sample_id,
            args.output,
            args.qc_json,
            args.min_anchor,
            args.min_unique,
            args.exclude_duplicates,
        )
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    print(
        f"Wrote {stats['written_rows']} uniquely supported junctions to {args.output}"
    )
    summary = stats["junction_assignment_summary"]
    denominator = summary["denominator"]
    if denominator:
        print(
            "Junction assignment QC (eligible observations): "
            f"{summary['assigned']:,}/{denominator:,} uniquely assigned "
            f"({summary['fraction_assigned']:.1%}); "
            f"{summary['ambiguous']:,} ambiguous "
            f"({summary['fraction_ambiguous']:.1%}); "
            f"{summary['unannotated']:,} unannotated "
            f"({summary['fraction_unannotated']:.1%}); "
            f"{summary['same_node']:,} same-node."
        )
        if summary["ambiguous"]:
            print(
                "Ambiguous mappings: "
                f"{summary['ambiguous_across_genes']:,} across genes, "
                f"{summary['ambiguous_within_one_gene']:,} within one gene; "
                "target multiplicities are recorded in the QC JSON."
            )
    orientation = stats["strand_orientation"]
    for mate in ("read1", "read2"):
        fraction = orientation[f"{mate}_antisense_fraction"]
        if fraction is not None:
            print(f"Strand check: {mate} antisense in {fraction:.2%} of spliced unique alignments.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
