#!/usr/bin/env python3
"""Build an annotated junction-to-exon map from a GTF file.

The script uses adjacent exons within each transcript. It intentionally keeps
the output explicit so later junction conversion never has to guess which gene
or annotated exon pair owns a splice junction.

Two node rules are available (see docs/data_methods.md):

- ``shared-splice-site`` (default, decided 2026-11-18): within one gene and
  strand, annotated exons that share a start or an end coordinate form one
  node. A junction read sees only those splice sites, so variants differing
  only at their far ends are one exon. An exon that strictly contains an
  annotated intron of its own gene (a retained-intron exon) does not join its
  start to its end, so it cannot merge the exons on either side of that
  intron when both sides already hold a normal exon. Every node therefore
  contains at least one whole exon, and its ID is ``chrom:lowest-highest:strand``
  over its exons; the lowest coordinate is a start site, which belongs to one
  node only, so IDs are unique within a gene. Junctions whose two ends fall in one node are kept in
  the map (so splice-site hints are unchanged) and excluded by the converters.
- ``exon-interval`` (frozen 2026-10-22, used by the 2026-11-17 pilot): one node
  per exact exon interval; a junction matching several exon pairs is
  ambiguous.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path
from typing import TextIO

JunctionKey = tuple[str, int, int, str]
JunctionTarget = tuple[str, str, str]


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


def open_gtf(path: Path) -> TextIO:
    """Open a plain or gzip-compressed GTF as UTF-8 text."""

    if path.suffix == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def exon_node(chromosome: str, start: int, end: int, strand: str) -> str:
    """Return the frozen Splice Girl ID for one annotated exon interval."""

    return f"{chromosome}:{start}-{end}:{strand}"


NODE_RULES = ("shared-splice-site", "exon-interval")


class _UnionFind:
    def __init__(self) -> None:
        self.parent: dict = {}

    def find(self, item):
        self.parent.setdefault(item, item)
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != root:
            self.parent[item], item = root, self.parent[item]
        return root

    def union(self, left, right) -> None:
        self.parent[self.find(left)] = self.find(right)


def _contains_intron(starts: list[int], suffix_min_end: list[int], exon) -> bool:
    """True when some intron lies strictly inside the exon (start < i0, i1 < end)."""
    index = bisect_right(starts, exon[0])
    return index < len(starts) and suffix_min_end[index] < exon[1]


def shared_site_nodes(chromosome: str, strand: str, exons, introns):
    """Map each exon start ('S', pos) and end ('E', pos) key to a node ID."""
    ordered = sorted(introns)
    starts = [start for start, _ in ordered]
    suffix_min_end = [0] * len(ordered)
    running = None
    for index in range(len(ordered) - 1, -1, -1):
        end = ordered[index][1]
        running = end if running is None else min(running, end)
        suffix_min_end[index] = running
    groups = _UnionFind()
    regular = set()
    retained_exons = []
    for exon in exons:
        if _contains_intron(starts, suffix_min_end, exon):
            retained_exons.append(exon)
        else:
            groups.union(("S", exon[0]), ("E", exon[1]))
    for exon in exons:
        if exon not in retained_exons:
            regular.add(groups.find(("S", exon[0])))
    # A retained-intron exon joins its start to its end unless that would merge
    # two different nodes that each already hold a normal exon.
    retained = 0
    for exon in retained_exons:
        start_root, end_root = groups.find(("S", exon[0])), groups.find(("E", exon[1]))
        if start_root != end_root and start_root in regular and end_root in regular:
            retained += 1
            continue
        keeps_regular = start_root in regular or end_root in regular
        groups.union(("S", exon[0]), ("E", exon[1]))
        if keeps_regular:
            regular.add(groups.find(("S", exon[0])))
    members: dict = defaultdict(list)
    for key in list(groups.parent):
        members[groups.find(key)].append(key[1])
    node_of_root = {}
    seen_ids = set()
    for root, positions in members.items():
        node = exon_node(chromosome, min(positions), max(positions), strand)
        if node in seen_ids:
            raise ValueError(f"node ID collision for {node}; refine the node naming rule")
        seen_ids.add(node)
        node_of_root[root] = node
    return {key: node_of_root[groups.find(key)] for key in groups.parent}, retained, len(members)


def build_gene_map(
    gtf_path: Path, output_path: Path, node_rule: str = "shared-splice-site"
) -> dict[str, int | str]:
    if node_rule not in NODE_RULES:
        raise ValueError(f"node_rule must be one of {', '.join(NODE_RULES)}")
    transcripts: dict[tuple[str, str, str, str], set[tuple[int, int]]] = defaultdict(
        set
    )
    exon_rows = 0

    with open_gtf(gtf_path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{gtf_path}:{line_number}: expected 9 GTF columns, "
                    f"found {len(fields)}"
                )
            chromosome, _, feature, start_text, end_text, _, strand, _, attrs = fields
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
                raise ValueError(f"{gtf_path}:{line_number}: invalid exon coordinates")

            transcripts[(chromosome, strand, gene_id, transcript_id)].add((start, end))
            exon_rows += 1

    # Every annotated intron, with the exon pair that flanks it in a transcript.
    gene_exons: dict[tuple[str, str, str], set[tuple[int, int]]] = defaultdict(set)
    gene_introns: dict[tuple[str, str, str], set[tuple[int, int]]] = defaultdict(set)
    pairs: list[tuple[str, str, str, tuple[int, int], tuple[int, int]]] = []
    skipped_overlaps = 0
    for (chromosome, strand, gene_id, _), exon_set in transcripts.items():
        exons = sorted(exon_set)
        gene_exons[(chromosome, strand, gene_id)].update(exons)
        for left, right in zip(exons, exons[1:]):
            if left[1] + 1 > right[0] - 1:
                skipped_overlaps += 1
                continue
            gene_introns[(chromosome, strand, gene_id)].add((left[1] + 1, right[0] - 1))
            pairs.append((chromosome, strand, gene_id, left, right))

    node_maps = {}
    retained_exons = nodes = 0
    if node_rule == "shared-splice-site":
        for gene_key, exons in gene_exons.items():
            node_maps[gene_key], retained, count = shared_site_nodes(
                gene_key[0], gene_key[1], sorted(exons), gene_introns[gene_key]
            )
            retained_exons += retained
            nodes += count
    else:
        nodes = sum(len(exons) for exons in gene_exons.values())

    junction_targets: dict[JunctionKey, set[JunctionTarget]] = defaultdict(set)
    for chromosome, strand, gene_id, left, right in pairs:
        if node_rule == "shared-splice-site":
            lookup = node_maps[(chromosome, strand, gene_id)]
            left_node, right_node = lookup[("E", left[1])], lookup[("S", right[0])]
        else:
            left_node = exon_node(chromosome, left[0], left[1], strand)
            right_node = exon_node(chromosome, right[0], right[1], strand)
        if strand == "+":
            strand_code = "1"
            from_node, to_node = left_node, right_node
        else:
            strand_code = "2"
            from_node, to_node = right_node, left_node
        junction_targets[(chromosome, left[1] + 1, right[0] - 1, strand_code)].add(
            (gene_id, from_node, to_node)
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "chromosome",
                "intron_start",
                "intron_end",
                "strand",
                "gene_id",
                "from_node",
                "to_node",
            ]
        )
        for key in sorted(junction_targets):
            for target in sorted(junction_targets[key]):
                writer.writerow([*key, *target])

    return {
        "node_rule": node_rule,
        "exon_rows": exon_rows,
        "transcripts": len(transcripts),
        "nodes": nodes,
        "retained_intron_exons_kept_apart": retained_exons,
        "junctions": len(junction_targets),
        "rows_written": sum(len(targets) for targets in junction_targets.values()),
        "same_node_junction_rows": sum(
            1
            for targets in junction_targets.values()
            for _, from_node, to_node in targets
            if from_node == to_node
        ),
        "skipped_overlapping_exons": skipped_overlaps,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("Build a junction-to-gene-and-exon map from adjacent GTF exons.")
    )
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--node-rule",
        choices=NODE_RULES,
        default="shared-splice-site",
        help="exon node rule; exon-interval reproduces the 2026-11-17 pilot",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.gtf.is_file():
        raise SystemExit(f"GTF file does not exist: {args.gtf}")
    try:
        stats = build_gene_map(args.gtf, args.output, args.node_rule)
    except (OSError, ValueError) as error:
        raise SystemExit(str(error)) from error

    print(f"Wrote {stats['rows_written']} junction annotation rows to {args.output}")
    print(
        "Source summary: "
        f"{stats['exon_rows']} exons, {stats['transcripts']} transcripts, "
        f"{stats['skipped_overlapping_exons']} overlapping exon pairs skipped; "
        f"node rule {stats['node_rule']}: {stats['nodes']} nodes, "
        f"{stats['same_node_junction_rows']} same-node junction rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
