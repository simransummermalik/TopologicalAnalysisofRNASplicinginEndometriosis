#!/usr/bin/env python3
"""Summarize a junction-to-node map for a gene list (annotation only, no reads).

For each map (built by preprocessing/build_junction_gene_map.py) this reports
how many annotated junctions touching the listed genes map to exactly one
usable node pair, how many are ambiguous or same-node, and how many genes have
at least one independent graph cycle if every usable annotated junction were
observed. Parallel junctions between the same node pair count as one edge,
matching the converter.

Example:
  python3 scripts/node_rule_review.py --genes data/gene_sets/candidate_genes.tsv \
    --map exon-interval=maps/exon_interval.tsv \
    --map shared-splice-site=maps/shared.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def read_map(path: Path):
    mapping = defaultdict(set)
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["chromosome"], row["intron_start"], row["intron_end"], row["strand"])
            mapping[key].add((row["gene_id"], row["from_node"], row["to_node"]))
    return mapping


def cycle_rank(edges: set[tuple[str, str]]) -> int:
    parent: dict[str, str] = {}

    def find(node: str) -> str:
        parent.setdefault(node, node)
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for left, right in edges:
        parent[find(left)] = find(right)
    nodes = set(parent)
    components = len({find(node) for node in nodes})
    return len(edges) - len(nodes) + components


def summarize(mapping, genes: set[str]) -> dict[str, float | int]:
    junctions = usable = ambiguous = same_node = 0
    edges: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for targets in mapping.values():
        if not {target[0] for target in targets} & genes:
            continue
        junctions += 1
        if len(targets) != 1:
            ambiguous += 1
            continue
        gene, from_node, to_node = next(iter(targets))
        if from_node == to_node:
            same_node += 1
            continue
        usable += 1
        edges[gene].add((from_node, to_node))
    ranks = sorted(cycle_rank(gene_edges) for gene_edges in edges.values())
    return {
        "annotated_junctions": junctions,
        "usable_one_node_pair": usable,
        "usable_percent": round(100 * usable / junctions, 1) if junctions else 0.0,
        "ambiguous": ambiguous,
        "same_node": same_node,
        "genes_with_edges": len(edges),
        "genes_with_cycle_rank_above_0": sum(rank > 0 for rank in ranks),
        "median_cycle_rank": ranks[len(ranks) // 2] if ranks else 0,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--genes", type=Path, required=True, help="TSV with a gene_id column")
    parser.add_argument(
        "--map", action="append", required=True, metavar="LABEL=PATH",
        help="a junction map to summarize; repeat for each rule",
    )
    args = parser.parse_args()
    with args.genes.open(encoding="utf-8") as handle:
        genes = {row["gene_id"] for row in csv.DictReader(handle, delimiter="\t")}
    report = {}
    for item in args.map:
        label, _, path = item.partition("=")
        if not path:
            parser.error("--map must be LABEL=PATH")
        report[label] = summarize(read_map(Path(path)), genes)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
