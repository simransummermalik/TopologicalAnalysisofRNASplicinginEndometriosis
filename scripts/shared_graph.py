#!/usr/bin/env python3
"""Build comparable sample graphs from complete, unfiltered cohort counts.

Thresholds are explicit options, not approved study defaults. Edge selection
ignores condition labels. A missing edge becomes zero only in a covered
sample/gene, after the caller declares the cohort count export complete.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "preprocessing"))
from hisat2_pipeline import checksum, save

JUNCTION_COLUMNS = ["sample_id", "gene_id", "from", "to", "count"]


def read_table(path, required, exact=False):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        if len(set(header)) != len(header) or (
            header != required if exact else not set(required).issubset(header)
        ):
            raise ValueError(f"{path}: expected unique columns {required}")
        rows = []
        for line, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path} line {line}: wrong number of fields")
            if any(
                not row[key].strip()
                or row[key] != row[key].strip()
                or any(char in row[key] for char in "\t\r\n")
                for key in required
            ):
                raise ValueError(f"{path} line {line}: empty or padded required value")
            rows.append(row)
    return header, rows


def write_table(path, header, rows):
    with Path(path).open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def positive_integer(value):
    try:
        number = int(value)
    except (ValueError, TypeError):
        raise argparse.ArgumentTypeError("must be a positive integer") from None
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def cycle_rank(edges):
    """Dimension m - n + c, ignoring directions only for connectivity."""
    parent = {}

    def find(node):
        parent.setdefault(node, node)
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for source, target in edges:
        a, b = find(source), find(target)
        if a != b:
            parent[a] = b
    components = len({find(node) for node in parent})
    return len(edges) - len(parent) + components


@dataclass
class Cohort:
    samples: list[str]
    counts: dict
    genes: list[str]
    family_header: list[str]
    family_rows: dict


def load_cohort(junctions, manifest, genes=None, complete_counts=False):
    if not complete_counts:
        raise ValueError(
            "--complete-counts is required: verify successful exports for every sample, "
            "with minimum junction support 1, before treating absent edges as measured zeros"
        )
    _, samples = read_table(manifest, ["sample_id", "condition"])
    ids = [row["sample_id"] for row in samples]
    if not ids or len(set(ids)) != len(ids):
        raise ValueError(
            "manifest must contain unique sample IDs and at least one sample"
        )
    conditions = {row["condition"].lower() for row in samples}
    if conditions != {"case", "control"}:
        raise ValueError("manifest must contain both case and control samples only")
    _, rows = read_table(junctions, JUNCTION_COLUMNS, exact=True)
    counts = defaultdict(dict)
    seen = set()
    observed_samples = set()
    sample_set = set(ids)
    for line, row in enumerate(rows, 2):
        sample, gene = row["sample_id"], row["gene_id"]
        edge = row["from"], row["to"]
        if sample not in sample_set:
            raise ValueError(
                f"junction line {line}: sample {sample!r} is absent from manifest"
            )
        if edge[0] == edge[1]:
            raise ValueError(f"junction line {line}: self-loop is not supported")
        key = sample, gene, edge
        if key in seen:
            raise ValueError(f"junction line {line}: duplicate sample/gene edge")
        seen.add(key)
        # Raw read-support counts must be integers, not normalized proportions.
        text = row["count"]
        if not text.isascii() or not text.isdecimal():
            raise ValueError(
                f"junction line {line}: count must be a nonnegative integer"
            )
        value = int(text)
        # Rust uses f64 counts: avoid silently rounding large input integers.
        if value > 2**53:
            raise ValueError(
                f"junction line {line}: count exceeds exact f64 integer range"
            )
        if value:
            observed_samples.add(sample)
        counts[gene].setdefault(edge, {})[sample] = value
    absent = sorted(sample_set - observed_samples)
    if absent:
        raise ValueError(
            f"manifest samples have no positive input counts: {', '.join(absent)}; "
            "check missing/failed sample exports"
        )
    if genes:
        family_header, family = read_table(genes, ["gene_id"])
        family_rows = {row["gene_id"]: row for row in family}
        if not family or len(family_rows) != len(family):
            raise ValueError(
                "gene family must contain unique gene IDs and at least one gene"
            )
    else:
        family_header = ["gene_id"]
        family_rows = {gene: {"gene_id": gene} for gene in counts}
    return Cohort(
        sorted(ids), dict(counts), sorted(family_rows), family_header, family_rows
    )


def build_shared(cohort, minimum_edge, minimum_samples, minimum_coverage):
    for value in (minimum_edge, minimum_samples, minimum_coverage):
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise ValueError("shared graph thresholds must be positive integers")
    if minimum_samples > len(cohort.samples):
        raise ValueError("minimum samples per edge exceeds the manifest cohort size")
    junctions, edge_audit, sample_audit, gene_audit, eligible = [], [], [], [], []
    for gene in cohort.genes:
        counts = cohort.counts.get(gene, {})
        edges = []
        for (source, target), values in sorted(counts.items()):
            supported = sum(value >= minimum_edge for value in values.values())
            selected = supported >= minimum_samples
            edge_audit.append(
                dict(
                    gene_id=gene,
                    **{"from": source, "to": target},
                    samples_meeting_support=supported,
                    selected=str(selected).lower(),
                )
            )
            if selected:
                edges.append((source, target))
        rank = cycle_rank(edges)
        covered = 0
        for sample in cohort.samples:
            total = sum(counts[edge].get(sample, 0) for edge in edges)
            if total > 2**53:
                raise ValueError(
                    f"{sample}/{gene}: retained support exceeds exact f64 integer range"
                )
            coverage_ok = total >= minimum_coverage
            covered += coverage_ok
            if not edges:
                status = "no_shared_edges"
            elif rank < 1:
                status = "no_cycle_in_shared_graph"
            elif not coverage_ok:
                status = "below_coverage"
            else:
                status = "included"
            zeros = sum(counts[edge].get(sample, 0) == 0 for edge in edges)
            sample_audit.append(
                dict(
                    gene_id=gene,
                    sample_id=sample,
                    shared_graph_support=total,
                    coverage_pass=str(coverage_ok).lower(),
                    cycle_rank=rank,
                    zero_edges=zeros,
                    status=status,
                )
            )
            if status == "included":
                for source, target in edges:
                    junctions.append(
                        dict(
                            sample_id=sample,
                            gene_id=gene,
                            **{"from": source, "to": target},
                            count=counts[(source, target)].get(sample, 0),
                        )
                    )
        if not edges:
            status = "no_shared_edges"
        elif rank < 1:
            status = "no_cycle_in_shared_graph"
        elif covered != len(cohort.samples):
            status = "incomplete_coverage"
        else:
            status = "eligible"
            eligible.append(cohort.family_rows[gene])
        gene_audit.append(
            dict(
                gene_id=gene,
                shared_edges=len(edges),
                cycle_rank=rank,
                covered_samples=covered,
                required_samples=len(cohort.samples),
                status=status,
            )
        )
    return dict(
        junctions=junctions,
        edges=edge_audit,
        samples=sample_audit,
        genes=gene_audit,
        eligible=eligible,
    )


def export_shared(output, cohort, result, manifest):
    """Only covered cyclic sample/gene graphs enter the Rust input table."""
    output.mkdir()
    write_table(output / "junctions.tsv", JUNCTION_COLUMNS, result["junctions"])
    write_table(
        output / "edges.tsv",
        ["gene_id", "from", "to", "samples_meeting_support", "selected"],
        result["edges"],
    )
    write_table(
        output / "sample_eligibility.tsv",
        [
            "gene_id",
            "sample_id",
            "shared_graph_support",
            "coverage_pass",
            "cycle_rank",
            "zero_edges",
            "status",
        ],
        result["samples"],
    )
    write_table(
        output / "gene_eligibility.tsv",
        [
            "gene_id",
            "shared_edges",
            "cycle_rank",
            "covered_samples",
            "required_samples",
            "status",
        ],
        result["genes"],
    )
    write_table(output / "eligible_genes.tsv", cohort.family_header, result["eligible"])
    shutil.copyfile(manifest, output / "samples.tsv")


def prepare(args):
    output = args.output_dir.resolve()
    if output.exists():
        raise ValueError(f"output folder already exists: {output}")
    cohort = load_cohort(
        args.junctions, args.manifest, args.genes, args.complete_counts
    )
    result = build_shared(
        cohort,
        args.minimum_edge_count,
        args.minimum_samples_per_edge,
        args.minimum_gene_support,
    )
    inputs = [
        checksum(path) for path in (args.junctions, args.manifest, args.genes) if path
    ]
    record = dict(
        status="running",
        mode="shared graph preparation",
        inputs=inputs,
        implementation=checksum(Path(__file__)),
        complete_counts_declared=True,
        approval="explicit caller settings; not team approval",
        settings=dict(
            minimum_edge_count=args.minimum_edge_count,
            minimum_samples_per_edge=args.minimum_samples_per_edge,
            minimum_gene_support=args.minimum_gene_support,
            minimum_cycle_rank=1,
            eligibility="all manifest samples covered",
        ),
        cohort=cohort.samples,
        outputs=[],
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        export_shared(output, cohort, result, args.manifest)
        record.update(
            status="complete",
            eligible_genes=len(result["eligible"]),
            retained_rows=len(result["junctions"]),
        )
        record["outputs"] = [checksum(path) for path in sorted(output.iterdir())]
    except BaseException as error:
        record.update(
            status="interrupted" if isinstance(error, KeyboardInterrupt) else "failed",
            error=str(error),
        )
        raise
    finally:
        if output.is_dir():
            save(output / "run_manifest.json", record)
    print(
        f"Shared graphs: {output}; {len(result['eligible'])} genes eligible in every sample"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("junctions", "manifest", "output-dir"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument(
        "--genes", type=Path, help="optional full gene family, including absent genes"
    )
    parser.add_argument(
        "--complete-counts",
        action="store_true",
        help="declare every sample export complete, with support cutoff 1",
    )
    for name in (
        "minimum-edge-count",
        "minimum-samples-per-edge",
        "minimum-gene-support",
    ):
        parser.add_argument(f"--{name}", type=positive_integer, required=True)
    try:
        prepare(parser.parse_args())
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
