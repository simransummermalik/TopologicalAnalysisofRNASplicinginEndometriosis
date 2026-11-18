#!/usr/bin/env python3
"""Select and independently verify one real gene under an explicit pilot plan.

This creates a one-sample technical report, never inferential results.
NumPy is required only for the independent SVD-based numerical check.
"""

from __future__ import annotations

import argparse
import csv
import html
import itertools
import json
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "preprocessing"))
from hisat2_pipeline import checksum, save


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, header, rows):
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def select_gene(
    rows, candidates, sample_id, support, minimum_edges, maximum_nodes, minimum_total
):
    groups = defaultdict(list)
    for row in rows:
        if row["sample_id"] != sample_id:
            raise ValueError("pilot input contains another sample; never pool samples")
        count = int(row["count"])
        if count < 1:
            raise ValueError("pilot junction counts must be positive integers")
        if row["gene_id"] in candidates and count >= support:
            groups[row["gene_id"]].append(row)
    eligible = []
    for gene, values in groups.items():
        nodes = {row[name] for row in values for name in ("from", "to")}
        total = sum(int(row["count"]) for row in values)
        if (
            len(values) >= minimum_edges
            and len(nodes) <= maximum_nodes
            and total >= minimum_total
        ):
            eligible.append((total, gene, values))
    if not eligible:
        raise ValueError(
            "no candidate gene meets the predeclared technical-pilot coverage/size rule"
        )
    eligible.sort(key=lambda item: (-item[0], item[1]))
    return eligible[0], len(eligible)


def gene_targets(path, gene):
    """Use the sorted annotation map, preserving ambiguity across other genes."""
    targets = {}
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        def key(row):
            return (
                row["chromosome"],
                int(row["intron_start"]),
                int(row["intron_end"]),
                row["strand"],
            )

        previous = None
        for junction, rows in itertools.groupby(reader, key):
            if previous is not None and junction <= previous:
                raise ValueError(
                    "annotation map must be in its original sorted key order"
                )
            previous = junction
            mappings = {
                (row["gene_id"], row["from_node"], row["to_node"]) for row in rows
            }
            if len(mappings) == 1:
                target = next(iter(mappings))
                if target[0] == gene:
                    targets[junction] = target[1:]
    return targets


def independent_counts(lines, targets, minimum_anchor):
    """An independent 0-based reference cursor; no converter functions imported."""
    counts = Counter()
    for line in lines:
        if line.startswith("@") or not line.strip():
            continue
        fields = line.rstrip().split("\t")
        flag = int(fields[1])
        if flag & (4 | 256 | 512 | 2048):
            continue
        tags = set(fields[11:])
        if "NH:i:1" not in tags:
            continue
        strand = "1" if "XS:A:+" in tags else "2" if "XS:A:-" in tags else None
        if strand is None:
            continue
        operations = [
            (int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", fields[5])
        ]
        cursor = int(fields[3]) - 1
        for index, (length, operation) in enumerate(operations):
            if operation == "N" and 0 < index < len(operations) - 1:
                left, right = operations[index - 1], operations[index + 1]
                if (
                    left[1] in "M=X"
                    and right[1] in "M=X"
                    and min(left[0], right[0]) >= minimum_anchor
                ):
                    edge = targets.get((fields[2], cursor + 1, cursor + length, strand))
                    if edge is not None:
                        counts[edge] += 1
            if operation in "MDN=X":
                cursor += length
    return counts


def projection(rows):
    import numpy as np

    nodes = sorted({row[name] for row in rows for name in ("from", "to")})
    node_index = {name: index for index, name in enumerate(nodes)}
    edges = sorted(rows, key=lambda row: (row["from"], row["to"]))
    incidence = np.zeros((len(nodes), len(edges)))
    signal = np.array([float(row["count"]) for row in edges])
    signal /= signal.sum()
    for index, row in enumerate(edges):
        incidence[node_index[row["from"]], index] = -1
        incidence[node_index[row["to"]], index] = 1
    potentials, _, rank, _ = np.linalg.lstsq(incidence.T, signal, rcond=None)
    gradient = incidence.T @ potentials
    cycle = signal - gradient
    result = {
        "algorithm": "NumPy SVD-based lstsq of the rectangular incidence transpose; no Rust solver reuse",
        "numpy_version": np.__version__,
        "nodes": len(nodes),
        "edges": len(edges),
        "cycle_rank": len(edges) - int(rank),
        "total_support": sum(int(row["count"]) for row in edges),
        "cycle_fraction": float(cycle @ cycle / (signal @ signal)),
        "reconstruction_error": float(np.linalg.norm(signal - gradient - cycle)),
        "kernel_error": float(np.linalg.norm(incidence @ cycle)),
        "orthogonality_error": float(abs(gradient @ cycle)),
    }
    return result, {
        (row["from"], row["to"]): (
            float(signal[i]),
            float(gradient[i]),
            float(cycle[i]),
        )
        for i, row in enumerate(edges)
    }


def run(args):
    if (
        min(
            args.min_support,
            args.min_anchor,
            args.min_edges,
            args.max_nodes,
            args.min_total,
        )
        < 1
    ):
        raise ValueError("pilot thresholds must be positive")
    if args.output_dir.exists():
        raise ValueError("output folder exists; choose a new folder")
    candidates = {row["gene_id"]: row for row in read_tsv(args.candidates)}
    rows = read_tsv(args.junctions)
    (_, gene, selected), eligible = select_gene(
        rows,
        candidates,
        args.sample_id,
        args.min_support,
        args.min_edges,
        args.max_nodes,
        args.min_total,
    )
    all_gene_rows = [row for row in rows if row["gene_id"] == gene]
    source_samples = [
        row
        for row in read_tsv(args.sample_manifest)
        if row["sample_id"] == args.sample_id
    ]
    if len(source_samples) != 1:
        raise ValueError("sample must have exactly one original metadata row")
    condition = source_samples[0]["condition"]
    if condition not in {"case", "control"}:
        raise ValueError("source sample must have an explicit case/control condition")
    import numpy

    subset = json.loads(args.subset_manifest.read_text())
    alignment_manifest = json.loads(args.alignment_manifest.read_text())
    reference_path = args.gene_map.parent / "reference_manifest.json"
    reference_manifest = json.loads(reference_path.read_text())
    if (
        subset.get("status") != "complete"
        or alignment_manifest.get("status") != "complete"
        or reference_manifest.get("status") != "complete"
    ):
        raise ValueError("subset and alignment manifests must both be complete")
    settings = alignment_manifest["settings"]
    if (
        settings["sample_id"] != args.sample_id
        or settings["minimum_anchor"] != args.min_anchor
        or settings["minimum_unique_reads"] != 1
        or settings["exclude_duplicates"]
    ):
        raise ValueError(
            "alignment/count settings disagree with the pilot recount policy"
        )
    aligned_inputs = {
        (item["path"], item["sha256"]) for item in alignment_manifest["inputs"]
    }
    if not all(
        (item["path"], item["sha256"]) in aligned_inputs for item in subset["outputs"]
    ):
        raise ValueError("alignment inputs do not match the verified subset outputs")

    output = args.output_dir.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()
    record = {
        "status": "running",
        "scope": "one biological sample; deterministic read-prefix technical pilot only; no disease inference",
        "sample_id": args.sample_id,
        "gene_id": gene,
        "gene_annotation": candidates[gene],
        "eligible_candidate_genes": eligible,
        "selection": {
            "order": "total retained support descending, gene_id ascending; no cycle-fraction selection",
            "minimum_junction_support": args.min_support,
            "minimum_anchor": args.min_anchor,
            "minimum_edges": args.min_edges,
            "maximum_nodes": args.max_nodes,
            "minimum_total_support": args.min_total,
        },
        "inputs": [],
        "outputs": [],
        "commands": [],
        "numpy_version": numpy.__version__,
        "subset_pairs": subset["written_pairs"],
    }
    save(output / "verification.json", record)
    try:
        paths = [
            args.junctions,
            args.candidates,
            args.gene_map,
            args.sample_manifest,
            args.bam,
            Path(str(args.bam) + ".bai"),
            args.subset_manifest,
            args.alignment_manifest,
            args.plan,
            reference_path,
            args.binary,
            Path(__file__),
        ]
        record["inputs"] = [checksum(path) for path in paths]
        for source, name in [
            (args.subset_manifest, "subset_manifest.json"),
            (args.alignment_manifest, "alignment_manifest.json"),
            (reference_path, "reference_manifest.json"),
        ]:
            shutil.copyfile(source, output / name)
        write_tsv(
            output / "junctions.tsv",
            ["sample_id", "gene_id", "from", "to", "count"],
            [
                [row[key] for key in ("sample_id", "gene_id", "from", "to", "count")]
                for row in selected
            ],
        )
        write_tsv(
            output / "samples.tsv",
            ["sample_id", "condition"],
            [[args.sample_id, condition]],
        )
        shutil.copyfile(args.plan, output / "plan.md")
        targets = gene_targets(args.gene_map, gene)
        positions = []
        chromosomes = set()
        for pair in targets.values():
            for node in pair:
                chromosome, interval, _ = node.split(":")
                chromosomes.add(chromosome)
                positions.extend(map(int, interval.split("-")))
        if len(chromosomes) != 1:
            raise ValueError(
                "pilot gene spans multiple reference sequences; manual review required"
            )
        region = f"{next(iter(chromosomes))}:{min(positions)}-{max(positions)}"
        command = [args.samtools, "view", "-F", "2820", str(args.bam.resolve()), region]
        record["commands"].append(command)
        alignment = subprocess.run(command, capture_output=True, text=True, check=True)
        counts = independent_counts(
            alignment.stdout.splitlines(), targets, args.min_anchor
        )
        expected = {
            (row["from"], row["to"]): int(row["count"]) for row in all_gene_rows
        }
        if dict(counts) != expected:
            raise ValueError(
                "independent BAM junction recount disagrees with converter output"
            )
        record["junction_recount"] = {
            "status": "PASS",
            "reference_region": region,
            "primary_region_records": len(alignment.stdout.splitlines()),
            "annotated_edges": len(counts),
            "supporting_read_junction_observations": sum(counts.values()),
        }
        quickcheck = [args.samtools, "quickcheck", "-v", str(args.bam.resolve())]
        record["commands"].append(quickcheck)
        subprocess.run(quickcheck, capture_output=True, text=True, check=True)
        flagstat = [args.samtools, "flagstat", str(args.bam.resolve())]
        record["commands"].append(flagstat)
        result = subprocess.run(flagstat, capture_output=True, text=True, check=True)
        (output / "alignment_flagstat.txt").write_text(result.stdout)
        for name in ["alignment_summary.txt", "junction_qc.json"]:
            shutil.copyfile(args.alignment_manifest.parent / name, output / name)
        for command in [
            [
                str(args.binary.resolve()),
                "analyze",
                str(output / "junctions.tsv"),
                "--output",
                str(output / "results"),
            ],
            [
                str(args.binary.resolve()),
                "visualize",
                args.sample_id,
                gene,
                str(output / "junctions.tsv"),
                str(output / "graph.svg"),
            ],
        ]:
            record["commands"].append(command)
            started = time.perf_counter()
            result = subprocess.run(command, capture_output=True, text=True, check=True)
            if command[1] == "analyze":
                record[
                    "rust_analysis_seconds_single_run_including_startup_and_export"
                ] = time.perf_counter() - started
            (
                output / ("analyze.log" if command[1] == "analyze" else "visualize.log")
            ).write_text(result.stdout + result.stderr)
        independent, components = projection(selected)
        rust_metrics = read_tsv(output / "results/sample_metrics.tsv")[0]
        max_delta = 0.0
        for row in read_tsv(output / "results/junction_components.tsv"):
            calculated = components[(row["from"], row["to"])]
            exported = [
                float(row[name])
                for name in (
                    "normalized_signal",
                    "gradient_component",
                    "cycle_component",
                )
            ]
            max_delta = max(
                max_delta, *(abs(a - b) for a, b in zip(calculated, exported))
            )
        fraction_delta = abs(
            float(rust_metrics["cycle_fraction"]) - independent["cycle_fraction"]
        )
        if (
            max_delta > 1e-8
            or fraction_delta > 1e-8
            or rust_metrics["checks_status"] != "PASS"
            or max(
                independent[name]
                for name in (
                    "reconstruction_error",
                    "kernel_error",
                    "orthogonality_error",
                )
            )
            > 1e-8
        ):
            raise ValueError(
                "independent SVD projection disagrees with Rust or fails math checks"
            )
        independent.update(
            status="PASS",
            maximum_exported_component_difference=max_delta,
            cycle_fraction_difference=fraction_delta,
        )
        record["independent_math"] = independent
        sensitivity = []
        for minimum in sorted({1, 3, 5, 10, args.min_support}):
            retained = [row for row in all_gene_rows if int(row["count"]) >= minimum]
            values = projection(retained)[0] if retained else {}
            if retained:
                directory = output / "sensitivity" / f"support_{minimum}"
                directory.mkdir(parents=True)
                table = directory / "junctions.tsv"
                write_tsv(
                    table,
                    ["sample_id", "gene_id", "from", "to", "count"],
                    [
                        [
                            row[key]
                            for key in ("sample_id", "gene_id", "from", "to", "count")
                        ]
                        for row in retained
                    ],
                )
                command = [
                    str(args.binary.resolve()),
                    "analyze",
                    str(table),
                    "--output",
                    str(directory / "results"),
                ]
                record["commands"].append(command)
                result = subprocess.run(
                    command, capture_output=True, text=True, check=True
                )
                (directory / "analysis.log").write_text(result.stdout + result.stderr)
                measured = read_tsv(directory / "results/sample_metrics.tsv")[0]
                if (
                    measured["checks_status"] != "PASS"
                    or abs(float(measured["cycle_fraction"]) - values["cycle_fraction"])
                    > 1e-8
                ):
                    raise ValueError(
                        f"Rust/SVD sensitivity disagreement at support {minimum}"
                    )
            sensitivity.append(
                [
                    minimum,
                    values.get("nodes", 0),
                    values.get("edges", 0),
                    values.get("total_support", 0),
                    values.get("cycle_rank", ""),
                    values.get("cycle_fraction", ""),
                    "measured" if retained else "no_retained_junctions",
                ]
            )
        write_tsv(
            output / "support_sensitivity.tsv",
            [
                "minimum_support",
                "nodes",
                "edges",
                "total_support",
                "cycle_rank",
                "cycle_fraction",
                "status",
            ],
            sensitivity,
        )
        symbol = html.escape(candidates[gene].get("gene_symbol", gene))
        escape = html.escape
        from string import Template

        template = Path(__file__).with_name("real_pilot_report.html")
        record["inputs"].append(checksum(template))
        report = Template(template.read_text()).substitute(
            symbol=symbol,
            gene=escape(gene),
            sample=escape(args.sample_id),
            condition=escape(condition),
            pairs=f"{subset['written_pairs']:,}",
            category=escape(
                candidates[gene].get("category", "not supplied").replace("_", " ")
            ),
            nodes=independent["nodes"],
            edges=independent["edges"],
            total_support=independent["total_support"],
            cycle_rank=independent["cycle_rank"],
            cycle_fraction=f"{independent['cycle_fraction']:.10f}",
        )
        (output / "report.html").write_text(report)
        record["status"] = "complete"
        record["outputs"] = [
            checksum(path)
            for path in sorted(output.rglob("*"))
            if path.is_file() and path.name != "verification.json"
        ]
    except Exception as error:
        record.update(status="failed", error=str(error))
        raise
    finally:
        save(output / "verification.json", record)
    print(f"Verified {gene}: {output / 'report.html'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    for name in (
        "junctions",
        "candidates",
        "gene-map",
        "bam",
        "sample-manifest",
        "subset-manifest",
        "alignment-manifest",
        "plan",
        "output-dir",
    ):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument(
        "--binary", type=Path, default=ROOT / "target/release/splice-girl"
    )
    parser.add_argument("--samtools", default="samtools")
    for name in ("min-support", "min-anchor", "min-edges", "max-nodes", "min-total"):
        parser.add_argument(f"--{name}", type=int, required=True)
    try:
        run(parser.parse_args())
    except (OSError, ValueError, subprocess.SubprocessError, ImportError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
