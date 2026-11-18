#!/usr/bin/env python3
"""Explain junction exclusions without changing the frozen assignment policy.

Full detail stays in the requested new run folder. The summary distinguishes
alternative exon pairs within one gene from cross-gene ambiguity. No ambiguous
observation is reassigned, duplicated into graphs, or used for inference.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from hisat2_pipeline import checksum, save
from sam_junctions import skipped_intervals
from star_junctions import load_gene_map

CLASSES = (
    "assigned",
    "within_gene_exon_pairs",
    "cross_gene",
    "cross_gene_and_exon_pairs",
    "same_node",
    "unannotated",
)
AMBIGUOUS = ("within_gene_exon_pairs", "cross_gene", "cross_gene_and_exon_pairs")


def classify(targets):
    if not targets:
        return "unannotated"
    if len(targets) == 1:
        _, from_node, to_node = next(iter(targets))
        return "same_node" if from_node == to_node else "assigned"
    by_gene = Counter(target[0] for target in targets)
    if len(by_gene) == 1:
        return "within_gene_exon_pairs"
    return "cross_gene_and_exon_pairs" if max(by_gene.values()) > 1 else "cross_gene"


def observe(lines, mapping, minimum_anchor, exclude_duplicates=False):
    observations = Counter()
    qc = Counter()
    strand = Counter()
    for number, raw in enumerate(lines, 1):
        if not raw.strip() or raw.startswith("@"):
            continue
        fields = raw.rstrip("\r\n").split("\t")
        if len(fields) < 11:
            raise ValueError(f"SAM line {number}: expected 11 or more fields")
        qc["sam_records"] += 1
        flag = int(fields[1])
        if flag < 0:
            raise ValueError(f"SAM line {number}: negative flag")
        if flag & (4 | 256 | 512 | 2048):
            qc["unmapped_secondary_qcfail_or_supplementary"] += 1
            continue
        if exclude_duplicates and flag & 1024:
            qc["duplicate_records_excluded"] += 1
            continue
        if "N" not in fields[5]:
            qc["unspliced_records"] += 1
            continue
        tags = {}
        for tag in fields[11:]:
            pieces = tag.split(":", 2)
            if len(pieces) == 3 and pieces[0] in {"NH", "XS"}:
                if pieces[0] in tags:
                    raise ValueError(f"SAM line {number}: duplicate {pieces[0]} tag")
                tags[pieces[0]] = tuple(pieces[1:])
        if tags.get("NH") != ("i", "1"):
            qc["nonunique_or_missing_NH_records"] += 1
            continue
        transcript = {("A", "+"): "1", ("A", "-"): "2"}.get(tags.get("XS"))
        if transcript is None:
            qc["unknown_transcript_strand_records"] += 1
            continue
        intervals = skipped_intervals(int(fields[3]), fields[5], minimum_anchor)
        if not intervals:
            qc["insufficient_anchor_records"] += 1
            continue
        mate = "read1" if flag & 64 else "read2" if flag & 128 else "unpaired"
        sense = (transcript == "1") == (not bool(flag & 16))
        strand[f"{mate}_{'sense' if sense else 'antisense'}"] += 1
        for start, end in intervals:
            key = (fields[2], start, end, transcript)
            observations[key] += 1
            kind = classify(mapping.get(key, set()))
            if kind == "assigned":
                qc["assigned_junction_observations"] += 1
            elif kind == "unannotated":
                qc["unannotated_junction_observations"] += 1
            elif kind == "same_node":
                qc["same_node_junction_observations"] += 1
            else:
                qc["ambiguous_junction_observations"] += 1
    return observations, qc, strand


def rows(path):
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def table(path, header, values):
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(values)


def run(args):
    if args.min_anchor < 1:
        raise ValueError("minimum anchor must be positive")
    if args.output_dir.exists():
        raise ValueError("output folder exists; choose a new folder")
    expected = json.loads(args.expected_qc.read_text()) if args.expected_qc else None
    if expected and (
        expected["sample_id"] != args.sample_id
        or expected["minimum_anchor"] != args.min_anchor
        or expected["exclude_duplicates"] != args.exclude_duplicates
    ):
        raise ValueError("audit settings disagree with the source QC")
    mapping = load_gene_map(args.gene_map)
    candidates = {row["gene_id"]: row for row in rows(args.candidates)}
    output = args.output_dir.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir()
    record = {
        "status": "running",
        "sample_id": args.sample_id,
        "minimum_anchor": args.min_anchor,
        "exclude_duplicates": args.exclude_duplicates,
        "policy": "diagnostic only; exact frozen exon-pair assignment remains unchanged",
        "inputs": [],
        "outputs": [],
    }
    save(output / "audit_manifest.json", record)
    try:
        source = args.bam or args.sam
        paths = [source, args.gene_map, args.candidates, Path(__file__)] + (
            [args.expected_qc] if args.expected_qc else []
        )
        record["inputs"] = [checksum(path) for path in paths]
        if args.bam:
            command = [args.samtools, "view", "-h", str(args.bam.resolve())]
            record["command"] = command
            with (output / "samtools.log").open("w") as log:
                process = subprocess.Popen(
                    command, stdout=subprocess.PIPE, stderr=log, text=True
                )
                try:
                    observed, qc, strand = observe(
                        process.stdout,
                        mapping,
                        args.min_anchor,
                        args.exclude_duplicates,
                    )
                except BaseException:
                    process.terminate()
                    process.wait()
                    raise
                finally:
                    process.stdout.close()
                code = process.wait()
                if code:
                    raise ValueError(f"samtools failed ({code}); see samtools.log")
        else:
            with args.sam.open(encoding="utf-8") as handle:
                observed, qc, strand = observe(
                    handle, mapping, args.min_anchor, args.exclude_duplicates
                )
        if expected:
            keys = {
                key
                for key in expected
                if key.endswith("_records") or key.endswith("_observations")
            } | {"sam_records", "unmapped_secondary_qcfail_or_supplementary"}
            for key in keys:
                if qc.get(key, 0) != expected.get(key, 0):
                    raise ValueError(
                        f"source QC mismatch for {key}: audit {qc.get(key, 0)} versus source {expected.get(key, 0)}"
                    )
        totals = Counter({kind: 0 for kind in CLASSES})
        key_totals = Counter({kind: 0 for kind in CLASSES})
        genes = defaultdict(Counter)
        details = []
        for key, count in sorted(observed.items()):
            targets = mapping.get(key, set())
            kind = classify(targets)
            ids = sorted({target[0] for target in targets})
            totals[kind] += count
            key_totals[kind] += 1
            for gene in ids:
                genes[gene][kind] += count
                genes[gene][kind + "_keys"] += 1
            details.append(
                [
                    *key,
                    kind,
                    count,
                    len(targets),
                    len(ids),
                    int(bool(set(ids) & candidates.keys())),
                    "|".join(ids),
                    "|".join(";".join(target) for target in sorted(targets)),
                ]
            )
        detail_header = [
            "chromosome",
            "intron_start",
            "intron_end",
            "strand",
            "assignment_class",
            "junction_observations",
            "exon_pair_targets",
            "gene_targets",
            "touches_candidate",
            "gene_ids",
            "possible_exon_pairs",
        ]
        table(output / "junction_assignments.tsv", detail_header, details)
        top = sorted(
            (row for row in details if row[4] in AMBIGUOUS),
            key=lambda row: (-row[5], tuple(row[:4])),
        )[:50]
        table(output / "top_ambiguous.tsv", detail_header, top)
        total = sum(totals.values())
        table(
            output / "assignment_summary.tsv",
            [
                "assignment_class",
                "observed_junction_keys",
                "junction_observations",
                "fraction_of_assignment_stage",
            ],
            [
                [
                    kind,
                    key_totals[kind],
                    totals[kind],
                    totals[kind] / total if total else "",
                ]
                for kind in CLASSES
            ],
        )
        candidate_rows = []
        for gene, annotation in sorted(candidates.items()):
            counts = genes[gene]
            any_observations = sum(counts[kind] for kind in CLASSES)
            candidate_rows.append(
                [
                    gene,
                    annotation.get("gene_symbol", ""),
                    annotation.get("category", ""),
                    "observed" if any_observations else "no_eligible_observation",
                    counts["assigned"],
                    *[counts[kind] for kind in AMBIGUOUS],
                    counts["assigned_keys"],
                    sum(counts[kind + "_keys"] for kind in AMBIGUOUS),
                    counts["same_node"],
                ]
            )
        table(
            output / "candidate_gene_coverage.tsv",
            [
                "gene_id",
                "gene_symbol",
                "category",
                "observation_status",
                "assigned_observations",
                "within_gene_ambiguous_potential_observations",
                "cross_gene_ambiguous_potential_observations",
                "both_ambiguous_potential_observations",
                "assigned_junction_keys",
                "ambiguous_junction_keys",
                "same_node_observations",
            ],
            candidate_rows,
        )
        summary = {
            "sample_id": args.sample_id,
            "scope": "input alignment only; no disease inference",
            "total_assignment_observations": total,
            "observations_by_class": dict(totals),
            "junction_keys_by_class": dict(key_totals),
            "filter_counts": dict(qc),
            "spliced_alignment_strand_counts": dict(strand),
            "source_qc_agreement": "PASS" if expected else "not_requested",
            "same_gene_assignment_upper_bound_if_exon_pair_ambiguity_were_ignored": totals[
                "assigned"
            ]
            + totals["within_gene_exon_pairs"],
            "limits": [
                "Potential per-gene ambiguous support is non-additive across genes.",
                "No ambiguous observations were rescued or assigned to an exon pair.",
                "Strand counts describe primary unique spliced alignments; verify the library protocol before choosing strandedness.",
            ],
        }
        save(output / "summary.json", summary)
        body = (
            "".join(
                f"<tr><td>{html.escape(kind.replace('_', ' '))}</td><td>{key_totals[kind]:,}</td><td>{totals[kind]:,}</td><td>{totals[kind] / total:.2%}</td></tr>"
                for kind in CLASSES
            )
            if total
            else ""
        )
        report = f'<!doctype html><html lang="en"><meta charset="utf-8"><title>Splice Girl junction-assignment audit</title><style>body{{max-width:65rem;margin:3rem auto;padding:0 1rem;font:17px/1.6 Georgia,serif;color:#111}}table{{border-collapse:collapse}}td,th{{padding:.4rem 1rem;border-bottom:1px solid #ccc;text-align:left}}a{{color:inherit}}</style><h1>Junction-assignment audit</h1><p>Sample: {html.escape(args.sample_id)}. This diagnoses exclusions from the input alignment; it does not change the production assignment policy.</p><table><tr><th>Assignment class</th><th>Junction keys</th><th>Junction observations</th><th>Share</th></tr>{body}</table><p>Within-gene ambiguity means the gene is identifiable but the full annotated exon pair is not. Cross-gene ambiguity means more than one gene is possible. Both means both conditions occur for a junction.</p><p>Counts are read-junction observations, not subjects or independent molecules. One read may support several junctions. Potential per-gene ambiguous support cannot be summed across genes.</p><p>Source QC agreement: {summary["source_qc_agreement"]}.</p><ul><li><a href="assignment_summary.tsv">Assignment totals</a></li><li><a href="top_ambiguous.tsv">Fifty most-supported ambiguous junctions</a></li><li><a href="candidate_gene_coverage.tsv">Candidate coverage</a></li><li><a href="summary.json">Filtering, strand counts and limits</a></li><li><a href="audit_manifest.json">Commands and checksums</a></li></ul></html>'
        (output / "report.html").write_text(report)
        record["status"] = "complete"
        record["outputs"] = [
            checksum(path)
            for path in sorted(output.iterdir())
            if path.name != "audit_manifest.json"
        ]
    except Exception as error:
        record.update(status="failed", error=str(error))
        raise
    finally:
        save(output / "audit_manifest.json", record)
    print(f"Assignment audit: {output / 'report.html'}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--sam", type=Path)
    source.add_argument("--bam", type=Path)
    for name in ("gene-map", "candidates", "output-dir"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--min-anchor", type=int, required=True)
    parser.add_argument("--expected-qc", type=Path)
    parser.add_argument("--exclude-duplicates", action="store_true")
    parser.add_argument("--samtools", default="samtools")
    try:
        run(parser.parse_args())
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
