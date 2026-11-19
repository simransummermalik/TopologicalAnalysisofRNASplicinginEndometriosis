#!/usr/bin/env python3
"""Recheck all observed candidate genes in the saved one-control biological pilot.

Run from the repository root with a Python containing NumPy. This deliberately
reuses the historical exon-interval/unstranded pilot, not the final RF/shared
cohort protocol. It never downloads reads, aligns again or computes inference.
Outputs go to a new directory supplied with --output-dir.
"""

import argparse
import csv
import hashlib
import html
import itertools
import json
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from validate_real_pilot import independent_counts, projection, read_tsv, write_tsv


def identity(path):
    with path.open("rb") as source:
        digest = hashlib.file_digest(source, "sha256").hexdigest()
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": digest}


def run(args):
    import numpy as np

    output = args.output_dir.resolve()
    if output.exists():
        raise ValueError("output already exists; choose a new directory")
    previous = ROOT / "data/pilots/SRR36710387_prefix1m_nov17/verified"
    recorded = json.loads((previous / "verification.json").read_text())
    expected = {str(Path(item["path"]).resolve()): item["sha256"] for item in recorded["inputs"]}
    paths = {
        "junctions": ROOT / "preprocessing/runs/SRR36710387_pilot_nov17/alignment/junctions.tsv",
        "gene_map": ROOT / "references/hisat2_v50_pilot_nov17/junction_exon_map.tsv",
        "bam": ROOT / "preprocessing/runs/SRR36710387_pilot_nov17/alignment/aligned.bam",
        "candidates": ROOT / "data/gene_sets/candidate_genes.tsv",
    }
    inputs = {role: identity(path) for role, path in paths.items()}
    for role, item in inputs.items():
        if expected.get(item["path"]) != item["sha256"]:
            raise ValueError(f"{role} differs from the independently verified pilot")
    family = {row["gene_id"]: row for row in read_tsv(paths["candidates"])}
    raw = read_tsv(paths["junctions"])
    if {row["sample_id"] for row in raw} != {"SRR36710387"}:
        raise ValueError("expected exactly the historical control sample")
    selected = [row for row in raw if row["gene_id"] in family]
    if not selected:
        raise ValueError("no observed candidate genes")
    observed = Counter()
    groups = defaultdict(list)
    for row in selected:
        key = row["gene_id"], row["from"], row["to"]
        if key in observed or not row["count"].isascii() or not row["count"].isdigit() or int(row["count"]) < 1:
            raise ValueError("duplicate edge or invalid raw count")
        observed[key] = int(row["count"])
        groups[row["gene_id"]].append(row)
    output.mkdir(parents=True)
    record = {
        "status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "SRR36710387 control; first 1,000,000 cleaned pairs (2.02% of one sample)",
        "node_rule": "historical exon-interval", "alignment_strandedness": "historical unstranded",
        "minimum_anchor": 8, "minimum_unique_support": 1, "exclude_duplicates": False,
        "inputs": inputs, "binary": identity(args.binary.resolve()),
        "script": identity(Path(__file__)), "prior_verification": identity(previous / "verification.json"),
        "numpy_version": np.__version__, "commands": [],
        "limitations": ["One control and no cases; no disease effect or p-value is calculated.",
                        "A deterministic prefix is not a full sample or random subsample.",
                        "Historical ambiguity exclusions remain; no discarded reads are rescued.",
                        "This does not validate the adopted RF/shared-splice-site cohort workflow.",
                        "Independent numerical checks and recount reuse the same underlying reads/reference."]}

    def call(command, name):
        record["commands"].append(list(map(str, command)))
        with (output / name).open("x") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=180)

    try:
        write_tsv(output / "junctions.tsv", ["sample_id", "gene_id", "from", "to", "count"],
                  ([row[name] for name in ["sample_id", "gene_id", "from", "to", "count"]] for row in selected))
        call([str(args.binary.resolve()), "validate", str(output / "junctions.tsv")], "validate.log")
        call([str(args.binary.resolve()), "analyze", str(output / "junctions.tsv"), "--output", str(output / "rust")], "analyze.log")
        call([str(args.samtools.resolve()), "quickcheck", "-v", str(paths["bam"])], "bam_integrity.log")

        # Resolve ambiguity across ALL genes before choosing candidate targets.
        # Filtering annotations to candidates first would falsely rescue reads.
        targets = {}
        with paths["gene_map"].open() as handle:
            def junction(row):
                return row["chromosome"], int(row["intron_start"]), int(row["intron_end"]), row["strand"]
            for key, rows in itertools.groupby(csv.DictReader(handle, delimiter="\t"), junction):
                mappings = {(row["gene_id"], row["from_node"], row["to_node"]) for row in rows}
                if len(mappings) == 1:
                    mapping = next(iter(mappings))
                    if mapping[0] in family and mapping[1] != mapping[2]:
                        targets[key] = mapping
        command = [str(args.samtools.resolve()), "view", "-F", "2820", str(paths["bam"])]
        record["commands"].append(command)
        with (output / "recount.log").open("x") as log:
            with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log, text=True) as process:
                try:
                    recounted = independent_counts(process.stdout, targets, 8)
                except BaseException:
                    process.kill()
                    process.wait()
                    raise
                if process.wait(timeout=30):
                    raise ValueError("BAM recount failed; see recount.log")
        differences = [(key, observed[key], recounted[key]) for key in sorted(observed.keys() | recounted.keys()) if observed[key] != recounted[key]]
        write_tsv(output / "recount_differences.tsv", ["gene_id", "from", "to", "exported", "recounted"],
                  ((*key, before, after) for key, before, after in differences))
        if differences:
            raise ValueError(f"{len(differences)} BAM count differences; inspect recount_differences.tsv")

        metrics = {row["gene_id"]: row for row in read_tsv(output / "rust/sample_metrics.tsv")}
        exported = defaultdict(dict)
        for row in read_tsv(output / "rust/junction_components.tsv"):
            exported[row["gene_id"]][row["from"], row["to"]] = row
        if set(metrics) != set(groups):
            raise ValueError("Rust did not return exactly the observed candidate genes")
        validations = []
        for gene, rows in sorted(groups.items()):
            check, components = projection(rows)
            if set(exported[gene]) != set(components):
                raise ValueError(f"edge mismatch for {gene}")
            maximum = max(abs(float(exported[gene][edge][name]) - calculated)
                          for edge, values in components.items()
                          for name, calculated in zip(["normalized_signal", "gradient_component", "cycle_component"], values))
            delta = abs(float(metrics[gene]["cycle_fraction"]) - check["cycle_fraction"])
            passed = metrics[gene]["checks_status"] == "PASS" and max(maximum, delta, check["reconstruction_error"], check["kernel_error"], check["orthogonality_error"]) <= 1e-8
            validations.append({"gene_id": gene, "gene_symbol": family[gene]["gene_symbol"],
                                **check, "maximum_component_difference": maximum,
                                "cycle_fraction_difference": delta, "status": "PASS" if passed else "FAIL"})
        header = list(validations[0])
        write_tsv(output / "independent_checks.tsv", header, ([row[name] for name in header] for row in validations))
        write_tsv(output / "candidate_coverage.tsv", ["gene_id", "gene_symbol", "status"],
                  ((gene, row["gene_symbol"], "observed" if gene in groups else "not_observed_in_prefix") for gene, row in sorted(family.items())))
        if any(row["status"] != "PASS" for row in validations):
            raise ValueError("independent SVD disagreement; inspect independent_checks.tsv")

        graphs = []
        for label, cyclic in [("cyclic", True), ("acyclic", False)]:
            eligible = [row for row in validations if (row["cycle_rank"] > 0) == cyclic]
            if eligible:
                gene = sorted(eligible, key=lambda row: (-row["total_support"], row["gene_id"]))[0]
                name = f"{label}.svg"
                call([str(args.binary.resolve()), "visualize", "SRR36710387", gene["gene_id"], str(output / "junctions.tsv"), str(output / name)], f"{label}_graph.log")
                graphs.append({"path": name, "gene_id": gene["gene_id"], "gene_symbol": gene["gene_symbol"], "selection": f"highest support among {label} graphs; gene ID breaks ties"})
        record.update(status="complete", candidate_genes=len(family), observed_genes=len(groups),
                      unobserved_genes=len(family) - len(groups), junction_rows=len(selected),
                      recounted_observations=sum(recounted.values()), recount_differences=0,
                      independent_checks_passed=len(validations),
                      cyclic_genes=sum(row["cycle_rank"] > 0 for row in validations),
                      maximum_component_difference=max(row["maximum_component_difference"] for row in validations),
                      maximum_cycle_fraction_difference=max(row["cycle_fraction_difference"] for row in validations),
                      example_graphs=graphs, biological_comparison_performed=False)
        lines = "".join(f'<tr><td>{html.escape(row["gene_symbol"])}</td><td>{html.escape(row["gene_id"])}</td><td>{row["total_support"]}</td><td>{row["cycle_rank"]}</td><td>{row["cycle_fraction"]:.6f}</td><td>{row["status"]}</td></tr>' for row in validations)
        links = " ".join(f'<a href="{row["path"]}">{html.escape(row["gene_symbol"])} graph</a>' for row in graphs)
        (output / "report.html").write_text(f'<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Splice Girl: real control validation</title><style>body{{max-width:1000px;margin:3rem auto;padding:0 1rem;font:18px Georgia,serif;color:#222}}table{{border-collapse:collapse;font-size:14px;width:100%}}td,th{{border-bottom:1px solid #ccc;padding:.5rem;text-align:left}}a{{color:#164f80}}</style><h1>Splice Girl: real control validation</h1><p>SRR36710387: first 1,000,000 cleaned read pairs; one control, no cases. Historical exon-interval nodes and unstranded alignment.</p><p><strong>{len(groups)} observed candidate genes passed independent SVD checks.</strong> All {len(selected)} exported junction edges matched the independent BAM recount ({sum(recounted.values())} read-junction observations).</p><p>{record["cyclic_genes"]} observed graphs have nonzero cycle rank. The {record["unobserved_genes"]} other candidates were not observed in this prefix; they are not assigned zero cycle fractions.</p><p>This is a technical validation, not evidence of a disease effect. It does not validate the final RF/shared-splice-site cohort protocol.</p><p>{links} · <a href="independent_checks.tsv">Numerical checks</a> · <a href="candidate_coverage.tsv">Coverage audit</a> · <a href="validation.json">Provenance</a></p><h2>All observed candidate genes</h2><table><thead><tr><th>Symbol</th><th>Gene ID</th><th>Support</th><th>Cycle rank</th><th>Cycle fraction</th><th>Checks</th></tr></thead><tbody>{lines}</tbody></table></html>')
        # Confirm that the source files remained unchanged during this check.
        for role, path in paths.items():
            if identity(path) != inputs[role]:
                raise ValueError(f"{role} changed during validation")
        record["outputs"] = [identity(path) for path in sorted(output.rglob("*")) if path.is_file()]
    except BaseException as error:
        record.update(status="failed", error=str(error))
        raise
    finally:
        (output / "validation.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({key: record[key] for key in ["status", "observed_genes", "junction_rows", "recounted_observations", "cyclic_genes", "maximum_component_difference"]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--samtools", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    run(parser.parse_args())
