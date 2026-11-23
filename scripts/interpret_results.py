#!/usr/bin/env python3
"""Join a configured Splice Girl run to the candidate-gene evidence tables.

Reads a finished `splice-girl run --config` results folder and the evidence
tables in data/gene_sets, and writes a new folder with:

  gene_interpretation.tsv  one row for every candidate gene: its category,
                           evidence sources, eligibility, primary result and
                           how many sensitivity scenarios tested it
  gene_junctions.tsv       per-sample shared-graph edge counts for every tested
                           gene, with condition, symbol and category
  category_summary.tsv     descriptive counts and effects per category
  set_tests.tsv            the confirmatory candidate-set test plus the
                           exploratory risk, environment and overlap set tests
  category_tests/          the Rust statistics output behind each exploratory
                           test
  interpretation_manifest.json  inputs with checksums and the commands run

The script stops if any gene's category flags disagree with the evidence
tables, so no gene can silently lose its evidence. Exploratory set tests use
the primary run's own shared-graph counts and eligible genes, the same test
settings, and report raw p-values only (2026-11-17 inference rules).
Nothing here selects genes by their results.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CATEGORIES = ("risk_only", "environment_only", "overlap")
SET_FAMILIES = (("risk", "is_risk"), ("environment", "is_environment"), ("overlap", "is_overlap"))
NO_LOST_EVIDENCE = {
    "is_risk": "endometriosis_risk.tsv",
    "is_environment": "environment_response.tsv",
    "is_overlap": "overlap.tsv",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, columns: list[str], rows) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, columns, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> dict:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return {"path": str(path.resolve()), "sha256": digest}


def number(value: str | None) -> float | None:
    try:
        return float(value) if value not in (None, "") else None
    except ValueError:
        return None


def fmt(value: float | None) -> str:
    return "" if value is None else f"{value:.10f}"


def load_evidence(gene_sets: Path):
    candidates = {row["gene_id"]: row for row in read_tsv(gene_sets / "candidate_genes.tsv")}
    risk = {row["gene_id"]: row for row in read_tsv(gene_sets / "endometriosis_risk.tsv")}
    environment = defaultdict(list)
    for row in read_tsv(gene_sets / "environment_response.tsv"):
        environment[row["gene_id"]].append(row)
    overlap = {row["gene_id"] for row in read_tsv(gene_sets / "overlap.tsv")}
    members = {"is_risk": set(risk), "is_environment": set(environment), "is_overlap": overlap}
    problems = []
    for gene_id, row in candidates.items():
        for flag, table in NO_LOST_EVIDENCE.items():
            if (row[flag] == "1") != (gene_id in members[flag]):
                problems.append(f"{gene_id} {row['gene_symbol']}: {flag}={row[flag]} but {table} disagrees")
    for flag, table in NO_LOST_EVIDENCE.items():
        for gene_id in sorted(members[flag] - set(candidates)):
            problems.append(f"{gene_id}: in {table} but not a candidate gene")
    if problems:
        raise ValueError("evidence tables disagree with candidate flags:\n  " + "\n  ".join(problems))
    return candidates, risk, environment


def check_result_annotations(rows: list[dict[str, str]], candidates, table: str) -> None:
    """Every result row must carry the candidate table's flags unchanged."""
    for row in rows:
        source = candidates.get(row["gene_id"])
        if source is None:
            raise ValueError(f"{table}: {row['gene_id']} is not a candidate gene")
        for column in ("gene_symbol", "is_risk", "is_environment", "is_overlap", "category"):
            if column in row and row[column] != source[column]:
                raise ValueError(
                    f"{table}: {row['gene_id']} {column} is {row[column]!r}, "
                    f"candidate table says {source[column]!r}"
                )


def environment_summary(rows: list[dict[str, str]]) -> tuple[str, str, str]:
    exposures = sorted({row["exposure"] for row in rows})
    publications = sorted(
        f"{row['exposure']}:{row['n_publications']}" for row in rows
    )
    evidence = sorted({e for row in rows for e in row["evidence_ids"].split("|") if e})
    return "|".join(exposures), "|".join(publications), "|".join(evidence)


def run_set_test(binary: Path, primary: Path, family: Path, label: str, output: Path) -> list[str]:
    command = [
        str(binary), "stats", str(primary / "junctions.tsv"), str(primary / "samples.tsv"),
        "--genes", str(family), "--family", label, "--alternative", "two-sided",
        "--design", "independent", "--min-per-group", "4", "--output", str(output),
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise ValueError(f"{label} set test failed: {result.stderr.strip() or result.stdout.strip()}")
    return command


def interpret(results: Path, gene_sets: Path, binary: Path, output: Path) -> dict:
    primary = results / "analysis" / "primary"
    required = [
        primary / name
        for name in (
            "junctions.tsv", "samples.tsv", "eligible_genes.tsv", "gene_eligibility.tsv",
            "statistics/statistics.tsv", "statistics/sample_metrics.tsv", "statistics/set_test.tsv",
        )
    ] + [results / "analysis" / "sensitivity.tsv"]
    for path in required:
        if not path.is_file():
            raise ValueError(f"not a finished configured run: missing {path}")
    if output.exists():
        raise ValueError(f"output folder already exists: {output}; choose a new folder")

    candidates, risk, environment = load_evidence(gene_sets)
    eligible = read_tsv(primary / "eligible_genes.tsv")
    stats_rows = read_tsv(primary / "statistics/statistics.tsv")
    check_result_annotations(eligible, candidates, "eligible_genes.tsv")
    check_result_annotations(stats_rows, candidates, "statistics.tsv")

    eligibility = {row["gene_id"]: row for row in read_tsv(primary / "gene_eligibility.tsv")}
    tests = {row["gene_id"]: row for row in stats_rows}
    condition = {row["sample_id"]: row["condition"] for row in read_tsv(primary / "samples.tsv")}
    cycle = defaultdict(lambda: {"case": [], "control": []})
    for row in read_tsv(primary / "statistics/sample_metrics.tsv"):
        value = number(row.get("cycle_fraction"))
        if value is not None and row.get("status") == "passed":
            cycle[row["gene_id"]][condition[row["sample_id"]]].append(value)
    sensitivity = defaultdict(set)
    scenarios = set()
    for row in read_tsv(results / "analysis" / "sensitivity.tsv"):
        scenarios.add(row["scenario"])
        if row["status"] == "eligible" and row.get("case_n") == "4" and row.get("control_n") == "4":
            sensitivity[row["gene_id"]].add(row["scenario"])

    output.mkdir(parents=True)
    genes = []
    for gene_id, candidate in sorted(candidates.items(), key=lambda item: item[1]["gene_symbol"]):
        gate = eligibility.get(gene_id)
        test = tests.get(gene_id, {})
        means = {
            group: statistics.fmean(values) if values else None
            for group, values in cycle[gene_id].items()
        } if gene_id in cycle else {"case": None, "control": None}
        risk_row = risk.get(gene_id, {})
        exposures, publications, env_ids = environment_summary(environment.get(gene_id, []))
        genes.append({
            "gene_id": gene_id,
            "gene_symbol": candidate["gene_symbol"],
            "category": candidate["category"],
            "eligibility": gate["status"] if gate else "not_in_eligibility_table",
            "shared_edges": gate["shared_edges"] if gate else "",
            "cycle_rank": gate["cycle_rank"] if gate else "",
            "covered_samples": gate["covered_samples"] if gate else "",
            "test_status": test.get("status", "not_tested"),
            "case_mean_cycle_fraction": fmt(means["case"]),
            "control_mean_cycle_fraction": fmt(means["control"]),
            "mean_difference_case_minus_control": test.get("mean_difference_case_minus_control", ""),
            "p_value_exploratory": test.get("p_value", ""),
            "bh_adjusted_p_exploratory": test.get("adjusted_p_value", ""),
            "untested_reason": test.get("untested_reason", ""),
            "sensitivity_scenarios_tested": f"{len(sensitivity[gene_id])} of {len(scenarios)}",
            "risk_sources": risk_row.get("source", ""),
            "risk_evidence_ids": risk_row.get("evidence_ids", ""),
            "risk_notes": risk_row.get("notes", ""),
            "environment_exposures": exposures,
            "environment_publications": publications,
            "environment_evidence_ids": env_ids,
        })
    write_tsv(output / "gene_interpretation.tsv", list(genes[0]), genes)

    tested = {row["gene_id"] for row in stats_rows if row["status"] == "tested"}
    junction_rows = []
    for row in read_tsv(primary / "junctions.tsv"):
        if row["gene_id"] in tested:
            candidate = candidates[row["gene_id"]]
            junction_rows.append({
                **row,
                "condition": condition[row["sample_id"]],
                "gene_symbol": candidate["gene_symbol"],
                "category": candidate["category"],
            })
    write_tsv(
        output / "gene_junctions.tsv",
        ["gene_id", "gene_symbol", "category", "sample_id", "condition", "from", "to", "count"],
        sorted(junction_rows, key=lambda r: (r["gene_symbol"], r["from"], r["to"], r["condition"], r["sample_id"])),
    )

    summary = []
    for category in CATEGORIES + ("all",):
        rows = [g for g in genes if category == "all" or g["category"] == category]
        effects = [number(g["mean_difference_case_minus_control"]) for g in rows if g["test_status"] == "tested"]
        effects = [e for e in effects if e is not None]
        p_values = [number(g["p_value_exploratory"]) for g in rows if g["test_status"] == "tested"]
        summary.append({
            "category": category,
            "candidate_genes": len(rows),
            "with_shared_edges": sum(number(g["shared_edges"]) not in (None, 0) for g in rows),
            "eligible": sum(g["eligibility"] == "eligible" for g in rows),
            "tested": len(effects),
            "mean_signed_effect": fmt(statistics.fmean(effects) if effects else None),
            "mean_absolute_effect": fmt(statistics.fmean(abs(e) for e in effects) if effects else None),
            "median_absolute_effect": fmt(statistics.median(abs(e) for e in effects) if effects else None),
            "raw_p_below_0.05_exploratory": sum(p is not None and p < 0.05 for p in p_values),
        })
    write_tsv(output / "category_summary.tsv", list(summary[0]), summary)

    commands = []
    set_rows = []
    primary_set = read_tsv(primary / "statistics/set_test.tsv")[0]
    set_rows.append({"role": "confirmatory", "set": "candidate", **primary_set})
    families = output / "category_tests"
    families.mkdir()
    columns = list(eligible[0]) if eligible else ["gene_id"]
    for label, flag in SET_FAMILIES:
        members = [row for row in eligible if row[flag] == "1"]
        if not members:
            set_rows.append({"role": "exploratory", "set": label, "genes_used": "0",
                             "reason": "no eligible genes in this set"})
            continue
        family = families / f"{label}_family.tsv"
        write_tsv(family, columns, members)
        commands.append(run_set_test(binary, primary, family, label, families / label))
        row = read_tsv(families / label / "set_test.tsv")[0]
        set_rows.append({"role": "exploratory", "set": label, **row})
    set_columns = ["role", "set"] + [c for c in primary_set if c not in ("role", "set")]
    write_tsv(output / "set_tests.tsv", set_columns, set_rows)

    manifest = {
        "tool": "scripts/interpret_results.py",
        "inputs": [sha256(path) for path in required]
        + [sha256(gene_sets / name) for name in
           ("candidate_genes.tsv", "endometriosis_risk.tsv", "environment_response.tsv", "overlap.tsv")],
        "binary": str(binary),
        "commands": commands,
        "candidate_genes": len(candidates),
        "tested_genes": len(tested),
        "note": "Only the candidate-set test is confirmatory; every other p-value is exploratory.",
    }
    (output / "interpretation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return {"summary": summary, "set_tests": set_rows}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--results", type=Path, required=True, help="configured run results folder")
    parser.add_argument("--gene-sets", type=Path, default=ROOT / "data/gene_sets")
    parser.add_argument("--binary", type=Path, default=ROOT / "target/release/splice-girl")
    parser.add_argument("--output", type=Path, required=True, help="new folder")
    args = parser.parse_args(argv)
    try:
        if not args.binary.is_file():
            raise ValueError(f"splice-girl binary not found: {args.binary} (build with cargo build --release)")
        outcome = interpret(args.results, args.gene_sets, args.binary, args.output)
    except (OSError, ValueError, KeyError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    for row in outcome["set_tests"]:
        p = row.get("p_value") or "n/a"
        print(f"{row['role']:>12} {row['set']:<12} genes {row.get('genes_used', '0'):>4}  p = {p}")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
