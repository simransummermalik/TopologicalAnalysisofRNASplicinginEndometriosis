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
  method_comparison.tsv    the Hodge cycle-fraction set test beside the same
                           test on a junction-usage entropy baseline and on
                           depth-adjusted gene junction expression (amendment
                           A3), with gene-level rank correlations
  gene_method_comparison.tsv  per tested gene: all three effects and cycle rank
  discussed_genes.tsv      the genes the write-up may name: the five tested
                           genes with the largest absolute effect
  junction_contributions.tsv  per shared edge of each discussed gene: case and
                           control cycle-energy shares; the share differences
                           add up to the gene's effect (amendment A2)
  overlap_comparison.tsv   whether overlap genes have larger effects than the
                           other tested candidates (amendment A1, exploratory)
  outcome.json             which pre-declared outcome the results fall in
                           (docs/interpretation_plan.md)
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
import itertools
import json
import math
import random
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CATEGORIES = ("risk_only", "environment_only", "overlap")
SET_FAMILIES = (("risk", "is_risk"), ("environment", "is_environment"), ("overlap", "is_overlap"))
TOLERANCE = 1e-9  # same as src/hodge.rs
ALPHA = 0.05
DISCUSSED = 5
# Amendment of 2026-11-24 (docs/interpretation_plan.md section 6).
OVERLAP_MINIMUM_TESTED = 3
OVERLAP_EXACT_LIMIT = 100_000
OVERLAP_DRAWS = 100_000
OVERLAP_SEED = 20261124
STRONG_CORRELATION = 0.5
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


def set_permutation(values: list[list[float]], is_case: list[bool]) -> dict:
    """Two-sided exact set test, identical to src/permutation.rs.

    The statistic is the mean over genes of |mean(case) - mean(control)|;
    every relabeling with the observed number of cases is applied to all
    genes at once.
    """
    n, k = len(is_case), sum(is_case)

    def effect(cases: set[int], gene: list[float]) -> float:
        case = [gene[i] for i in range(n) if i in cases]
        control = [gene[i] for i in range(n) if i not in cases]
        return statistics.fmean(case) - statistics.fmean(control)

    observed_cases = {i for i in range(n) if is_case[i]}
    observed = statistics.fmean(abs(effect(observed_cases, g)) for g in values)
    signed = statistics.fmean(effect(observed_cases, g) for g in values)
    tolerance = TOLERANCE * max(abs(observed), 1.0)
    relabelings = list(itertools.combinations(range(n), k))
    extreme = sum(
        statistics.fmean(abs(effect(set(c), g)) for g in values) >= observed - tolerance
        for c in relabelings
    )
    return {
        "genes_used": len(values),
        "set_statistic": observed,
        "signed_effect_case_minus_control": signed,
        "p_value": extreme / len(relabelings),
        "extreme_count": extreme,
        "permutations_considered": len(relabelings),
    }


def usage_entropy(counts: list[float]) -> float | None:
    """Normalized Shannon entropy of read shares over a gene's shared edges.

    0 means every read is on one edge, 1 means reads are spread evenly. It
    uses the same edges and counts as the cycle fraction but no topology.
    """
    total = sum(counts)
    if total <= 0 or len(counts) < 2:
        return None
    shares = [c / total for c in counts if c > 0]
    return -sum(s * math.log(s) for s in shares) / math.log(len(counts))


def spearman(x: list[float], y: list[float]) -> float | None:
    def ranks(values):
        order = sorted(range(len(values)), key=values.__getitem__)
        result = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            for position in order[i : j + 1]:
                result[position] = (i + j) / 2 + 1
            i = j + 1
        return result

    if len(x) < 3:
        return None
    rx, ry = ranks(x), ranks(y)
    try:
        return statistics.correlation(rx, ry)
    except statistics.StatisticsError:
        return None


def sample_depths(junctions: Path) -> dict[str, float]:
    """Total assigned junction reads per sample in the run's input table."""
    depth: dict[str, float] = defaultdict(float)
    with junctions.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        sample, count = header.index("sample_id"), header.index("count")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > count:
                depth[fields[sample]] += float(fields[count])
    return dict(depth)


def expression_log2(gene_reads: float, depth: float) -> float | None:
    """log2 junction reads per million assigned junction reads in the sample."""
    if gene_reads <= 0 or depth <= 0:
        return None
    return math.log2(gene_reads / depth * 1e6)


def overlap_comparison(effects: dict[str, tuple[str, float]]) -> dict:
    """Amendment A1: is the mean |effect| larger in overlap genes than in the rest?

    `effects` maps tested gene ID to (category, case-minus-control effect).
    One-sided, because the proposal's hypothesis names a direction. Gene
    labels are permuted: every set of tested genes the size of the overlap
    group is a possible "overlap" group. Exact when there are at most
    OVERLAP_EXACT_LIMIT such sets, otherwise OVERLAP_DRAWS seeded draws.
    """
    by_category = defaultdict(list)
    for category, effect in effects.values():
        by_category[category].append(abs(effect))
    means = {c: (statistics.fmean(v) if v else None) for c, v in
             ((c, by_category.get(c, [])) for c in CATEGORIES)}
    counts = {c: len(by_category.get(c, [])) for c in CATEGORIES}
    result = {"tested_" + c: counts[c] for c in CATEGORIES}
    result.update({"mean_abs_effect_" + c: means[c] for c in CATEGORIES})
    overlap = by_category.get("overlap", [])
    rest = by_category.get("risk_only", []) + by_category.get("environment_only", [])
    if len(overlap) < OVERLAP_MINIMUM_TESTED or not rest:
        result.update(verdict="not_assessable", difference=None, p_value=None, method="",
                      permutations_considered=0,
                      reason=f"fewer than {OVERLAP_MINIMUM_TESTED} overlap genes tested"
                      if len(overlap) < OVERLAP_MINIMUM_TESTED else "no other tested genes")
        return result
    pool = overlap + rest
    k, n = len(overlap), len(pool)
    total = sum(pool)

    def difference(chosen_sum: float) -> float:
        return chosen_sum / k - (total - chosen_sum) / (n - k)

    observed = difference(sum(overlap))
    tolerance = TOLERANCE * max(abs(observed), 1.0)
    if math.comb(n, k) <= OVERLAP_EXACT_LIMIT:
        sums = [sum(pool[i] for i in c) for c in itertools.combinations(range(n), k)]
        extreme = sum(difference(s) >= observed - tolerance for s in sums)
        p_value, considered, method = extreme / len(sums), len(sums), "exact"
    else:
        rng = random.Random(OVERLAP_SEED)
        extreme = sum(
            difference(sum(rng.sample(pool, k))) >= observed - tolerance
            for _ in range(OVERLAP_DRAWS)
        )
        p_value = (extreme + 1) / (OVERLAP_DRAWS + 1)
        considered, method = OVERLAP_DRAWS, f"monte_carlo_seed_{OVERLAP_SEED}"
    larger = all(
        means[c] is None or means["overlap"] > means[c] for c in ("risk_only", "environment_only")
    )
    result.update(
        difference=observed, p_value=p_value, method=method, permutations_considered=considered,
        verdict="consistent_with_hypothesis" if larger and p_value < ALPHA else "not_supported",
        reason="",
    )
    return result


def junction_contributions(primary: Path, discussed: list[dict]) -> list[dict]:
    """Amendment A2: where in each discussed gene the cycle-fraction change sits.

    A sample's cycle fraction is the sum over edges of
    cycle_component^2 / signal_energy, so the case-minus-control differences
    of these per-edge shares add up exactly to the gene's effect.
    """
    wanted = {g["gene_id"]: g for g in discussed}
    condition = {row["sample_id"]: row["condition"] for row in read_tsv(primary / "samples.tsv")}
    energy = {
        (row["sample_id"], row["gene_id"]): float(row["signal_energy"])
        for row in read_tsv(primary / "statistics/sample_metrics.tsv")
        if row["gene_id"] in wanted and row["status"] == "passed"
    }
    values = defaultdict(lambda: defaultdict(lambda: {"case": [], "control": []}))
    for row in read_tsv(primary / "statistics/junction_components.tsv"):
        key = (row["sample_id"], row["gene_id"])
        if key not in energy:
            continue
        cycle = float(row["cycle_component"])
        edge = values[row["gene_id"]][(row["from"], row["to"])]
        group = condition[row["sample_id"]]
        edge[group].append((float(row["normalized_signal"]), cycle, cycle * cycle / energy[key]))
        edge.setdefault("count_" + group, []).append(float(row["raw_count"]))

    rows = []
    for gene_id, gene in wanted.items():
        gene_rows = []
        for (start, end), edge in values[gene_id].items():
            means = {
                group: [statistics.fmean(v[i] for v in edge[group]) for i in range(3)]
                for group in ("case", "control")
            }
            gene_rows.append({
                "discussed_rank": gene["rank"],
                "gene_id": gene_id,
                "gene_symbol": gene["gene_symbol"],
                "from": start,
                "to": end,
                "case_mean_count": fmt(statistics.fmean(edge["count_case"])),
                "control_mean_count": fmt(statistics.fmean(edge["count_control"])),
                "case_mean_normalized_signal": fmt(means["case"][0]),
                "control_mean_normalized_signal": fmt(means["control"][0]),
                "case_mean_cycle_component": fmt(means["case"][1]),
                "control_mean_cycle_component": fmt(means["control"][1]),
                "case_mean_cycle_energy_share": fmt(means["case"][2]),
                "control_mean_cycle_energy_share": fmt(means["control"][2]),
                "cycle_energy_share_difference": means["case"][2] - means["control"][2],
            })
        total = sum(r["cycle_energy_share_difference"] for r in gene_rows)
        effect = float(gene["mean_difference_case_minus_control"])
        if not gene_rows or abs(total - effect) > 1e-8:
            raise ValueError(f"{gene_id}: junction shares do not add up to the gene's effect")
        gene_rows.sort(key=lambda r: (-abs(r["cycle_energy_share_difference"]), r["from"], r["to"]))
        for position, row in enumerate(gene_rows, 1):
            row["rank_within_gene"] = position
            row["cycle_energy_share_difference"] = fmt(row["cycle_energy_share_difference"])
        rows.extend(gene_rows)
    return rows


def classify(primary_p: float | None, baseline_p: float | None) -> dict:
    """The pre-declared outcome from docs/interpretation_plan.md."""
    if primary_p is None:
        return {"outcome": "not_tested",
                "statement": "The primary candidate-set test could not be run; no case/control conclusion."}
    if primary_p >= ALPHA:
        return {"outcome": "null",
                "statement": "No difference in cycle fraction between endometriosis and control samples "
                             "was detected across the candidate genes (primary test)."}
    if baseline_p is not None and baseline_p < ALPHA:
        return {"outcome": "difference_also_in_baseline",
                "statement": "Cycle fraction differed between groups across the candidate genes, and so "
                             "did junction-usage entropy; the Hodge summary is not shown to add "
                             "information beyond junction usage."}
    return {"outcome": "difference_not_in_baseline",
            "statement": "Cycle fraction differed between groups across the candidate genes, while the "
                         "junction-usage entropy baseline did not reach p < 0.05 (exploratory "
                         "comparison; it does not show the Hodge summary is better in general)."}


def expression_check(outcome: str, expression_p: float | None, rho: float | None) -> dict:
    """Amendment A3: the wording required about expression (interpretation plan section 6)."""
    notes = []
    if outcome.startswith("difference") and expression_p is not None and expression_p < ALPHA:
        notes.append("Gene junction expression also differed across the candidate genes "
                     "(exploratory); the cycle-fraction difference is not shown to be separate "
                     "from an expression difference.")
    if rho is not None and rho >= STRONG_CORRELATION:
        notes.append("Genes with larger cycle-fraction effects also tended to have larger "
                     "expression effects; the cycle-fraction effects largely track expression.")
    return {"spearman": rho, "strong_correlation_threshold": STRONG_CORRELATION,
            "required_statements": notes}


def compare_methods(primary: Path, output: Path, primary_set: dict, depth: dict[str, float]) -> dict:
    """Hodge cycle fraction versus the junction-usage baseline and expression, same genes and samples."""
    samples = read_tsv(primary / "samples.tsv")
    order = [row["sample_id"] for row in samples]
    is_case = [row["condition"] == "case" for row in samples]
    stats_rows = [r for r in read_tsv(primary / "statistics/statistics.tsv") if r["status"] == "tested"]
    tested = {r["gene_id"]: r for r in stats_rows}
    cycle = defaultdict(dict)
    for row in read_tsv(primary / "statistics/sample_metrics.tsv"):
        if row["gene_id"] in tested and row["status"] == "passed":
            cycle[row["gene_id"]][row["sample_id"]] = float(row["cycle_fraction"])
    edges = defaultdict(lambda: defaultdict(dict))
    for row in read_tsv(primary / "junctions.tsv"):
        if row["gene_id"] in tested:
            edges[row["gene_id"]][(row["from"], row["to"])][row["sample_id"]] = float(row["count"])
    rank = {r["gene_id"]: r["cycle_rank"] for r in read_tsv(primary / "gene_eligibility.tsv")}

    missing_depth = [s for s in order if depth.get(s, 0) <= 0]
    if missing_depth:
        raise ValueError(f"no junction reads in the input table for: {', '.join(missing_depth)}")
    genes, hodge, baseline, expression = [], [], [], []
    for gene_id in sorted(tested):
        if len(cycle[gene_id]) != len(order):
            continue
        entropy, level = [], []
        for sample in order:
            # Every selected shared edge is present for a covered sample
            # (a measured zero is written as 0).
            counts = [by_sample.get(sample) for by_sample in edges[gene_id].values()]
            entropy.append(None if None in counts else usage_entropy(counts))
            level.append(None if None in counts else expression_log2(sum(counts), depth[sample]))
        if None in entropy or None in level:
            continue
        genes.append(gene_id)
        hodge.append([cycle[gene_id][s] for s in order])
        baseline.append(entropy)
        expression.append(level)

    if not genes:
        result = {"hodge": None, "baseline": None, "expression": None}
    else:
        result = {
            "hodge": set_permutation(hodge, is_case),
            "baseline": set_permutation(baseline, is_case),
            "expression": set_permutation(expression, is_case),
        }
        rust = number(primary_set.get("set_statistic"))
        rust_p = number(primary_set.get("p_value"))
        mine = result["hodge"]
        # Cross-check: this reimplementation must reproduce the Rust primary test.
        if (
            str(mine["genes_used"]) != primary_set.get("genes_used")
            or rust is None
            or abs(mine["set_statistic"] - rust) > 1e-8
            or str(mine["extreme_count"]) != primary_set.get("extreme_count")
            or rust_p is None
            or abs(mine["p_value"] - rust_p) > 1e-9  # Rust writes 10 decimals
        ):
            raise ValueError("method comparison does not reproduce the Rust primary set test")

    def mean_diff(values):
        case = [v for v, c in zip(values, is_case) if c]
        control = [v for v, c in zip(values, is_case) if not c]
        return statistics.fmean(case) - statistics.fmean(control)

    gene_rows = [
        {
            "gene_id": gene_id,
            "gene_symbol": tested[gene_id]["gene_symbol"],
            "category": tested[gene_id]["category"],
            "cycle_rank": rank.get(gene_id, ""),
            "cycle_fraction_difference": fmt(mean_diff(h)),
            "usage_entropy_difference": fmt(mean_diff(b)),
            "expression_log2_difference": fmt(mean_diff(e)),
        }
        for gene_id, h, b, e in zip(genes, hodge, baseline, expression)
    ]
    write_tsv(output / "gene_method_comparison.tsv",
              ["gene_id", "gene_symbol", "category", "cycle_rank", "cycle_fraction_difference",
               "usage_entropy_difference", "expression_log2_difference"], gene_rows)

    abs_cf = [abs(float(r["cycle_fraction_difference"])) for r in gene_rows]
    abs_ent = [abs(float(r["usage_entropy_difference"])) for r in gene_rows]
    abs_expr = [abs(float(r["expression_log2_difference"])) for r in gene_rows]
    ranks = [float(r["cycle_rank"]) for r in gene_rows]
    correlations = {
        "spearman_abs_cycle_fraction_vs_abs_entropy_effect": spearman(abs_cf, abs_ent),
        "spearman_abs_cycle_fraction_effect_vs_cycle_rank": spearman(abs_cf, ranks),
        "spearman_abs_cycle_fraction_vs_abs_expression_effect": spearman(abs_cf, abs_expr),
    }
    rows = []
    for metric, role in (("hodge", "cycle_fraction (primary)"),
                         ("baseline", "junction_usage_entropy (baseline)"),
                         ("expression", "gene_junction_expression_log2 (descriptive)")):
        test = result[metric] or {}
        rows.append({"metric": role, **{k: (fmt(v) if isinstance(v, float) else v) for k, v in test.items()}})
    for name, value in correlations.items():
        rows.append({"metric": name, "set_statistic": fmt(value), "genes_used": len(gene_rows)})
    write_tsv(output / "method_comparison.tsv",
              ["metric", "genes_used", "set_statistic", "signed_effect_case_minus_control",
               "p_value", "extreme_count", "permutations_considered"], rows)
    outcome = classify(
        result["hodge"]["p_value"] if result["hodge"] else None,
        result["baseline"]["p_value"] if result["baseline"] else None,
    )
    expression_p = result["expression"]["p_value"] if result["expression"] else None
    expression_rho = correlations["spearman_abs_cycle_fraction_vs_abs_expression_effect"]
    outcome.update(
        primary_p_value=result["hodge"]["p_value"] if result["hodge"] else None,
        baseline_p_value=result["baseline"]["p_value"] if result["baseline"] else None,
        expression_p_value=expression_p,
        expression_check=expression_check(outcome["outcome"], expression_p, expression_rho),
        alpha=ALPHA,
        rules="docs/interpretation_plan.md",
        **correlations,
    )
    (output / "outcome.json").write_text(json.dumps(outcome, indent=2, sort_keys=True) + "\n",
                                         encoding="utf-8")
    return outcome


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
    required.append(results / "resolved_configuration.json")
    required.append(primary / "statistics/junction_components.tsv")
    for path in required:
        if not path.is_file():
            raise ValueError(f"not a finished configured run: missing {path}")
    if output.exists():
        raise ValueError(f"output folder already exists: {output}; choose a new folder")

    configuration = json.loads((results / "resolved_configuration.json").read_text(encoding="utf-8"))
    input_junctions = Path(configuration["junctions"])
    if not input_junctions.is_file():
        raise ValueError(f"the run's input junction table is missing: {input_junctions}")
    required.append(input_junctions)
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

    # docs/interpretation_plan.md section 4: chosen by rule, not by eye.
    ranked = sorted(
        (g for g in genes if g["test_status"] == "tested"),
        key=lambda g: (-abs(float(g["mean_difference_case_minus_control"])),
                       float(g["p_value_exploratory"]), g["gene_id"]),
    )
    discussed = [{"rank": i, **g} for i, g in enumerate(ranked[:DISCUSSED], 1)]
    write_tsv(output / "discussed_genes.tsv", ["rank"] + list(genes[0]), discussed)
    contributions = junction_contributions(primary, discussed)
    write_tsv(output / "junction_contributions.tsv",
              list(contributions[0]) if contributions else ["discussed_rank", "gene_id"], contributions)

    overlap = overlap_comparison({
        g["gene_id"]: (g["category"], float(g["mean_difference_case_minus_control"]))
        for g in genes if g["test_status"] == "tested"
    })
    write_tsv(output / "overlap_comparison.tsv", list(overlap),
              [{k: fmt(v) if isinstance(v, float) else ("" if v is None else v) for k, v in overlap.items()}])

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

    outcome = compare_methods(primary, output, primary_set, sample_depths(input_junctions))
    outcome["overlap_comparison"] = overlap
    (output / "outcome.json").write_text(json.dumps(outcome, indent=2, sort_keys=True) + "\n",
                                         encoding="utf-8")

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
    return {"summary": summary, "set_tests": set_rows, "outcome": outcome}


def outcome_line(name: str) -> str:
    return name.replace("_", " ")


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
    decided = outcome["outcome"]
    print(f"Pre-declared outcome: {outcome_line(decided['outcome'])}. {decided['statement']}")
    for note in decided["expression_check"]["required_statements"]:
        print(f"Also required: {note}")
    overlap = decided["overlap_comparison"]
    print(f"Overlap comparison (exploratory): {outcome_line(overlap['verdict'])}"
          + (f", p = {overlap['p_value']:.4f}" if overlap["p_value"] is not None else
             f" ({overlap['reason']})"))
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
