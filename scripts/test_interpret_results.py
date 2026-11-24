import contextlib
import csv
import io
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from interpret_results import (
    ROOT,
    classify,
    expression_check,
    expression_log2,
    load_evidence,
    main,
    overlap_comparison,
    set_permutation,
    spearman,
    usage_entropy,
)

GENE_SETS = ROOT / "data/gene_sets"
BINARY = ROOT / "target/release/splice-girl"
EXAMPLE = ROOT / "examples/shared_graph"
# The synthetic triangle gene, copied onto one real candidate per category.
STAND_INS = {
    "ENSG00000001617.14": "risk_only",  # SEMA3F
    "ENSG00000001084.15": "environment_only",  # GCLC
    "ENSG00000082175.17": "overlap",  # PGR
}


def rows(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class EvidenceTests(unittest.TestCase):
    def test_committed_tables_agree_with_candidate_flags(self):
        candidates, risk, environment = load_evidence(GENE_SETS)
        self.assertEqual(len(candidates), 661)
        self.assertEqual(len(risk), 103)
        self.assertEqual(len(environment), 565)

    def test_a_gene_that_lost_its_evidence_stops_the_join(self):
        with tempfile.TemporaryDirectory() as folder:
            copy = Path(folder) / "gene_sets"
            shutil.copytree(GENE_SETS, copy)
            path = copy / "endometriosis_risk.tsv"
            lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
            path.write_text("".join(line for line in lines if "SEMA3F" not in line), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "SEMA3F: is_risk=1 but endometriosis_risk.tsv"):
                load_evidence(copy)


class MethodComparisonTests(unittest.TestCase):
    def test_usage_entropy_runs_from_one_edge_to_even_spread(self):
        self.assertEqual(usage_entropy([10, 0, 0]), 0)
        self.assertAlmostEqual(usage_entropy([5, 5, 5]), 1.0)
        self.assertIsNone(usage_entropy([0, 0, 0]))
        self.assertIsNone(usage_entropy([7]))

    def test_set_permutation_matches_the_rust_rules(self):
        is_case = [True] * 4 + [False] * 4
        perfect = set_permutation([[1, 1, 1, 1, 0, 0, 0, 0]], is_case)
        self.assertEqual((perfect["extreme_count"], perfect["permutations_considered"]), (2, 70))
        # Opposite directions do not cancel in the two-sided set statistic.
        opposite = set_permutation(
            [[1, 1, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1]], is_case
        )
        self.assertEqual(opposite["set_statistic"], 1.0)
        self.assertEqual(opposite["signed_effect_case_minus_control"], 0.0)
        self.assertAlmostEqual(opposite["p_value"], 2 / 70)
        flat = set_permutation([[0.5] * 8], is_case)
        self.assertEqual(flat["p_value"], 1.0)

    def test_outcome_rules(self):
        self.assertEqual(classify(None, None)["outcome"], "not_tested")
        self.assertEqual(classify(0.5, 0.01)["outcome"], "null")
        self.assertEqual(classify(0.05, None)["outcome"], "null")
        self.assertEqual(classify(2 / 70, 2 / 70)["outcome"], "difference_also_in_baseline")
        self.assertEqual(classify(2 / 70, 0.4)["outcome"], "difference_not_in_baseline")

    def test_spearman_handles_ties_and_small_inputs(self):
        self.assertAlmostEqual(spearman([1, 2, 3, 4], [10, 20, 30, 40]), 1.0)
        self.assertAlmostEqual(spearman([1, 2, 3, 4], [4, 3, 2, 1]), -1.0)
        self.assertIsNone(spearman([1, 2], [1, 2]))
        self.assertIsNone(spearman([1, 1, 1], [1, 2, 3]))


@unittest.skipUnless(BINARY.is_file(), "needs cargo build --release")
class AmendmentRuleTests(unittest.TestCase):
    def test_overlap_comparison_needs_three_tested_overlap_genes(self):
        effects = {"a": ("overlap", 0.5), "b": ("overlap", 0.4), "c": ("risk_only", 0.1)}
        result = overlap_comparison(effects)
        self.assertEqual(result["verdict"], "not_assessable")
        self.assertIsNone(result["p_value"])

    def test_overlap_comparison_is_exact_and_one_sided(self):
        effects = {f"o{i}": ("overlap", -0.9 + i * 0.01) for i in range(3)}
        effects.update({f"r{i}": ("risk_only", 0.01 * i) for i in range(4)})
        effects.update({f"e{i}": ("environment_only", 0.02) for i in range(3)})
        result = overlap_comparison(effects)
        # The three largest |effects| are the overlap genes: 1 of C(10, 3) = 120 sets.
        self.assertEqual(result["method"], "exact")
        self.assertEqual(result["permutations_considered"], 120)
        self.assertAlmostEqual(result["p_value"], 1 / 120)
        self.assertEqual(result["verdict"], "consistent_with_hypothesis")
        # Smaller overlap effects are not evidence for the hypothesis (one-sided).
        flipped = {g: (c, 0.001 if c == "overlap" else e) for g, (c, e) in effects.items()}
        self.assertEqual(overlap_comparison(flipped)["verdict"], "not_supported")

    def test_overlap_comparison_samples_large_pools_reproducibly(self):
        effects = {f"o{i}": ("overlap", 0.3) for i in range(7)}
        effects.update({f"e{i}": ("environment_only", (i % 10) / 100) for i in range(300)})
        first, second = overlap_comparison(effects), overlap_comparison(effects)
        self.assertTrue(first["method"].startswith("monte_carlo_seed_"))
        self.assertEqual(first["p_value"], second["p_value"])
        self.assertLess(first["p_value"], 0.001)

    def test_expression_rules(self):
        self.assertAlmostEqual(expression_log2(100, 1e8), 0.0)
        self.assertIsNone(expression_log2(0, 1e8))
        quiet = expression_check("null", 0.01, 0.1)
        self.assertEqual(quiet["required_statements"], [])
        both = expression_check("difference_not_in_baseline", 0.01, 0.6)
        self.assertEqual(len(both["required_statements"]), 2)


class InterpretRunTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        root = Path(cls.temporary.name)
        lines = (EXAMPLE / "junctions.tsv").read_text(encoding="utf-8").splitlines()
        body = [line.split("\t") for line in lines[1:]]
        out = [lines[0]]
        for number, gene_id in enumerate(STAND_INS):
            for fields in body:
                if fields[1] == "TRIANGLE":
                    # Vary the counts a little so the three genes differ.
                    count = int(fields[4]) + number * (1 if fields[0] in ("S0", "S1") else 0)
                    out.append("\t".join([fields[0], gene_id, fields[2], fields[3], str(count)]))
        (root / "junctions.tsv").write_text("\n".join(out) + "\n", encoding="utf-8")
        shutil.copy(EXAMPLE / "samples.tsv", root / "samples.tsv")
        config = json.loads((EXAMPLE / "analysis.json").read_text(encoding="utf-8"))
        config.update(
            genes=str(GENE_SETS / "candidate_genes.tsv"),
            output="results",
        )
        config["statistics"]["family"] = "candidate_genes"
        (root / "analysis.json").write_text(json.dumps(config), encoding="utf-8")
        subprocess.run(
            [str(BINARY), "run", "--config", str(root / "analysis.json")],
            check=True, capture_output=True,
        )
        cls.results = root / "results"
        cls.output = root / "interpretation"
        with contextlib.redirect_stdout(io.StringIO()) as cls.printed:
            cls.code = main(["--results", str(cls.results), "--binary", str(BINARY),
                             "--output", str(cls.output)])

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_every_candidate_has_one_row_with_its_evidence(self):
        self.assertEqual(self.code, 0)
        genes = rows(self.output / "gene_interpretation.tsv")
        self.assertEqual(len(genes), 661)
        self.assertEqual(len({g["gene_id"] for g in genes}), 661)
        by_id = {g["gene_id"]: g for g in genes}
        for gene_id, category in STAND_INS.items():
            gene = by_id[gene_id]
            self.assertEqual(gene["category"], category)
            self.assertEqual(gene["test_status"], "tested")
            self.assertTrue(gene["mean_difference_case_minus_control"])
            self.assertEqual(bool(gene["risk_evidence_ids"]), category != "environment_only")
            self.assertEqual(bool(gene["environment_evidence_ids"]), category != "risk_only")
        untested = [g for g in genes if g["gene_id"] not in STAND_INS]
        self.assertTrue(all(g["test_status"] == "not_tested" for g in untested))

    def test_set_tests_label_one_confirmatory_and_three_exploratory(self):
        sets = rows(self.output / "set_tests.tsv")
        self.assertEqual([(s["role"], s["set"]) for s in sets], [
            ("confirmatory", "candidate"), ("exploratory", "risk"),
            ("exploratory", "environment"), ("exploratory", "overlap"),
        ])
        # Overlap genes are in both the risk and environment sets.
        self.assertEqual([s["genes_used"] for s in sets], ["3", "2", "2", "1"])
        self.assertTrue(all(s["case_n"] == "4" and s["control_n"] == "4" for s in sets))
        self.assertIn("confirmatory candidate", self.printed.getvalue())

    def test_category_summary_and_junction_links(self):
        summary = {r["category"]: r for r in rows(self.output / "category_summary.tsv")}
        self.assertEqual(summary["all"]["candidate_genes"], "661")
        self.assertEqual(summary["all"]["tested"], "3")
        for category in STAND_INS.values():
            self.assertEqual(summary[category]["tested"], "1")
        links = rows(self.output / "gene_junctions.tsv")
        self.assertEqual({r["gene_id"] for r in links}, set(STAND_INS))
        self.assertEqual({r["condition"] for r in links}, {"case", "control"})
        discussed = rows(self.output / "discussed_genes.tsv")
        self.assertEqual([d["rank"] for d in discussed], ["1", "2", "3"])
        effects = [abs(float(d["mean_difference_case_minus_control"])) for d in discussed]
        self.assertEqual(effects, sorted(effects, reverse=True))
        manifest = json.loads((self.output / "interpretation_manifest.json").read_text())
        self.assertEqual(len(manifest["commands"]), 3)

    def test_junction_shares_add_up_to_each_discussed_gene_effect(self):
        discussed = {d["gene_id"]: d for d in rows(self.output / "discussed_genes.tsv")}
        shares = rows(self.output / "junction_contributions.tsv")
        self.assertEqual({r["gene_id"] for r in shares}, set(discussed))
        for gene_id, gene in discussed.items():
            edges = [r for r in shares if r["gene_id"] == gene_id]
            self.assertEqual(len(edges), 3)  # the triangle's three shared edges
            self.assertEqual([r["rank_within_gene"] for r in edges], ["1", "2", "3"])
            total = sum(float(r["cycle_energy_share_difference"]) for r in edges)
            self.assertAlmostEqual(total, float(gene["mean_difference_case_minus_control"]), places=8)

    def test_overlap_comparison_is_not_assessable_with_one_overlap_gene(self):
        overlap = rows(self.output / "overlap_comparison.tsv")[0]
        self.assertEqual(overlap["verdict"], "not_assessable")
        self.assertEqual(overlap["tested_overlap"], "1")
        self.assertIn("Overlap comparison (exploratory): not assessable", self.printed.getvalue())

    def test_method_comparison_reproduces_primary_and_records_outcome(self):
        comparison = {r["metric"]: r for r in rows(self.output / "method_comparison.tsv")}
        primary = rows(self.output / "set_tests.tsv")[0]
        hodge = comparison["cycle_fraction (primary)"]
        self.assertEqual(hodge["genes_used"], "3")
        self.assertEqual(hodge["extreme_count"], primary["extreme_count"])
        self.assertEqual(comparison["junction_usage_entropy (baseline)"]["genes_used"], "3")
        self.assertIn("spearman_abs_cycle_fraction_vs_abs_entropy_effect", comparison)
        self.assertIn("spearman_abs_cycle_fraction_vs_abs_expression_effect", comparison)
        self.assertEqual(comparison["gene_junction_expression_log2 (descriptive)"]["genes_used"], "3")
        per_gene = rows(self.output / "gene_method_comparison.tsv")
        self.assertEqual(len(per_gene), 3)
        self.assertTrue(all(r["expression_log2_difference"] for r in per_gene))
        outcome = json.loads((self.output / "outcome.json").read_text())
        self.assertIn(outcome["outcome"], {
            "null", "difference_also_in_baseline", "difference_not_in_baseline"})
        self.assertEqual(outcome["rules"], "docs/interpretation_plan.md")
        self.assertIn("required_statements", outcome["expression_check"])
        self.assertEqual(outcome["overlap_comparison"]["verdict"], "not_assessable")
        self.assertIn("Pre-declared outcome", self.printed.getvalue())

    def test_changed_category_in_results_is_refused(self):
        statistics = self.results / "analysis/primary/statistics/statistics.tsv"
        original = statistics.read_text(encoding="utf-8")
        try:
            statistics.write_text(original.replace("\trisk_only", "\toverlap", 1), encoding="utf-8")
            with contextlib.redirect_stderr(io.StringIO()) as err:
                code = main(["--results", str(self.results), "--binary", str(BINARY),
                             "--output", str(self.output.parent / "again")])
            self.assertEqual(code, 1)
            self.assertIn("candidate table says 'risk_only'", err.getvalue())
        finally:
            statistics.write_text(original, encoding="utf-8")


if __name__ == "__main__":
    unittest.main()
