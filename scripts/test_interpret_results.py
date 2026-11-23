import contextlib
import csv
import io
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from interpret_results import ROOT, load_evidence, main

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


@unittest.skipUnless(BINARY.is_file(), "needs cargo build --release")
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
        manifest = json.loads((self.output / "interpretation_manifest.json").read_text())
        self.assertEqual(len(manifest["commands"]), 3)

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
