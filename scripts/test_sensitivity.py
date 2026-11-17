import argparse
import csv
import json
import tempfile
import unittest
from pathlib import Path

from sensitivity import ROOT, retained_rows, run, thresholds


class SensitivityTests(unittest.TestCase):
    def test_threshold_grid_is_explicit_finite_and_coverage_uses_retained_edges(self):
        self.assertEqual(thresholds("3,1,3"), [1, 3])
        for invalid in ["nan", "-1", "inf", ""]:
            with self.assertRaises(ValueError):
                thresholds(invalid)
        rows = [
            dict(sample_id="S", gene_id="G", count=str(value)) for value in [10, 1, 0]
        ]
        self.assertEqual(retained_rows(rows, 2, 11), [])
        self.assertEqual(len(retained_rows(rows, 2, 10)), 1)

    def test_real_cli_recomputes_demo_and_keeps_disappearing_genes_missing(self):
        with tempfile.TemporaryDirectory() as directory:
            args = argparse.Namespace(
                junctions=ROOT / "examples/demo/all_genes.tsv",
                manifest=ROOT / "examples/demo/manifest.tsv",
                output_dir=Path(directory) / "sensitivity",
                minimum_junction_counts="1,100000",
                minimum_gene_support="0",
                binary=ROOT / "target/release/splice-girl",
            )
            run(args)
            with (args.output_dir / "sensitivity.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 4)
            self.assertEqual(
                [row["status"] for row in rows[:2]], ["complete", "complete"]
            )
            self.assertTrue(
                all(
                    row["status"] == "no_retained_junctions"
                    and row["mean_difference_case_minus_control"] == ""
                    for row in rows[2:]
                )
            )
            self.assertEqual(
                json.loads((args.output_dir / "run_manifest.json").read_text())[
                    "status"
                ],
                "complete",
            )
            with self.assertRaisesRegex(ValueError, "already exists"):
                run(args)


if __name__ == "__main__":
    unittest.main()
