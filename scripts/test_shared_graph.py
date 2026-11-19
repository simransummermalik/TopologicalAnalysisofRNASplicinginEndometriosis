"""Scientific acceptance checks for shared topology, missingness and CLI handoff."""

import argparse
import csv
import json
import hashlib
import subprocess
import tempfile
import unittest
from pathlib import Path

from shared_graph import (
    ROOT,
    build_shared,
    cycle_rank,
    load_cohort,
    prepare,
    read_table,
    write_table,
)
from sensitivity import run, shared_scenarios


class SharedGraphTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="splice-girl-shared-")
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.samples = [f"S{number}" for number in range(8)]
        self.manifest = self.root / "samples.tsv"
        write_table(
            self.manifest,
            ["sample_id", "condition", "subject_id"],
            [
                dict(
                    sample_id=sample,
                    condition="case" if i < 4 else "control",
                    subject_id=sample,
                )
                for i, sample in enumerate(self.samples)
            ],
        )
        self.genes = self.root / "family.tsv"
        write_table(
            self.genes,
            ["gene_id", "gene_symbol"],
            [
                dict(gene_id=gene, gene_symbol=gene + "_symbol")
                for gene in ("TRIANGLE", "TREE", "LOW_COVERAGE", "ABSENT")
            ],
        )
        self.rows = []
        for i, sample in enumerate(self.samples):
            for gene, source, target, count in [
                ("TRIANGLE", "A", "B", 30),
                ("TRIANGLE", "B", "C", 30),
                ("TRIANGLE", "A", "C", 5 if i < 4 else 0),
                ("TREE", "A", "B", 80),
                ("TREE", "B", "C", 80),
                ("LOW_COVERAGE", "A", "B", 20 if i < 7 else 1),
                ("LOW_COVERAGE", "B", "C", 20 if i < 7 else 0),
                ("LOW_COVERAGE", "A", "C", 10 if i < 7 else 0),
            ]:
                # Absent input rows deliberately exercise measured-zero insertion.
                if count:
                    self.rows.append(
                        dict(
                            sample_id=sample,
                            gene_id=gene,
                            **{"from": source, "to": target},
                            count=count,
                        )
                    )
        self.rows.append(
            dict(
                sample_id="S0",
                gene_id="LOW_COVERAGE",
                **{"from": "X", "to": "Y"},
                count=10000,
            )
        )
        self.junctions = self.root / "raw.tsv"
        self.write_rows(self.junctions, self.rows)

    def write_rows(self, path, rows):
        write_table(path, ["sample_id", "gene_id", "from", "to", "count"], rows)

    def cohort(self):
        return load_cohort(self.junctions, self.manifest, self.genes, True)

    def options(self):
        return argparse.Namespace(
            junctions=self.junctions,
            manifest=self.manifest,
            genes=self.genes,
            output_dir=self.root / "prepared",
            complete_counts=True,
            minimum_edge_count=5,
            minimum_samples_per_edge=4,
            minimum_gene_support=50,
        )

    def sensitivity_options(self):
        return argparse.Namespace(
            junctions=self.junctions,
            manifest=self.manifest,
            genes=self.genes,
            output_dir=self.root / "sensitivity",
            complete_counts=True,
            mode="shared",
            minimum_junction_counts="3,5,10",
            minimum_samples_per_edge="2,4,6",
            minimum_gene_support="20,50,100",
            primary_edge_count=5,
            primary_sample_count=4,
            primary_gene_support=50,
            binary=ROOT / "target/release/splice-girl",
            run_stats=True,
            family="synthetic-shared",
            alternative="two-sided",
            design="independent",
            min_per_group=4,
        )

    def test_fixed_topology_keeps_measured_zeros_and_coverage_uses_only_shared_edges(
        self,
    ):
        result = build_shared(self.cohort(), 5, 4, 50)
        triangle = [row for row in result["junctions"] if row["gene_id"] == "TRIANGLE"]
        self.assertEqual(len(triangle), 24)
        expected = {("A", "B"), ("B", "C"), ("A", "C")}
        for sample in self.samples:
            self.assertEqual(
                {
                    (row["from"], row["to"])
                    for row in triangle
                    if row["sample_id"] == sample
                },
                expected,
            )
        self.assertEqual(sum(row["count"] == 0 for row in triangle), 4)
        audit = {(row["gene_id"], row["sample_id"]): row for row in result["samples"]}
        self.assertEqual(audit[("LOW_COVERAGE", "S0")]["shared_graph_support"], 50)
        self.assertEqual(audit[("LOW_COVERAGE", "S7")]["status"], "below_coverage")
        self.assertFalse(
            any(
                row["sample_id"] == "S7" and row["gene_id"] == "LOW_COVERAGE"
                for row in result["junctions"]
            )
        )
        self.assertEqual([row["gene_id"] for row in result["eligible"]], ["TRIANGLE"])
        self.assertEqual(result["eligible"][0]["gene_symbol"], "TRIANGLE_symbol")
        status = {row["gene_id"]: row["status"] for row in result["genes"]}
        self.assertEqual(
            status,
            dict(
                ABSENT="no_shared_edges",
                LOW_COVERAGE="incomplete_coverage",
                TREE="no_cycle_in_shared_graph",
                TRIANGLE="eligible",
            ),
        )

    def test_filter_is_label_and_row_order_invariant(self):
        expected = build_shared(self.cohort(), 5, 4, 50)
        flipped = self.root / "flipped.tsv"
        write_table(
            flipped,
            ["sample_id", "condition"],
            [
                dict(sample_id=sample, condition="control" if i < 4 else "case")
                for i, sample in enumerate(self.samples)
            ],
        )
        reversed_input = self.root / "reversed.tsv"
        self.write_rows(reversed_input, list(reversed(self.rows)))
        actual = build_shared(
            load_cohort(reversed_input, flipped, self.genes, True), 5, 4, 50
        )
        self.assertEqual(expected, actual)

    def test_cycle_rank_handles_dags_disconnected_graphs_and_opposite_edges(self):
        self.assertEqual(cycle_rank([]), 0)
        self.assertEqual(cycle_rank([("A", "B"), ("B", "C")]), 0)
        # An exon-skipping DAG has a one-dimensional cycle space.
        self.assertEqual(
            cycle_rank([("A", "B"), ("B", "C"), ("A", "C"), ("D", "E")]), 1
        )
        self.assertEqual(cycle_rank([("A", "B"), ("B", "A")]), 1)

    def test_rejects_missing_exports_unacknowledged_zeros_and_invalid_raw_rows(self):
        with self.assertRaisesRegex(ValueError, "complete-counts"):
            load_cohort(self.junctions, self.manifest)
        for index, (change, message) in enumerate(
            [
                (
                    lambda rows: [row for row in rows if row["sample_id"] != "S7"],
                    "no positive",
                ),
                (lambda rows: rows + [rows[0]], "duplicate"),
                (lambda rows: [dict(rows[0], count="NaN")] + rows[1:], "integer"),
                (lambda rows: [dict(rows[0], count="0.5")] + rows[1:], "integer"),
                (
                    lambda rows: [dict(rows[0], sample_id="unknown")] + rows[1:],
                    "absent from manifest",
                ),
                (lambda rows: [dict(rows[0], to="A")] + rows[1:], "self-loop"),
            ]
        ):
            path = self.root / f"bad-{index}.tsv"
            self.write_rows(path, change(self.rows))
            with self.assertRaisesRegex(ValueError, message):
                load_cohort(path, self.manifest, self.genes, True)
        for values in ((0, 4, 50), (5, 9, 50), (5, 4, 0)):
            with self.assertRaises(ValueError):
                build_shared(self.cohort(), *values)

    def test_prepare_exports_provenance_and_protects_existing_outputs(self):
        args = self.options()
        prepare(args)
        record = json.loads((args.output_dir / "run_manifest.json").read_text())
        self.assertEqual(record["status"], "complete")
        self.assertEqual(record["eligible_genes"], 1)
        self.assertTrue(record["complete_counts_declared"])
        self.assertEqual(len(record["inputs"]), 3)
        self.assertEqual(
            (args.output_dir / "samples.tsv").read_bytes(), self.manifest.read_bytes()
        )
        with self.assertRaisesRegex(ValueError, "already exists"):
            prepare(args)

    def test_one_factor_grid_has_seven_unique_scenarios_and_rejects_bad_settings(self):
        args = self.sensitivity_options()
        scenarios = shared_scenarios(args)
        self.assertEqual(len(scenarios), 7)
        self.assertEqual(scenarios[0], ("primary", (5, 4, 50)))
        self.assertEqual(len(set(setting for _, setting in scenarios)), 7)
        self.assertTrue(
            all(
                sum(a != b for a, b in zip(setting, (5, 4, 50))) == 1
                for _, setting in scenarios[1:]
            )
        )
        args.minimum_junction_counts = "1.5"
        with self.assertRaisesRegex(ValueError, "integers"):
            shared_scenarios(args)

    def test_rust_handoff_preserves_zero_edge_topology_and_excludes_ineligible_tests(
        self,
    ):
        args = self.sensitivity_options()
        run(args)
        _, summary = read_table(
            args.output_dir / "sensitivity.tsv", ["scenario", "gene_id", "status"]
        )
        self.assertEqual(len(summary), 7 * 4)
        self.assertEqual(
            sum(
                row["scenario"] == "primary" and row["status"] == "eligible"
                for row in summary
            ),
            1,
        )
        self.assertTrue(
            all(
                not row["mean_difference_case_minus_control"]
                for row in summary
                if row["status"] != "eligible"
            )
        )
        _, metrics = read_table(
            args.output_dir / "primary/results/sample_metrics.tsv",
            ["gene_id", "sample_id", "cycle_fraction"],
        )
        control = next(
            row
            for row in metrics
            if row["gene_id"] == "TRIANGLE" and row["sample_id"] == "S7"
        )
        # On the common triangle, signal (1/2, 0, 1/2) projects to cycle energy 1/3.
        # Total energy is 1/2, hence fraction 2/3, despite that sample's zero shortcut.
        self.assertAlmostEqual(float(control["cycle_fraction"]), 2 / 3, places=8)
        _, tests = read_table(
            args.output_dir / "primary/statistics/statistics.tsv", ["gene_id", "status"]
        )
        self.assertEqual(
            [(row["gene_id"], row["status"]) for row in tests], [("TRIANGLE", "tested")]
        )
        _, family = read_table(
            args.output_dir / "primary/eligible_genes.tsv", ["gene_id", "gene_symbol"]
        )
        self.assertEqual(family[0]["gene_symbol"], "TRIANGLE_symbol")
        record = json.loads((args.output_dir / "run_manifest.json").read_text())
        self.assertEqual(record["status"], "complete")
        self.assertTrue(
            any(
                row["statistics_status"] == "no_eligible_genes"
                for row in record["scenarios"]
            )
        )
        with self.assertRaisesRegex(ValueError, "already exists"):
            run(args)

    def test_sensitivity_failure_is_recorded_without_complete_status(self):
        args = self.sensitivity_options()
        args.binary = Path("/usr/bin/false")
        with self.assertRaisesRegex(ValueError, "analysis failed"):
            run(args)
        record = json.loads((args.output_dir / "run_manifest.json").read_text())
        self.assertEqual(record["status"], "failed")
        self.assertNotEqual(record["scenarios"][0].get("status"), "complete")

    def test_native_cli_matches_python_scenarios_and_records_valid_checksums(self):
        args = self.sensitivity_options()
        run(args)
        native = self.root / "native"
        process = subprocess.run(
            [
                str(args.binary),
                "sensitivity",
                "--junctions",
                str(args.junctions),
                "--manifest",
                str(args.manifest),
                "--genes",
                str(args.genes),
                "--complete-counts",
                "--primary-edge-count",
                "5",
                "--primary-sample-count",
                "4",
                "--primary-gene-support",
                "50",
                "--minimum-junction-counts",
                "3,5,10",
                "--minimum-samples-per-edge",
                "2,4,6",
                "--minimum-gene-support",
                "20,50,100",
                "--run-stats",
                "--family",
                args.family,
                "--alternative",
                "two-sided",
                "--design",
                "independent",
                "--min-per-group",
                "4",
                "--output",
                str(native),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(
            (native / "sensitivity.tsv").read_text(),
            (args.output_dir / "sensitivity.tsv").read_text(),
        )
        record = json.loads((native / "run_manifest.json").read_text())
        self.assertEqual(record["status"], "complete")
        self.assertEqual(record["implementation"], "native Rust")
        self.assertEqual(
            record["tool"]["executable_sha256"],
            hashlib.sha256(args.binary.read_bytes()).hexdigest(),
        )
        for scenario in record["scenarios"]:
            name = scenario["scenario"]
            for file in (
                "junctions.tsv",
                "samples.tsv",
                "edges.tsv",
                "sample_eligibility.tsv",
                "gene_eligibility.tsv",
                "eligible_genes.tsv",
            ):
                self.assertEqual(
                    (native / name / file).read_text(),
                    (args.output_dir / name / file).read_text(),
                )
            if scenario["statistics_status"] == "complete":
                # Run IDs differ; every scientific value and supplied annotation must agree.
                for file in ("statistics.tsv", "set_test.tsv", "set_scores.tsv"):
                    _, expected = read_table(
                        args.output_dir / name / "statistics" / file, ["run_id"]
                    )
                    _, actual = read_table(
                        native / name / "statistics" / file, ["run_id"]
                    )
                    for rows in (expected, actual):
                        for row in rows:
                            row.pop("run_id")
                    self.assertEqual(actual, expected)
        for manifest in [native / "run_manifest.json"] + list(
            native.glob("*/run_manifest.json")
        ):
            metadata = json.loads(manifest.read_text())
            for output in metadata["outputs"]:
                self.assertEqual(
                    output["sha256"],
                    hashlib.sha256(
                        (manifest.parent / output["path"]).read_bytes()
                    ).hexdigest(),
                )


if __name__ == "__main__":
    unittest.main()
