from __future__ import annotations

import argparse
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from assemble_runs import assemble, main  # noqa: E402
from hisat2_pipeline import checksum  # noqa: E402

COUNT_UNIT = "primary uniquely mapped read alignments with NH:i:1 and XS strand"
PASS = {"read1_antisense": 990, "read1_sense": 10, "read2_sense": 985, "read2_antisense": 15}
HOLD = {"read1_antisense": 600, "read1_sense": 400}


class Cohort:
    """Eight fake finished runs sharing one fake reference."""

    def __init__(self, root: Path):
        self.root = root
        self.reference = root / "reference"
        self.reference.mkdir()
        for name in ["map.tsv", "splice_sites.tsv"] + [f"genome.{n}.ht2" for n in range(1, 9)]:
            (self.reference / name).write_text(name)
        self.samples = [(f"S{n}", "case" if n < 4 else "control") for n in range(8)]
        for sample_id, _ in self.samples:
            self.make_run(sample_id)

    def make_run(self, sample_id, strand_check=None, **changes):
        folder = self.root / "runs" / sample_id
        folder.mkdir(parents=True, exist_ok=True)
        counts = folder / "junctions.tsv"
        counts.write_text(
            "sample_id\tgene_id\tfrom\tto\tcount\n"
            f"{sample_id}\tG1\tA\tB\t7\n{sample_id}\tG1\tB\tC\t3\n"
        )
        reads = folder / "reads.fq"
        reads.write_text(sample_id)
        ref = self.reference
        run = {
            "pipeline": "Splice Girl Cleaver/HISAT2",
            "status": "complete",
            "settings": {
                "sample_id": sample_id,
                "strandedness": "RF",
                "minimum_anchor": 8,
                "minimum_unique_reads": 1,
                "exclude_duplicates": False,
                "count_unit": COUNT_UNIT,
            },
            "inputs": [checksum(reads), checksum(ref / "map.tsv"), checksum(ref / "splice_sites.tsv")]
            + [checksum(ref / f"genome.{n}.ht2") for n in range(1, 9)],
            "outputs": [checksum(counts)],
            "commands": [
                ["hisat2", "-x", str(ref / "genome"), "--known-splicesite-infile",
                 str((ref / "splice_sites.tsv").resolve())],
                ["python", "sam_junctions.py", "--gene-map", str((ref / "map.tsv").resolve())],
            ],
            "stages": [{"status": "complete"}],
        }
        if strand_check is not None:
            run["strand_check"] = strand_check
        for key, value in changes.items():
            if key in run["settings"]:
                run["settings"][key] = value
            else:
                run[key] = value
        (folder / "run_manifest.json").write_text(json.dumps(run))
        return folder

    def strand_file(self, sample_id, counts, expected="RF"):
        path = self.root / "strand" / f"{sample_id}.json"
        path.parent.mkdir(exist_ok=True)
        path.write_text(json.dumps({**counts, "expected": expected, "sample_id": sample_id}))
        return path

    def sheet(self, strand=None):
        strand = strand if strand is not None else {
            s: self.strand_file(s, PASS) for s, _ in self.samples
        }
        lines = ["sample_id\tcondition\trun_dir\tstrand_check"]
        for sample_id, condition in self.samples:
            path = strand.get(sample_id)
            lines.append(
                f"{sample_id}\t{condition}\truns/{sample_id}\t"
                f"{path.relative_to(self.root) if path else ''}"
            )
        sheet = self.root / "runs.tsv"
        sheet.write_text("\n".join(lines) + "\n")
        return sheet

    def args(self, sheet, **overrides):
        values = dict(
            samples=sheet,
            output_dir=self.root / "cohort",
            strandedness="RF",
            min_anchor=8,
            min_unique=1,
            exclude_duplicates=False,
        )
        values.update(overrides)
        return argparse.Namespace(**values)


class AssembleTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.cohort = Cohort(Path(self.temporary.name))

    def tearDown(self):
        self.temporary.cleanup()

    def test_eight_runs_become_one_batch_record(self):
        record = assemble(self.cohort.args(self.cohort.sheet()))
        output = self.cohort.root / "cohort"
        saved = json.loads((output / "batch_manifest.json").read_text())
        self.assertEqual(saved, record)
        self.assertEqual(saved["status"], "complete")
        self.assertEqual(saved["pipeline"], "Splice Girl assembled Cleaver/HISAT2 runs")
        self.assertEqual(len(saved["samples"]), 8)
        self.assertEqual(saved["excluded_samples"], [])
        first = saved["samples"][0]
        self.assertEqual(first["strand_check"]["status"], "pass")
        self.assertEqual(first["strand_check"]["source"], "strand_orientation.py")
        self.assertTrue(first["run_manifest"]["path"].endswith("S0/run_manifest.json"))
        table = (output / "junctions.tsv").read_text().splitlines()
        self.assertEqual(len(table), 1 + 16)
        self.assertEqual(
            (output / "samples.tsv").read_text().splitlines()[1:3],
            ["S0\tcase", "S1\tcase"],
        )
        outputs = {Path(item["path"]).name: item["sha256"] for item in saved["outputs"]}
        self.assertEqual(outputs["junctions.tsv"], checksum(output / "junctions.tsv")["sha256"])

    def test_strand_check_recorded_in_run_manifest_needs_no_file(self):
        for sample_id, _ in self.cohort.samples:
            self.cohort.make_run(sample_id, strand_check={**PASS, "status": "pass"})
        record = assemble(self.cohort.args(self.cohort.sheet(strand={})))
        self.assertEqual(record["samples"][0]["strand_check"]["source"], "run_manifest")

    def test_held_sample_is_excluded_and_status_is_recomputed(self):
        strand = {s: self.cohort.strand_file(s, PASS) for s, _ in self.cohort.samples}
        # A stored "pass" is ignored: the call is recomputed from the counts.
        strand["S2"] = self.cohort.strand_file("S2", {**HOLD, "status": "pass"})
        record = assemble(self.cohort.args(self.cohort.sheet(strand)))
        self.assertEqual([s["sample_id"] for s in record["excluded_samples"]], ["S2"])
        self.assertEqual(
            record["excluded_samples"][0]["exclusion_reason"],
            "strand_check_hold_for_review",
        )
        output = self.cohort.root / "cohort"
        self.assertNotIn("S2\t", (output / "junctions.tsv").read_text())
        self.assertNotIn("S2\t", (output / "samples.tsv").read_text())
        self.assertEqual(len(record["samples"]), 7)

    def test_missing_strand_check_or_wrong_settings_stop_before_writing(self):
        strand = {s: self.cohort.strand_file(s, PASS) for s, _ in self.cohort.samples}
        no_check = dict(strand)
        del no_check["S5"]
        cases = [
            (no_check, "S5: no strand check", {}),
            (strand, "minimum_anchor is 8, expected 10", {"min_anchor": 10}),
        ]
        for checks, message, overrides in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(ValueError, message):
                    assemble(self.cohort.args(self.cohort.sheet(checks), **overrides))
                self.assertFalse((self.cohort.root / "cohort").exists())

    def test_runs_must_share_one_reference_and_be_unchanged(self):
        sheet = self.cohort.sheet()
        other = self.cohort.root / "reference" / "map.tsv"
        folder = self.cohort.root / "runs" / "S3"
        run = json.loads((folder / "run_manifest.json").read_text())
        for item in run["inputs"]:
            if item["path"] == str(other.resolve()):
                item["sha256"] = "0" * 64
        (folder / "run_manifest.json").write_text(json.dumps(run))
        with self.assertRaisesRegex(ValueError, "S3: gene_map differ"):
            assemble(self.cohort.args(sheet))
        self.cohort.make_run("S3")
        with (folder / "junctions.tsv").open("a") as handle:
            handle.write("S3\tG1\tC\tD\t1\n")
        with self.assertRaisesRegex(ValueError, "S3: junctions.tsv changed"):
            assemble(self.cohort.args(sheet))

    def test_incomplete_run_and_wrong_strand_file_are_refused(self):
        self.cohort.make_run("S1", stages=[{"status": "complete"}, {"status": "failed"}])
        with self.assertRaisesRegex(ValueError, "S1: run is not a completed"):
            assemble(self.cohort.args(self.cohort.sheet()))
        self.cohort.make_run("S1")
        strand = {s: self.cohort.strand_file(s, PASS) for s, _ in self.cohort.samples}
        strand["S4"] = strand["S6"]
        with self.assertRaisesRegex(ValueError, "S4: strand_check file is for S6"):
            assemble(self.cohort.args(self.cohort.sheet(strand)))

    def test_command_line_reports_exclusions_and_refuses_existing_output(self):
        strand = {s: self.cohort.strand_file(s, PASS) for s, _ in self.cohort.samples}
        strand["S7"] = self.cohort.strand_file("S7", HOLD)
        sheet = self.cohort.sheet(strand)
        argv = ["--samples", str(sheet), "--output-dir", str(self.cohort.root / "cohort"),
                "--strandedness", "RF", "--min-anchor", "8", "--min-unique", "1"]
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(main(argv), 0)
        self.assertIn("Assembled 7 samples", out.getvalue())
        self.assertIn("S7 excluded by the strand check (read 1 agreement 60.00%)", out.getvalue())
        with contextlib.redirect_stderr(io.StringIO()) as err:
            self.assertEqual(main(argv), 1)
        self.assertIn("already exists", err.getvalue())


if __name__ == "__main__":
    unittest.main()
