import argparse
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from batch_preprocessing import batch, load_samples


class BatchTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        for name in ["a.fastq", "b.fastq", "map.tsv", "hints.tsv"]:
            (self.root / name).write_text("fixture")
        for number in range(1, 9):
            (self.root / f"genome.{number}.ht2").write_text("fixture")
        self.sheet = self.root / "samples.tsv"
        self.sheet.write_text(
            "sample_id\tcondition\tread1\tread2\nA\tcase\ta.fastq\t\nB\tcontrol\tb.fastq\t\n"
        )
        self.args = argparse.Namespace(
            command="run",
            samples=self.sheet,
            index=self.root / "genome",
            gene_map=self.root / "map.tsv",
            splice_sites=self.root / "hints.tsv",
            output_dir=self.root / "output",
            strandedness="unstranded",
            min_anchor=8,
            min_unique=1,
            threads=2,
            skip_cleaning=False,
            exclude_duplicates=False,
        )

    def tearDown(self):
        self.temp.cleanup()

    @staticmethod
    def fake_run(job):
        job.output_dir.mkdir()
        (job.output_dir / "junctions.tsv").write_text(
            f"sample_id\tgene_id\tfrom\tto\tcount\n{job.sample_id}\tG\tE1\tE2\t1\n"
        )
        (job.output_dir / "run_manifest.json").write_text('{"status":"complete"}')

    def test_merge_keeps_samples_separate_and_protects_existing_outputs(self):
        with (
            patch("batch_preprocessing.require_tools"),
            patch("batch_preprocessing.run", side_effect=self.fake_run),
            contextlib.redirect_stdout(io.StringIO()),
        ):
            batch(self.args)
        output = self.args.output_dir
        self.assertEqual(
            json.loads((output / "batch_manifest.json").read_text())["status"],
            "complete",
        )
        self.assertEqual((output / "junctions.tsv").read_text().count("\tG\t"), 2)
        self.assertIn("B\tcontrol", (output / "samples.tsv").read_text())
        with self.assertRaisesRegex(ValueError, "already exists"):
            batch(self.args)

    def test_failed_sample_never_publishes_a_complete_merged_table(self):
        def fail_second(job):
            if job.sample_id == "B":
                raise ValueError("synthetic failure")
            self.fake_run(job)

        with (
            patch("batch_preprocessing.require_tools"),
            patch("batch_preprocessing.run", side_effect=fail_second),
            self.assertRaisesRegex(ValueError, "synthetic failure"),
        ):
            batch(self.args)
        output = self.args.output_dir
        manifest = json.loads((output / "batch_manifest.json").read_text())
        self.assertEqual(manifest["status"], "failed")
        self.assertEqual(
            [r["status"] for r in manifest["samples"]], ["complete", "failed"]
        )
        self.assertFalse((output / "junctions.tsv").exists())

    def test_all_paths_are_checked_before_running_and_reused_reads_are_rejected(self):
        self.sheet.write_text(
            self.sheet.read_text().replace("b.fastq", "missing.fastq")
        )
        with self.assertRaisesRegex(ValueError, "missing read1"):
            batch(self.args)
        self.assertFalse(self.args.output_dir.exists())
        self.sheet.write_text(
            self.sheet.read_text().replace("missing.fastq", "a.fastq")
        )
        with self.assertRaisesRegex(ValueError, "reused"):
            load_samples(self.sheet)

    def test_plan_creates_no_outputs(self):
        self.args.command = "plan"
        with contextlib.redirect_stdout(io.StringIO()) as captured:
            batch(self.args)
        self.assertIn("hisat2 --dta-cufflinks", captured.getvalue())
        self.assertFalse(self.args.output_dir.exists())


if __name__ == "__main__":
    unittest.main()
