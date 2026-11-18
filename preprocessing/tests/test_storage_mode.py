"""Storage/process regressions use tiny inputs; real formats are checked by smoke()."""

import argparse
import contextlib
import io
import json
import os
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from batch_preprocessing import batch
from build_junction_gene_map import build_gene_map
from hisat2_pipeline import command_groups, configuration, run
from pipeline_runtime import disk_budget, run_stage, storage_settings
from sam_junctions import convert_lines, convert_sam

HERE = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).parent / "fixtures"
SAM = "read\t0\tchr1\t150\t60\t50M100N50M\t*\t0\t0\t*\t*\tNH:i:1\tXS:A:+\n"


class StorageTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.read = self.root / "reads.fastq"
        self.read.write_bytes(b"x" * 100)
        self.map = self.root / "map.tsv"
        build_gene_map(FIXTURES / "annotation.gtf", self.map)
        for number in range(1, 9):
            (self.root / f"genome.{number}.ht2").write_text("fixture")
        self.args = argparse.Namespace(
            command="run",
            sample_id="S",
            read1=self.read,
            read2=None,
            index=self.root / "genome",
            gene_map=self.map,
            splice_sites=self.map,
            output_dir=self.root / "output",
            strandedness="unstranded",
            threads=2,
            min_anchor=8,
            min_unique=1,
            skip_cleaning=True,
            exclude_duplicates=False,
            storage_mode="stream-bam",
            sort_memory_mb=256,
            min_free_gb=0,
        )

    def tearDown(self):
        self.temporary.cleanup()

    def test_stream_plan_has_real_pipe_groups_no_sam_file_and_explicit_sort_limits(
        self,
    ):
        output, commands = configuration(self.args)
        groups = dict(command_groups(self.args, commands))
        self.assertEqual([len(group) for group in groups.values()], [2, 1, 2, 1])
        align, sort = groups["align_sort"]
        self.assertNotIn("-S", align)
        self.assertIn("256M", sort)
        self.assertEqual(groups["junctions"][1][3:5], ["-", "--gene-map"])
        self.assertFalse(
            any(str(output / "aligned.sam") in command for command in commands)
        )
        self.assertFalse(output.exists())
        self.args.storage_mode = "retain-sam"
        _, legacy = configuration(self.args)
        self.assertIn(str(output / "aligned.sam"), legacy[0])

    def test_disk_budget_counts_only_new_work_and_reserve(self):
        streaming = disk_budget(self.args, self.args.output_dir)
        self.assertEqual(streaming["estimated_additional_bytes"], 300)
        self.args.skip_cleaning = False
        self.assertEqual(
            disk_budget(self.args, self.args.output_dir)["estimated_additional_bytes"],
            400,
        )
        self.args.storage_mode = "retain-sam"
        self.assertEqual(
            disk_budget(self.args, self.args.output_dir)["estimated_additional_bytes"],
            1300,
        )
        for value in [float("nan"), float("inf"), -1]:
            self.args.min_free_gb = value
            with self.assertRaisesRegex(ValueError, "min-free-gb"):
                storage_settings(self.args)

    def test_insufficient_disk_refuses_before_tools_or_output_creation(self):
        original = self.read.read_bytes()
        with (
            patch(
                "pipeline_runtime.shutil.disk_usage",
                return_value=SimpleNamespace(free=1),
            ),
            patch("hisat2_pipeline.require_tools") as tools,
        ):
            with self.assertRaisesRegex(ValueError, "No processing started"):
                run(self.args)
            tools.assert_not_called()
        self.assertFalse(self.args.output_dir.exists())
        self.assertEqual(self.read.read_bytes(), original)

    def test_batch_budget_includes_all_samples_before_starting_any(self):
        read2 = self.root / "other.fastq"
        read2.write_bytes(b"x" * 100)
        sheet = self.root / "samples.tsv"
        sheet.write_text(
            "sample_id\tcondition\tread1\tread2\nA\tcase\treads.fastq\t\nB\tcontrol\tother.fastq\t\n"
        )
        self.args.samples = sheet
        with (
            patch(
                "pipeline_runtime.shutil.disk_usage",
                return_value=SimpleNamespace(free=400),
            ),
            patch("batch_preprocessing.run") as execute,
        ):
            with self.assertRaisesRegex(ValueError, "not enough disk"):
                batch(self.args)
            execute.assert_not_called()
        self.assertFalse(self.args.output_dir.exists())

    def test_file_and_stream_converter_have_identical_counts_qc_and_protection(self):
        path = self.root / "reads.sam"
        path.write_text(SAM + SAM.replace("NH:i:1", "NH:i:2"))
        convert_sam(
            path, self.map, "S", self.root / "file.tsv", self.root / "file.json", 8, 1
        )
        convert_lines(
            io.StringIO(path.read_text()),
            self.map,
            "S",
            self.root / "stream.tsv",
            self.root / "stream.json",
            8,
            1,
        )
        self.assertEqual(
            (self.root / "file.tsv").read_bytes(),
            (self.root / "stream.tsv").read_bytes(),
        )
        self.assertEqual(
            (self.root / "file.json").read_bytes(),
            (self.root / "stream.json").read_bytes(),
        )
        with self.assertRaisesRegex(ValueError, "overwrite"):
            convert_lines(
                io.StringIO(SAM),
                self.map,
                "S",
                self.root / "stream.tsv",
                self.root / "stream.json",
                8,
                1,
            )

    def test_stage_checks_producer_exit_even_when_consumer_succeeds(self):
        stages = []
        producer = [
            sys.executable,
            "-c",
            "import sys; print('partial data'); sys.exit(7)",
        ]
        consumer = [sys.executable, "-c", "import sys; sys.stdin.read()"]
        with (
            contextlib.redirect_stdout(io.StringIO()),
            self.assertRaisesRegex(ValueError, "exit codes"),
        ):
            run_stage("pipe", [producer, consumer], self.root, 0, stages, lambda: None)
        self.assertEqual(stages[0]["status"], "failed")
        self.assertEqual(stages[0]["return_codes"][0], 7)
        self.assertTrue(stages[0]["logs"])

    def test_consumer_failure_stops_the_producer_without_hanging(self):
        stages = []
        producer = [
            sys.executable,
            "-c",
            "import os\nwhile True: os.write(1, b'x' * 65536)",
        ]
        consumer = [sys.executable, "-c", "import sys; sys.exit(9)"]
        started = time.monotonic()
        with (
            contextlib.redirect_stdout(io.StringIO()),
            self.assertRaisesRegex(ValueError, "stage pipe failed"),
        ):
            run_stage("pipe", [producer, consumer], self.root, 0, stages, lambda: None)
        self.assertLess(time.monotonic() - started, 10)
        self.assertEqual(stages[0]["return_codes"][1], 9)
        self.assertTrue(all(code is not None for code in stages[0]["return_codes"]))

    def test_keyboard_interrupt_stops_child_and_records_interruption(self):
        stages = []
        interrupted = False

        def persist():
            nonlocal interrupted
            if stages[0]["pids"] and not interrupted:
                interrupted = True
                raise KeyboardInterrupt()

        with (
            contextlib.redirect_stdout(io.StringIO()),
            self.assertRaises(KeyboardInterrupt),
        ):
            run_stage(
                "slow",
                [[sys.executable, "-c", "import time; time.sleep(30)"]],
                self.root,
                0,
                stages,
                persist,
            )
        self.assertEqual(stages[0]["status"], "interrupted")
        self.assertIsNotNone(stages[0]["return_codes"][0])

    def test_disk_reserve_is_rechecked_before_starting_a_stage(self):
        stages = []
        with (
            patch(
                "pipeline_runtime.shutil.disk_usage",
                return_value=SimpleNamespace(free=1),
            ),
            patch("pipeline_runtime.subprocess.Popen") as process,
        ):
            with self.assertRaisesRegex(ValueError, "reserve reached"):
                run_stage(
                    "next", [["not-started"]], self.root, 100, stages, lambda: None
                )
            process.assert_not_called()
        self.assertEqual(stages[0]["status"], "failed")

    def fake_commands(self, failed=None):
        """Exercise orchestration with fake format bytes, not BAM-format validation."""
        output = self.args.output_dir
        partial = output / "aligned.partial.bam"
        producer = [
            sys.executable,
            "-c",
            f"import sys; sys.stdout.write({SAM!r}); sys.exit({7 if failed == 'align' else 0})",
        ]
        sorter = [
            sys.executable,
            "-c",
            f"import sys; from pathlib import Path; Path({str(partial)!r}).write_text(sys.stdin.read())",
        ]
        quickcheck = [sys.executable, "-c", "pass"]
        view = [
            sys.executable,
            "-c",
            f"import sys; from pathlib import Path; sys.stdout.write(Path({str(partial)!r}).read_text()); sys.exit({8 if failed == 'view' else 0})",
        ]
        converter = [
            sys.executable,
            str(HERE / "sam_junctions.py"),
            "--sam",
            "-",
            "--gene-map",
            str(self.map),
            "--sample-id",
            "S",
            "--output",
            str(output / "junctions.partial.tsv"),
            "--qc-json",
            str(output / "junction_qc.json"),
            "--min-anchor",
            "8",
            "--min-unique",
            "1",
        ]
        indexer = [
            sys.executable,
            "-c",
            f"import sys; from pathlib import Path; Path({str(partial) + '.bai'!r}).write_bytes(b'fixture'); sys.exit({9 if failed == 'index' else 0})",
        ]
        return output, [producer, sorter, quickcheck, view, converter, indexer]

    def test_partial_bam_and_counts_are_not_published_after_any_pipeline_failure(self):
        for failed in ["align", "view", "index"]:
            with self.subTest(failed=failed):
                self.args.output_dir = self.root / failed
                with (
                    patch(
                        "hisat2_pipeline.configuration",
                        return_value=self.fake_commands(failed),
                    ),
                    patch("hisat2_pipeline.require_tools"),
                    patch("hisat2_pipeline.versions", return_value={}),
                    contextlib.redirect_stdout(io.StringIO()),
                ):
                    with self.assertRaisesRegex(ValueError, "failed"):
                        run(self.args)
                output = self.args.output_dir
                record = json.loads((output / "run_manifest.json").read_text())
                self.assertEqual(record["status"], "failed")
                self.assertEqual(record["stages"][-1]["status"], "failed")
                self.assertFalse((output / "aligned.bam").exists())
                self.assertFalse((output / "junctions.tsv").exists())
                self.assertTrue((output / "aligned.partial.bam").is_file())
                self.assertEqual(self.read.read_bytes(), b"x" * 100)

    @unittest.skipUnless(os.name == "posix", "POSIX signal/process-group regression")
    def test_sigterm_records_interrupted_run_and_stops_streaming_children(self):
        # A separate Python parent receives the real signal, never this test runner.
        output, commands = self.fake_commands()
        commands[:2] = [
            [sys.executable, "-c", "import time; time.sleep(30)"] for _ in range(2)
        ]
        script = self.root / "interrupt_case.py"
        values = vars(self.args)
        script.write_text(
            "import argparse,sys\nfrom pathlib import Path\n"
            f"sys.path.insert(0,{str(HERE)!r})\n"
            "import hisat2_pipeline as pipeline\nfrom pipeline_runtime import interruption_handlers,RunInterrupted\n"
            f"args=argparse.Namespace(**{values!r})\n"
            f"pipeline.configuration=lambda args:(Path({str(output)!r}),{commands!r})\n"
            "pipeline.require_tools=lambda names:None\npipeline.versions=lambda names:{}\n"
            "try:\n with interruption_handlers(): pipeline.run(args)\n"
            "except RunInterrupted as error: sys.exit(128+error.signum)\n"
        )
        # Path reprs are PosixPath(...); bind that name in the child namespace.
        script.write_text(
            script.read_text().replace(
                "from pathlib import Path", "from pathlib import Path,PosixPath"
            )
        )
        process = subprocess.Popen(
            [sys.executable, str(script)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        pids = []
        try:
            manifest = output / "run_manifest.json"
            deadline = time.monotonic() + 8
            while time.monotonic() < deadline:
                if manifest.exists():
                    record = json.loads(manifest.read_text())
                    if record["stages"] and len(record["stages"][-1]["pids"]) == 2:
                        pids = record["stages"][-1]["pids"]
                        break
                if process.poll() is not None:
                    self.fail(str(process.communicate()))
                time.sleep(0.02)
            self.assertEqual(len(pids), 2)
            process.send_signal(signal.SIGTERM)
            stdout, stderr = process.communicate(timeout=8)
            self.assertEqual(process.returncode, 143, stdout + stderr)
            record = json.loads(manifest.read_text())
            self.assertEqual(record["status"], "interrupted")
            self.assertEqual(record["stages"][-1]["status"], "interrupted")
            self.assertFalse((output / "junctions.tsv").exists())
            for pid in pids:
                with self.assertRaises(ProcessLookupError):
                    os.kill(pid, 0)
        finally:
            if process.poll() is None:
                process.kill()
            process.communicate()
            for pid in pids:
                try:
                    os.killpg(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass


if __name__ == "__main__":
    unittest.main()
