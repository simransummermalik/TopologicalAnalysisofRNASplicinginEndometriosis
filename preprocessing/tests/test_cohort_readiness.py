import argparse
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cohort_readiness import pair_counts, run


class CohortReadinessTests(unittest.TestCase):
    def test_mate_counts_must_agree_and_millions_are_converted_to_pairs(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "qc.tsv"
            path.write_text(
                "Sample\tfastqc-total_sequences\nS_forward\t1.234567\nS_reverse\t1.234567\n"
            )
            self.assertEqual(pair_counts([path]), {"S": 1234567})
            path.write_text(
                path.read_text().replace("S_reverse\t1.234567", "S_reverse\t1.2")
            )
            with self.assertRaisesRegex(ValueError, "unequal"):
                pair_counts([path])

    def test_missing_inputs_never_become_ready_and_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, qc, subset, sam, bam = [
                root / name
                for name in [
                    "metadata.tsv",
                    "qc.tsv",
                    "subset.json",
                    "pilot.sam",
                    "pilot.bam",
                ]
            ]
            metadata.write_text("sample_id\tcondition\nS\tcontrol\n")
            qc.write_text(
                "Sample\tfastqc-total_sequences\nS_forward\t0.000002\nS_reverse\t0.000002\n"
            )
            subset.write_text(
                json.dumps(
                    {"written_pairs": 1, "outputs": [{"bytes": 10}, {"bytes": 10}]}
                )
            )
            sam.write_bytes(b"x" * 100)
            bam.write_bytes(b"x" * 10)
            args = argparse.Namespace(
                metadata=metadata,
                qc=[qc],
                samples=None,
                storage=root,
                pilot_manifest=subset,
                pilot_sam=sam,
                pilot_bam=bam,
                output=root / "check.json",
            )
            with contextlib.redirect_stdout(io.StringIO()):
                run(args)
            result = json.loads(args.output.read_text())
            self.assertEqual(
                result["estimated_working_bytes_with_25_percent_margin"], 375
            )
            self.assertFalse(result["all_paired_inputs_present"])
            with self.assertRaisesRegex(ValueError, "exists"):
                run(args)
            sheet = root / "samples.tsv"
            sheet.write_text("sample_id\tcondition\tread1\tread2\nS\tcase\tR1\tR2\n")
            args.samples = sheet
            with self.assertRaisesRegex(ValueError, "condition disagrees"):
                run(args)


if __name__ == "__main__":
    unittest.main()
