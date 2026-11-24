from __future__ import annotations

import json
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from build_junction_gene_map import build_gene_map  # noqa: E402
from hisat2_pipeline import strand_check  # noqa: E402
from sam_junctions import count_junctions  # noqa: E402
from strand_orientation import count_sam, summarize  # noqa: E402

FIXTURES = Path(__file__).parent / "fixtures"
READ1, READ2, PAIRED, REVERSE = 64, 128, 1, 16


def sam(flag, xs="+", cigar="50M100N50M", nh="NH:i:1"):
    tags = "\t".join(tag for tag in (nh, f"XS:A:{xs}" if xs else "") if tag)
    return f"r\t{flag}\tchr1\t150\t60\t{cigar}\t*\t0\t0\t*\t*\t{tags}\n"


def rf_pair(xs="+"):
    """An RF fragment: read 1 on the opposite strand to the transcript."""
    if xs == "+":
        return [sam(PAIRED | READ1 | REVERSE, xs), sam(PAIRED | READ2, xs)]
    return [sam(PAIRED | READ1, xs), sam(PAIRED | READ2 | REVERSE, xs)]


class StrandOrientationTests(unittest.TestCase):
    def test_rf_pairs_on_both_transcript_strands(self):
        counts = count_sam(rf_pair("+") + rf_pair("-"))
        summary = summarize(counts, "RF")
        self.assertEqual(summary["read1_antisense"], 2)
        self.assertEqual(summary["read2_sense"], 2)
        self.assertEqual(summary["read1_antisense_fraction"], 1.0)
        self.assertEqual(summary["status"], "pass")
        self.assertEqual(summarize(counts, "FR")["status"], "hold_for_review")

    def test_ineligible_records_are_not_counted(self):
        lines = [
            sam(PAIRED | READ1, cigar="100M"),  # unspliced
            sam(PAIRED | READ1, nh="NH:i:2"),  # multimapped
            sam(PAIRED | READ1, xs=None),  # no transcript strand
            sam(PAIRED | READ1 | 256),  # secondary
            sam(PAIRED | READ1 | 4),  # unmapped
        ]
        summary = summarize(count_sam(lines), "RF")
        self.assertIsNone(summary["read1_antisense_fraction"])
        self.assertEqual(summary["status"], "no_read1_records")

    def test_threshold_is_95_percent_of_read1(self):
        counts = Counter(read1_antisense=95, read1_sense=5)
        self.assertEqual(summarize(counts, "RF")["status"], "pass")
        counts = Counter(read1_antisense=94, read1_sense=6)
        self.assertEqual(summarize(counts, "RF")["status"], "hold_for_review")
        self.assertNotIn("status", summarize(counts))

    def test_junction_counter_and_pipeline_record_the_check(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            gene_map = root / "map.tsv"
            build_gene_map(FIXTURES / "annotation.gtf", gene_map)
            lines = rf_pair("+") * 19 + [sam(PAIRED | READ1, "+")]
            _, stats = count_junctions(lines, gene_map, "sample", 8, 1)
            orientation = stats["strand_orientation"]
            self.assertEqual(orientation["read1_antisense"], 19)
            self.assertEqual(orientation["read1_sense"], 1)
            qc = root / "junction_qc.json"
            qc.write_text(json.dumps(stats), encoding="utf-8")
            check = strand_check(qc, "RF")
            self.assertEqual(check["status"], "pass")
            self.assertAlmostEqual(check["read1_agreement"], 0.95)
            self.assertEqual(strand_check(qc, "unstranded")["status"], "not_checked")

    def test_missing_or_empty_paired_orientation_counts_do_not_pass(self):
        with tempfile.TemporaryDirectory() as folder:
            qc = Path(folder) / "junction_qc.json"
            for record in ({}, {"strand_orientation": {}}):
                qc.write_text(json.dumps(record), encoding="utf-8")
                self.assertEqual(strand_check(qc, "RF")["status"], "no_read1_records")
                self.assertEqual(strand_check(qc, "FR")["status"], "no_read1_records")

    def test_invalid_agreement_threshold_is_rejected(self):
        for minimum in (-0.1, 1.1, float("nan"), float("inf")):
            with self.subTest(minimum=minimum), self.assertRaisesRegex(ValueError, "between 0 and 1"):
                summarize(Counter(read1_antisense=100), "RF", minimum)


if __name__ == "__main__":
    unittest.main()
