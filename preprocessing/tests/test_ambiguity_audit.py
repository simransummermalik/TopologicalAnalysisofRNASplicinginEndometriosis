import argparse
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ambiguity_audit import classify, observe, run


class AmbiguityAuditTests(unittest.TestCase):
    def test_distinguishes_exon_pair_and_gene_ambiguity_without_rescuing_counts(self):
        self.assertEqual(classify(set()), "unannotated")
        self.assertEqual(classify({("G", "A", "B")}), "assigned")
        self.assertEqual(
            classify({("G", "A", "B"), ("G", "C", "B")}), "within_gene_exon_pairs"
        )
        self.assertEqual(classify({("G", "A", "B"), ("H", "A", "B")}), "cross_gene")
        self.assertEqual(
            classify({("G", "A", "B"), ("G", "C", "B"), ("H", "A", "B")}),
            "cross_gene_and_exon_pairs",
        )

    def test_multijunction_reads_count_strandedness_once_and_apply_primary_filters(
        self,
    ):
        def sam(flag=80, tags="NH:i:1\tXS:A:+", cigar="10M20N10M20N10M"):
            return f"read\t{flag}\tchr1\t100\t60\t{cigar}\t*\t0\t0\t*\t*\t{tags}\n"

        mapping = {
            ("chr1", 110, 129, "1"): {("G", "A", "B")},
            ("chr1", 140, 159, "1"): {("G", "A", "B"), ("G", "C", "B")},
        }
        observed, qc, strand = observe(
            [
                sam(),
                sam(flag=256),
                sam(tags="NH:i:2\tXS:A:+"),
                sam(tags="NH:i:1"),
                sam(cigar="5M20N10M"),
            ],
            mapping,
            8,
        )
        self.assertEqual(sum(observed.values()), 2)
        self.assertEqual(qc["assigned_junction_observations"], 1)
        self.assertEqual(qc["ambiguous_junction_observations"], 1)
        self.assertEqual(strand["read1_antisense"], 1)
        self.assertEqual(qc["insufficient_anchor_records"], 1)

    def test_saved_summary_matches_expected_qc_and_preserves_cross_gene_potential(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            gene_map, sam, candidates, expected = [
                root / name
                for name in ["map.tsv", "reads.sam", "candidates.tsv", "qc.json"]
            ]
            gene_map.write_text(
                "chromosome\tintron_start\tintron_end\tstrand\tgene_id\tfrom_node\tto_node\nchr1\t110\t129\t1\tG\tA\tB\nchr1\t110\t129\t1\tH\tA\tB\n"
            )
            sam.write_text(
                "read\t80\tchr1\t100\t60\t10M20N10M\t*\t0\t0\t*\t*\tNH:i:1\tXS:A:+\n"
            )
            candidates.write_text(
                "gene_id\tgene_symbol\tcategory\nG\tGeneG\trisk_only\nH\tGeneH\tenvironment_only\nNONE\tMissing\toverlap\n"
            )
            expected.write_text(
                json.dumps(
                    dict(
                        sample_id="S",
                        minimum_anchor=8,
                        exclude_duplicates=False,
                        sam_records=1,
                        ambiguous_junction_observations=1,
                    )
                )
            )
            args = argparse.Namespace(
                sam=sam,
                bam=None,
                gene_map=gene_map,
                candidates=candidates,
                expected_qc=expected,
                sample_id="S",
                min_anchor=8,
                exclude_duplicates=False,
                output_dir=root / "audit",
            )
            with contextlib.redirect_stdout(io.StringIO()):
                run(args)
            summary = json.loads((args.output_dir / "summary.json").read_text())
            self.assertEqual(summary["source_qc_agreement"], "PASS")
            self.assertEqual(summary["observations_by_class"]["cross_gene"], 1)
            self.assertEqual(
                summary[
                    "same_gene_assignment_upper_bound_if_exon_pair_ambiguity_were_ignored"
                ],
                0,
            )
            self.assertIn(
                "NONE\tMissing\toverlap\tno_eligible_observation",
                (args.output_dir / "candidate_gene_coverage.tsv").read_text(),
            )
            with self.assertRaisesRegex(ValueError, "folder exists"):
                run(args)
