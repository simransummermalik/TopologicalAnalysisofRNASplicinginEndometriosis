import tempfile
import unittest
from pathlib import Path

from validate_real_pilot import gene_targets, independent_counts, select_gene


class RealPilotTests(unittest.TestCase):
    def test_selection_uses_coverage_then_id_and_never_unlisted_genes(self):
        rows = []
        for gene in ["B", "A", "UNLISTED"]:
            for source, target in [("E1", "E2"), ("E2", "E3"), ("E1", "E3")]:
                rows.append(
                    dict(
                        sample_id="S",
                        gene_id=gene,
                        **{"from": source, "to": target},
                        count="1000" if gene == "UNLISTED" else "10",
                    )
                )
        (total, gene, selected), eligible = select_gene(
            rows, {"A": {}, "B": {}}, "S", 3, 3, 12, 20
        )
        self.assertEqual((total, gene, len(selected), eligible), (30, "A", 3, 2))
        with self.assertRaisesRegex(ValueError, "no candidate"):
            select_gene(rows, {"A": {}}, "S", 11, 3, 12, 20)
        with self.assertRaisesRegex(ValueError, "another sample"):
            select_gene(rows, {"A": {}}, "OTHER", 3, 3, 12, 20)

    def test_annotation_recount_preserves_cross_gene_ambiguity(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "map.tsv"
            path.write_text(
                "chromosome\tintron_start\tintron_end\tstrand\tgene_id\tfrom_node\tto_node\nchr1\t110\t129\t1\tG\tE1\tE2\nchr1\t110\t129\t1\tH\tE1\tE2\nchr1\t210\t229\t2\tG\tE3\tE2\n"
            )
            self.assertEqual(
                gene_targets(path, "G"), {("chr1", 210, 229, "2"): ("E3", "E2")}
            )

    def test_independent_recount_checks_coordinates_anchor_strand_and_uniqueness(self):
        def sam(flag=0, tags="NH:i:1\tXS:A:+", position=100, cigar="10M20N10M"):
            return f"read\t{flag}\tchr1\t{position}\t60\t{cigar}\t*\t0\t0\t*\t*\t{tags}"

        targets = {("chr1", 110, 129, "1"): ("E1", "E2")}
        records = [
            sam(flag=16),
            sam(cigar="5S10M20N10M"),
            sam(flag=2048),
            sam(tags="NH:i:2\tXS:A:+"),
            sam(tags="NH:i:1"),
            sam(tags="NH:i:1\tXS:A:-"),
            sam(position=105, cigar="5M20N10M"),
        ]
        self.assertEqual(independent_counts(records, targets, 8), {("E1", "E2"): 2})


if __name__ == "__main__":
    unittest.main()
