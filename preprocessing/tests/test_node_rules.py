"""Junction-to-node rules decided 2026-11-18 (docs/data_methods.md)."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from build_junction_gene_map import build_gene_map  # noqa: E402
from hisat2_pipeline import splice_hints  # noqa: E402
from sam_junctions import count_junctions  # noqa: E402
from star_junctions import convert_star_junctions  # noqa: E402

# One plus-strand gene, five transcripts:
#   T1  A 100-199 | B 300-399 | C 500-599        full-length
#   T2  A'150-199 | B 300-399 | C'500-650        truncated start, longer end
#   T3  A 100-199 | C 500-599                    exon skipping
#   T4  A 100-199 | B'300-420 | C 500-599        alternative donor on B
#   T5  A 100-199 | R 300-599                    retained intron (R spans B..C)
# Junction reads see only splice sites, so A/A', B/B' and C/C' are one exon each.
EXONS = {
    "T1": [(100, 199), (300, 399), (500, 599)],
    "T2": [(150, 199), (300, 399), (500, 650)],
    "T3": [(100, 199), (500, 599)],
    "T4": [(100, 199), (300, 420), (500, 599)],
    "T5": [(100, 199), (300, 599)],
}
A, B, C = "chr1:100-199:+", "chr1:300-420:+", "chr1:500-650:+"


def write_gtf(path: Path) -> None:
    lines = []
    for transcript, exons in EXONS.items():
        for start, end in exons:
            lines.append(
                f'chr1\ttest\texon\t{start}\t{end}\t.\t+\t.\tgene_id "G"; transcript_id "{transcript}";\n'
            )
    # A minus-strand gene checks that from/to follow transcript order.
    for start, end in [(100, 149), (250, 299)]:
        lines.append(
            f'chr2\ttest\texon\t{start}\t{end}\t.\t-\t.\tgene_id "M"; transcript_id "TM";\n'
        )
    path.write_text("".join(lines), encoding="utf-8")


def load(path: Path) -> dict[tuple[str, str, str, str], set[tuple[str, str, str]]]:
    mapping = defaultdict(set)
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["chromosome"], row["intron_start"], row["intron_end"], row["strand"])
            mapping[key].add((row["gene_id"], row["from_node"], row["to_node"]))
    return mapping


def sam(start: int, gap: int, position: int) -> str:
    """One unique plus-strand read with 20-base anchors spanning one intron."""
    return f"r\t0\tchr1\t{position}\t60\t20M{gap}N20M\t*\t0\t0\t*\t*\tNH:i:1\tXS:A:+\n"


class NodeRuleTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.gtf = self.root / "annotation.gtf"
        write_gtf(self.gtf)
        self.shared = self.root / "shared.tsv"
        self.exact = self.root / "exact.tsv"
        self.shared_stats = build_gene_map(self.gtf, self.shared)
        self.exact_stats = build_gene_map(self.gtf, self.exact, "exon-interval")

    def tearDown(self) -> None:
        self.temp.cleanup()

    def test_shared_splice_site_nodes_make_every_junction_unique(self) -> None:
        mapping = load(self.shared)
        self.assertEqual(self.shared_stats["node_rule"], "shared-splice-site")
        self.assertEqual(
            {key: targets for key, targets in mapping.items() if key[0] == "chr1"},
            {
                ("chr1", "200", "299", "1"): {("G", A, B)},
                ("chr1", "400", "499", "1"): {("G", B, C)},
                ("chr1", "421", "499", "1"): {("G", B, C)},  # alternative donor
                ("chr1", "200", "499", "1"): {("G", A, C)},  # exon skipping
            },
        )
        self.assertEqual(
            mapping[("chr2", "150", "249", "2")], {("M", "chr2:250-299:-", "chr2:100-149:-")}
        )
        # R contains intron 400-499, so it may not merge B and C.
        self.assertEqual(self.shared_stats["retained_intron_exons_kept_apart"], 1)
        self.assertEqual(self.shared_stats["same_node_junction_rows"], 0)

    def test_exon_interval_rule_is_unchanged_and_ambiguous_here(self) -> None:
        mapping = load(self.exact)
        self.assertEqual(self.exact_stats["node_rule"], "exon-interval")
        self.assertEqual(len(mapping[("chr1", "200", "299", "1")]), 4)
        self.assertEqual(
            mapping[("chr1", "200", "499", "1")],
            {("G", "chr1:100-199:+", "chr1:500-599:+")},
        )

    def test_both_rules_give_identical_splice_site_hints(self) -> None:
        shared_hints, exact_hints = self.root / "h1.tsv", self.root / "h2.tsv"
        splice_hints(self.shared, shared_hints)
        splice_hints(self.exact, exact_hints)
        self.assertEqual(shared_hints.read_text(), exact_hints.read_text())

    def test_reads_count_and_parallel_junctions_sum_to_one_edge(self) -> None:
        lines = [
            sam(200, 100, 180),  # 200-299
            sam(200, 100, 180),
            sam(400, 100, 380),  # 400-499
            sam(421, 79, 401),  # 421-499, alternative donor, same node pair
            sam(200, 300, 180),  # 200-499, skipping
        ]
        counts, stats = count_junctions(lines, self.shared, "s", 8, 1)
        self.assertEqual(
            dict(counts), {("G", A, B): 2, ("G", B, C): 2, ("G", A, C): 1}
        )
        self.assertEqual(stats["assigned_junction_observations"], 5)
        self.assertNotIn("ambiguous_junction_observations", stats)
        _, exact_stats = count_junctions(lines, self.exact, "s", 8, 1)
        # Old rule: 200-299 (two reads) and 400-499 match several exon pairs.
        self.assertEqual(exact_stats["ambiguous_junction_observations"], 3)

    def test_same_node_junctions_are_excluded_with_their_own_label(self) -> None:
        gene_map = self.root / "loop.tsv"
        gene_map.write_text(
            "chromosome\tintron_start\tintron_end\tstrand\tgene_id\tfrom_node\tto_node\n"
            f"chr1\t200\t299\t1\tG\t{A}\t{A}\n"
            f"chr1\t400\t499\t1\tG\t{B}\t{C}\n",
            encoding="utf-8",
        )
        counts, stats = count_junctions(
            [sam(200, 100, 180), sam(400, 100, 380)], gene_map, "s", 8, 1
        )
        self.assertEqual(dict(counts), {("G", B, C): 1})
        self.assertEqual(stats["same_node_junction_observations"], 1)
        sj = self.root / "SJ.out.tab"
        sj.write_text("chr1\t200\t299\t1\t1\t1\t5\t0\t20\nchr1\t400\t499\t1\t1\t1\t4\t0\t20\n")
        star = convert_star_junctions(
            sj_path=sj, gene_map_path=gene_map, output_path=self.root / "out.tsv", sample_id="s"
        )
        self.assertEqual(star["same_node_junctions"], 1)
        self.assertEqual(star["written_rows"], 1)

    def test_unknown_rule_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            build_gene_map(self.gtf, self.root / "bad.tsv", "nearest")


if __name__ == "__main__":
    unittest.main()
