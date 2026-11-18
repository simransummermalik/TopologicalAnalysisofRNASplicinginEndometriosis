from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from strand_kmer_check import (  # noqa: E402
    call_library,
    classify_reads,
    reverse_complement,
    strand_kmers,
)

random.seed(7)
MRNA = "".join(random.choice("ACGT") for _ in range(600))


def fastq(reads: list[str]) -> list[str]:
    lines = []
    for index, read in enumerate(reads):
        lines += [f"@r{index}\n", read + "\n", "+\n", "I" * len(read) + "\n"]
    return lines


def fragments(count: int = 20) -> list[tuple[str, str]]:
    """Read 1 = antisense, read 2 = sense, as in an RF library."""
    pairs = []
    for start in range(0, count * 20, 20):
        fragment = MRNA[start : start + 200]
        pairs.append((reverse_complement(fragment)[:100], fragment[:100]))
    return pairs


class StrandCheckTests(unittest.TestCase):
    def setUp(self) -> None:
        self.sense, self.antisense = strand_kmers([MRNA])

    def counts(self, reads):
        return classify_reads(fastq(reads), self.sense, self.antisense)

    def test_rf_fr_and_unstranded_libraries(self) -> None:
        pairs = fragments()
        rf = [self.counts([a for a, _ in pairs]), self.counts([b for _, b in pairs])]
        self.assertEqual(rf[0]["antisense"], 20)
        self.assertEqual(rf[1]["sense"], 20)
        self.assertEqual(call_library(*rf), "RF")
        self.assertEqual(call_library(rf[1], rf[0]), "FR")
        mixed = self.counts([a for a, _ in pairs[:10]] + [b for _, b in pairs[10:]])
        self.assertEqual(call_library(mixed, mixed), "unresolved")

    def test_unrelated_reads_do_not_count(self) -> None:
        other = "".join(random.choice("ACGT") for _ in range(100))
        counts = self.counts([other])
        self.assertEqual((counts["reads"], counts["sense"], counts["antisense"]), (1, 0, 0))
        self.assertEqual(call_library(counts, counts), "unresolved")


if __name__ == "__main__":
    unittest.main()
