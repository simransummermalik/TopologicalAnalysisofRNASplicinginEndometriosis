import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hisat2_pipeline import checksum
from subsample_pairs import extract


class SubsampleTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.left, self.right = self.root / "left.fq", self.root / "right.fq"
        for path, mate in [(self.left, 1), (self.right, 2)]:
            path.write_text(f"@A/{mate}\nACGT\n+\nIIII\n@B/{mate}\nTTTT\n+\nJJJJ\n")

    def tearDown(self):
        self.temp.cleanup()

    def test_exact_prefix_preserves_mates_and_hashes_without_overwriting(self):
        output = self.root / "subset"
        record = extract(
            self.left,
            self.right,
            output,
            1,
            checksum(self.left)["sha256"],
            checksum(self.right)["sha256"],
        )
        self.assertEqual(record["status"], "complete")
        self.assertEqual(record["written_pairs"], 1)
        with gzip.open(output / "R1.fastq.gz", "rt") as handle:
            self.assertEqual(handle.read(), "@A/1\nACGT\n+\nIIII\n")
        self.assertEqual(
            record["outputs"][0]["sha256"], checksum(output / "R1.fastq.gz")["sha256"]
        )
        with self.assertRaises(FileExistsError):
            extract(self.left, self.right, output, 1)

    def test_bad_checksum_unsynchronized_mates_and_short_input_leave_failed_manifests(
        self,
    ):
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            extract(self.left, self.right, self.root / "bad-hash", 1, "0" * 64)
        with self.assertRaisesRegex(ValueError, "ended before"):
            extract(self.left, self.right, self.root / "short", 3)
        self.right.write_text(self.right.read_text().replace("@A/2", "@X/2"))
        with self.assertRaisesRegex(ValueError, "names disagree"):
            extract(self.left, self.right, self.root / "bad-mates", 1)
        for folder in ["bad-hash", "short", "bad-mates"]:
            self.assertEqual(
                json.loads((self.root / folder / "subset_manifest.json").read_text())[
                    "status"
                ],
                "failed",
            )


if __name__ == "__main__":
    unittest.main()
