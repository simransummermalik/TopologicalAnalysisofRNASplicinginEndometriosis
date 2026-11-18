from __future__ import annotations

import argparse
import contextlib
import io
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_junction_gene_map import build_gene_map
from hisat2_pipeline import configuration, index_files, reference, splice_hints
from sam_junctions import convert_sam, count_junctions, skipped_intervals

FIXTURES = Path(__file__).parent / "fixtures"


def sam(cigar="50M100N50M", position=150, flag=0, tags="NH:i:1\tXS:A:+", chrom="chr1"):
    return f"read\t{flag}\t{chrom}\t{position}\t60\t{cigar}\t*\t0\t0\t*\t*\t{tags}\n"


class HisatTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.map = self.root / "map.tsv"
        build_gene_map(FIXTURES / "annotation.gtf", self.map)

    def tearDown(self):
        self.temporary.cleanup()

    def test_cigar_reference_coordinates_include_deletions_not_insertions_or_clipping(
        self,
    ):
        self.assertEqual(skipped_intervals(150, "5S10M2I40M100N50M3S", 8), [(200, 299)])
        self.assertEqual(skipped_intervals(145, "10M5D40M100N50M", 8), [(200, 299)])
        self.assertEqual(
            skipped_intervals(150, "50M100N100M100N50M", 8), [(200, 299), (400, 499)]
        )
        self.assertEqual(skipped_intervals(195, "5M100N50M", 8), [])
        with self.assertRaisesRegex(ValueError, "invalid"):
            skipped_intervals(150, "50M100N50Mbad", 8)

    def test_plus_and_minus_direction_comes_from_transcript_not_read_strand(self):
        lines = [
            sam(flag=16),
            sam("30M100N30M", 120, tags="NH:i:1\tXS:A:-", chrom="chr2"),
        ]
        counts, stats = count_junctions(lines, self.map, "sample", 8, 1)
        self.assertEqual(counts[("GENE_A", "chr1:100-199:+", "chr1:300-399:+")], 1)
        self.assertEqual(counts[("GENE_B", "chr2:250-299:-", "chr2:100-149:-")], 1)
        self.assertEqual(stats["assigned_junction_observations"], 2)

    def test_nonprimary_multimapped_unknown_strand_and_low_anchor_are_explicit(self):
        lines = [sam(flag=flag) for flag in [4, 256, 512, 2048]]
        lines += [
            sam(tags="NH:i:2\tXS:A:+"),
            sam(tags="XS:A:+"),
            sam(tags="NH:i:1"),
            sam("5M100N50M", 195),
            sam(),
        ]
        counts, stats = count_junctions(lines, self.map, "sample", 8, 1)
        self.assertEqual(sum(counts.values()), 1)
        self.assertEqual(stats["unmapped_secondary_qcfail_or_supplementary"], 4)
        self.assertEqual(stats["nonunique_or_missing_NH_records"], 2)
        self.assertEqual(stats["unknown_transcript_strand_records"], 1)
        self.assertEqual(stats["insufficient_anchor_records"], 1)

    def test_support_threshold_and_duplicate_policy_are_applied_after_counting(self):
        counts, _ = count_junctions([sam(), sam(flag=1024)], self.map, "sample", 8, 2)
        self.assertEqual(sum(counts.values()), 2)
        counts, stats = count_junctions(
            [sam(), sam(flag=1024)], self.map, "sample", 8, 2, True
        )
        self.assertEqual(counts, {})
        self.assertEqual(stats["below_minimum_unique_edges"], 1)
        self.assertEqual(stats["duplicate_records_excluded"], 1)

    def test_converter_protects_outputs_and_leaves_qc_on_no_eligible_junctions(self):
        source, output, qc = [
            self.root / name for name in ["reads.sam", "junctions.tsv", "qc.json"]
        ]
        source.write_text(sam())
        convert_sam(source, self.map, "sample", output, qc, 8, 1)
        with self.assertRaisesRegex(ValueError, "overwrite"):
            convert_sam(source, self.map, "sample", output, qc, 8, 1)
        missing, failed_qc = self.root / "none.tsv", self.root / "none.json"
        with self.assertRaisesRegex(ValueError, "no uniquely"):
            convert_sam(source, self.map, "sample", missing, failed_qc, 8, 2)
        self.assertFalse(missing.exists())
        self.assertEqual(json.loads(failed_qc.read_text())["written_rows"], 0)

    def test_hints_convert_to_zero_based_flanking_bases(self):
        hints = self.root / "hints.tsv"
        splice_hints(self.map, hints)
        self.assertIn("chr1\t198\t299\t+", hints.read_text())
        self.assertIn("chr2\t148\t249\t-", hints.read_text())

    def test_complete_index_and_explicit_paired_plan(self):
        prefix = self.root / "genome"
        with self.assertRaisesRegex(ValueError, "incomplete"):
            index_files(prefix)
        for number in range(1, 9):
            Path(f"{prefix}.{number}.ht2").write_text("fixture")
        read2 = self.root / "R2.fastq"
        read2.write_text("fixture")
        args = argparse.Namespace(
            sample_id="sample",
            read1=FIXTURES / "reads_R1.fastq",
            read2=read2,
            index=prefix,
            gene_map=self.map,
            splice_sites=self.map,
            output_dir=self.root / "run",
            threads=2,
            min_anchor=8,
            min_unique=1,
            strandedness="RF",
            skip_cleaning=False,
            exclude_duplicates=False,
        )
        _, commands = configuration(args)
        self.assertEqual(commands[0][:2], ["cleaver", "fastp"])
        self.assertIn("-2", commands[0])
        self.assertIn("--known-splicesite-infile", commands[1])
        self.assertIn("--dta-cufflinks", commands[1])
        self.assertIn("RF", commands[1])
        self.assertEqual(commands[-2][:3], ["cleaver", "samtools", "sort"])
        args.strandedness = "F"
        with self.assertRaisesRegex(ValueError, "paired reads"):
            configuration(args)

    def test_reference_cleanup_removes_only_its_generated_copy_after_success(self):
        output = self.root / "reference"
        args = argparse.Namespace(
            command="prepare-reference",
            fasta=FIXTURES / "tiny_genome.fa",
            gtf=FIXTURES / "tiny_annotation.gtf",
            output_dir=output,
            threads=2,
            build_bmaxdivn=16,
            discard_expanded_reference=True,
        )
        source = args.fasta.read_bytes()

        def build(command, **kwargs):
            self.assertIn("--bmaxdivn", command)
            for number in range(1, 9):
                (output / f"genome.{number}.ht2").write_text("fixture")

        with (
            patch("hisat2_pipeline.require_tools"),
            patch("hisat2_pipeline.versions", return_value={}),
            patch("hisat2_pipeline.subprocess.run", side_effect=build),
            contextlib.redirect_stdout(io.StringIO()),
        ):
            reference(args)
        self.assertFalse((output / "reference.fa").exists())
        self.assertEqual(args.fasta.read_bytes(), source)
        self.assertFalse(
            json.loads((output / "reference_manifest.json").read_text())[
                "expanded_reference_retained"
            ]
        )

    def test_partial_or_still_running_reference_is_not_accepted_for_alignment(self):
        prefix = self.root / "building"
        for number in range(1, 9):
            Path(f"{prefix}.{number}.ht2").write_text("fixture")
        Path(f"{prefix}.8.ht2").write_text("")
        with self.assertRaisesRegex(ValueError, "incomplete"):
            index_files(prefix)
        Path(f"{prefix}.8.ht2").write_text("fixture")
        manifest = self.root / "reference_manifest.json"
        manifest.write_text('{"status":"running"}')
        with self.assertRaisesRegex(ValueError, "not complete"):
            index_files(prefix)
        manifest.write_text('{"status":"complete"}')
        self.assertEqual(len(index_files(prefix)), 8)

    def test_reference_build_failure_preserves_diagnostics_and_expanded_copy(self):
        output = self.root / "failed-reference"
        args = argparse.Namespace(
            command="prepare-reference",
            fasta=FIXTURES / "tiny_genome.fa",
            gtf=FIXTURES / "tiny_annotation.gtf",
            output_dir=output,
            threads=2,
            build_bmaxdivn=16,
            discard_expanded_reference=True,
        )
        with (
            patch("hisat2_pipeline.require_tools"),
            patch("hisat2_pipeline.versions", return_value={}),
            patch(
                "hisat2_pipeline.subprocess.run",
                side_effect=OSError("index build failed"),
            ),
            self.assertRaisesRegex(OSError, "build failed"),
        ):
            reference(args)
        self.assertTrue((output / "reference.fa").exists())
        self.assertEqual(
            json.loads((output / "reference_manifest.json").read_text())["status"],
            "failed",
        )


if __name__ == "__main__":
    unittest.main()
