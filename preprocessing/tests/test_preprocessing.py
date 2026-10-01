from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

PREPROCESSING_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PREPROCESSING_DIR))

from build_junction_gene_map import build_gene_map  # noqa: E402
from run_preprocessing import (  # noqa: E402
    PipelineConfig,
    build_commands,
    output_paths,
    validate_star_outputs,
)
from star_junctions import convert_star_junctions  # noqa: E402


FIXTURES = Path(__file__).resolve().parent / "fixtures"


class PreprocessingTests(unittest.TestCase):
    def test_gtf_to_splice_girl_tsv(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            temp = Path(temp_directory)
            gene_map = temp / "junction_gene_map.tsv"
            junction_tsv = temp / "sample.junctions.tsv"
            qc_json = temp / "sample.junctions.json"

            map_stats = build_gene_map(FIXTURES / "annotation.gtf", gene_map)
            self.assertEqual(map_stats["junctions"], 3)

            stats = convert_star_junctions(
                sj_path=FIXTURES / "sample.SJ.out.tab",
                gene_map_path=gene_map,
                output_path=junction_tsv,
                sample_id="SMOKE",
                minimum_unique_reads=1,
                qc_path=qc_json,
            )
            self.assertEqual(stats["written_rows"], 3)
            self.assertEqual(stats["unannotated_junctions"], 1)
            self.assertEqual(stats["below_minimum_unique_reads"], 1)
            self.assertTrue(qc_json.is_file())

            with junction_tsv.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(
                rows,
                [
                    {
                        "sample_id": "SMOKE",
                        "gene_id": "GENE_A",
                        "from": "chr1:199:+",
                        "to": "chr1:300:+",
                        "count": "12",
                    },
                    {
                        "sample_id": "SMOKE",
                        "gene_id": "GENE_A",
                        "from": "chr1:399:+",
                        "to": "chr1:500:+",
                        "count": "8",
                    },
                    {
                        "sample_id": "SMOKE",
                        "gene_id": "GENE_B",
                        "from": "chr2:250:-",
                        "to": "chr2:149:-",
                        "count": "7",
                    },
                ],
            )

    def test_command_plan_supports_single_end_reads(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            temp = Path(temp_directory)
            gene_map = temp / "map.tsv"
            gene_map.write_text(
                "chromosome\tintron_start\tintron_end\tstrand\tgene_id\n",
                encoding="utf-8",
            )
            config = PipelineConfig(
                sample_id="SMOKE",
                read1=FIXTURES / "reads_R1.fastq",
                read2=None,
                star_index=FIXTURES / "star_index",
                gene_map=gene_map,
                output_dir=temp / "run",
                threads=2,
                phix_index=None,
                minimum_unique_reads=1,
                include_multimapping=False,
            )
            config.validate()
            commands = build_commands(config)
            self.assertEqual([command[0] for command in commands[:3]], ["cleaver", "STAR", "cleaver"])
            self.assertEqual(commands[0][1], "fastp")
            self.assertEqual(commands[2][1:3], ["samtools", "index"])
            self.assertIn("SortedByCoordinate", commands[1])
            self.assertEqual(commands[-1][1], str(PREPROCESSING_DIR / "star_junctions.py"))

    def test_command_plan_uses_cleaver_paired_end_options(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            temp = Path(temp_directory)
            config = PipelineConfig(
                sample_id="PAIRED",
                read1=FIXTURES / "reads_R1.fastq",
                read2=FIXTURES / "reads_R1.fastq",
                star_index=FIXTURES / "star_index",
                gene_map=FIXTURES / "annotation.gtf",
                output_dir=temp / "run",
                threads=2,
                phix_index=None,
                minimum_unique_reads=1,
                include_multimapping=False,
            )

            commands = build_commands(config)

            self.assertEqual(commands[0][:2], ["cleaver", "fastp"])
            self.assertIn("-I", commands[0])
            self.assertIn("-O", commands[0])
            self.assertIn("-2", commands[0])
            self.assertEqual(commands[1].count("--readFilesIn"), 1)
            self.assertTrue(
                any(value.endswith("PAIRED.R2.fastq.gz") for value in commands[1])
            )

    def test_star_zero_input_reads_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp_directory:
            temp = Path(temp_directory)
            config = PipelineConfig(
                sample_id="SMOKE",
                read1=FIXTURES / "reads_R1.fastq",
                read2=None,
                star_index=FIXTURES / "star_index",
                gene_map=temp / "map.tsv",
                output_dir=temp / "run",
                threads=2,
                phix_index=None,
                minimum_unique_reads=1,
                include_multimapping=False,
            )
            paths = output_paths(config)
            paths["star_log"].parent.mkdir(parents=True)
            paths["star_log"].write_text(
                "Number of input reads | 0\n", encoding="utf-8"
            )

            with self.assertRaisesRegex(ValueError, "zero input reads"):
                validate_star_outputs(config)


if __name__ == "__main__":
    unittest.main()
