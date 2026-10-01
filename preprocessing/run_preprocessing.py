#!/usr/bin/env python3
"""Run one raw-read preprocessing workflow for Splice Girl."""

from __future__ import annotations

import argparse
import json
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
CONVERTER = SCRIPT_DIR / "star_junctions.py"
REQUIRED_TOOLS = ("cleaver", "STAR")


@dataclass(frozen=True)
class PipelineConfig:
    sample_id: str
    read1: Path
    read2: Path | None
    star_index: Path
    gene_map: Path
    output_dir: Path
    threads: int
    phix_index: str | None
    minimum_unique_reads: int
    include_multimapping: bool

    def validate(self) -> None:
        if not re.fullmatch(r"[A-Za-z0-9._-]+", self.sample_id):
            raise ValueError(
                "sample ID may contain only letters, numbers, dot, underscore, "
                "and hyphen"
            )
        if not self.read1.is_file():
            raise ValueError(f"read 1 file does not exist: {self.read1}")
        if self.read2 is not None and not self.read2.is_file():
            raise ValueError(f"read 2 file does not exist: {self.read2}")
        if not self.star_index.is_dir():
            raise ValueError(f"STAR index directory does not exist: {self.star_index}")
        if not self.gene_map.is_file():
            raise ValueError(f"junction gene map does not exist: {self.gene_map}")
        if self.threads < 1:
            raise ValueError("threads must be at least 1")
        if self.minimum_unique_reads < 0:
            raise ValueError("minimum unique reads must be zero or greater")


def output_paths(config: PipelineConfig) -> dict[str, Path]:
    sample = config.sample_id
    output = config.output_dir
    return {
        "trimmed_r1": output / "trimmed" / f"{sample}.R1.fastq.gz",
        "trimmed_r2": output / "trimmed" / f"{sample}.R2.fastq.gz",
        "clean_r1": output / "cleaned" / f"{sample}.R1.fastq.gz",
        "clean_r2": output / "cleaned" / f"{sample}.R2.fastq.gz",
        "fastp_json": output / "qc" / f"{sample}.fastp.json",
        "star_prefix": output / "aligned" / f"{sample}.",
        "star_log": output / "aligned" / f"{sample}.Log.final.out",
        "bam": output / "aligned" / f"{sample}.Aligned.sortedByCoord.out.bam",
        "sj": output / "aligned" / f"{sample}.SJ.out.tab",
        "junction_tsv": output / "splice_girl" / f"{sample}.junctions.tsv",
        "junction_qc": output / "qc" / f"{sample}.junctions.json",
        "manifest": output / "run_manifest.json",
    }


def build_commands(config: PipelineConfig) -> list[list[str]]:
    paths = output_paths(config)
    cleaver_fastp = [
        "cleaver",
        "fastp",
        "-i",
        str(config.read1),
        "-o",
        str(paths["trimmed_r1"]),
        "-j",
        str(paths["fastp_json"]),
    ]
    if config.read2 is not None:
        cleaver_fastp.extend(
            [
                "-I",
                str(config.read2),
                "-O",
                str(paths["trimmed_r2"]),
                "-2",
            ]
        )

    commands = [cleaver_fastp]
    star_r1 = paths["trimmed_r1"]
    star_r2 = paths["trimmed_r2"] if config.read2 is not None else None

    if config.phix_index:
        if config.read2 is not None:
            phix = [
                "bowtie2",
                "-x",
                config.phix_index,
                "-1",
                str(paths["trimmed_r1"]),
                "-2",
                str(paths["trimmed_r2"]),
                "--very-sensitive-local",
                "--un-conc-gz",
                str(config.output_dir / "cleaned" / f"{config.sample_id}.R%.fastq.gz"),
                "-S",
                "/dev/null",
                "-p",
                str(config.threads),
            ]
            star_r2 = paths["clean_r2"]
        else:
            phix = [
                "bowtie2",
                "-x",
                config.phix_index,
                "-U",
                str(paths["trimmed_r1"]),
                "--very-sensitive-local",
                "--un-gz",
                str(paths["clean_r1"]),
                "-S",
                "/dev/null",
                "-p",
                str(config.threads),
            ]
        commands.append(phix)
        star_r1 = paths["clean_r1"]

    star = [
        "STAR",
        "--runThreadN",
        str(config.threads),
        "--genomeDir",
        str(config.star_index),
        "--readFilesIn",
        str(star_r1),
    ]
    if star_r2 is not None:
        star.append(str(star_r2))
    star.extend(
        [
            "--readFilesCommand",
            "gzip -cd",
            "--outFileNamePrefix",
            str(paths["star_prefix"]),
            "--outSAMtype",
            "BAM",
            "SortedByCoordinate",
        ]
    )
    commands.append(star)
    commands.append(["cleaver", "samtools", "index", str(paths["bam"])])

    converter = [
        sys.executable,
        str(CONVERTER),
        "--sj",
        str(paths["sj"]),
        "--gene-map",
        str(config.gene_map),
        "--sample-id",
        config.sample_id,
        "--output",
        str(paths["junction_tsv"]),
        "--qc-json",
        str(paths["junction_qc"]),
        "--min-unique",
        str(config.minimum_unique_reads),
    ]
    if config.include_multimapping:
        converter.append("--include-multimapping")
    commands.append(converter)
    return commands


def tool_status(include_phix: bool) -> dict[str, str | None]:
    names = [*REQUIRED_TOOLS]
    if include_phix:
        names.append("bowtie2")
    return {name: shutil.which(name) for name in names}


def tool_details(include_phix: bool) -> dict[str, dict[str, str | None]]:
    version_arguments = {
        "cleaver": ["--version"],
        "STAR": ["--version"],
        "bowtie2": ["--version"],
    }
    details: dict[str, dict[str, str | None]] = {}
    for name, path in tool_status(include_phix).items():
        version = None
        if path is not None:
            try:
                result = subprocess.run(
                    [path, *version_arguments[name]],
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=10,
                )
                combined = "\n".join(
                    part.strip() for part in (result.stdout, result.stderr) if part.strip()
                )
                version = combined.splitlines()[0] if combined else "unknown"
            except (OSError, subprocess.TimeoutExpired):
                version = "unknown"
        details[name] = {"path": path, "version": version}
    return details


def command_plan(config: PipelineConfig) -> str:
    return "\n".join(f"{index}. {shlex.join(command)}" for index, command in enumerate(build_commands(config), 1))


def validate_star_outputs(config: PipelineConfig) -> None:
    """Reject a nominally successful STAR run that did not process reads."""
    paths = output_paths(config)
    log_path = paths["star_log"]
    if not log_path.is_file():
        raise ValueError(f"STAR did not create its final log: {log_path}")

    input_reads: int | None = None
    for line in log_path.read_text(encoding="utf-8").splitlines():
        if "Number of input reads" not in line:
            continue
        _, separator, value = line.partition("|")
        if not separator:
            continue
        try:
            input_reads = int(value.strip().replace(",", ""))
        except ValueError as error:
            raise ValueError(
                f"could not read STAR's input-read count from {log_path}"
            ) from error
        break

    if input_reads is None:
        raise ValueError(f"STAR's final log has no input-read count: {log_path}")
    if input_reads == 0:
        raise ValueError(
            "STAR reported zero input reads; stop before BAM indexing and inspect "
            f"{log_path}"
        )
    for required_path in (paths["bam"], paths["sj"]):
        if not required_path.is_file():
            raise ValueError(f"STAR did not create expected output: {required_path}")


def write_manifest(
    config: PipelineConfig,
    commands: list[list[str]],
    status: str,
    error: str | None = None,
) -> None:
    paths = output_paths(config)
    serializable_config = {
        key: str(value) if isinstance(value, Path) else value
        for key, value in asdict(config).items()
    }
    manifest = {
        "pipeline": "Splice Girl RNA-seq preprocessing",
        "status": status,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "configuration": serializable_config,
        "tools": tool_details(config.phix_index is not None),
        "commands": commands,
        "outputs": {key: str(value) for key, value in paths.items()},
    }
    if error:
        manifest["error"] = error
    paths["manifest"].parent.mkdir(parents=True, exist_ok=True)
    with paths["manifest"].open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")


def execute(config: PipelineConfig) -> None:
    statuses = tool_status(config.phix_index is not None)
    missing = [name for name, path in statuses.items() if path is None]
    if missing:
        raise ValueError(
            "missing required command(s): " + ", ".join(missing) + ". See preprocessing/README.md."
        )

    paths = output_paths(config)
    protected_outputs = [
        paths["trimmed_r1"],
        paths["fastp_json"],
        paths["star_log"],
        paths["bam"],
        paths["sj"],
        paths["junction_tsv"],
    ]
    if config.read2 is not None:
        protected_outputs.append(paths["trimmed_r2"])
    existing = [path for path in protected_outputs if path.exists()]
    if existing:
        raise ValueError(
            "refusing to overwrite existing output(s): "
            + ", ".join(str(path) for path in existing)
        )

    for directory in ("trimmed", "cleaned", "aligned", "qc", "splice_girl"):
        (config.output_dir / directory).mkdir(parents=True, exist_ok=True)

    commands = build_commands(config)
    write_manifest(config, commands, "running")
    try:
        for number, command in enumerate(commands, start=1):
            print(f"[{number}/{len(commands)}] {shlex.join(command)}", flush=True)
            subprocess.run(command, check=True)
            if command[0] == "STAR":
                validate_star_outputs(config)
    except (OSError, subprocess.CalledProcessError, ValueError) as error:
        write_manifest(config, commands, "failed", str(error))
        raise ValueError(f"preprocessing stopped: {error}") from error
    write_manifest(config, commands, "complete")
    print(f"Splice Girl input: {paths['junction_tsv']}")


def add_pipeline_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--read1", required=True, type=Path)
    parser.add_argument("--read2", type=Path)
    parser.add_argument("--star-index", required=True, type=Path)
    parser.add_argument("--gene-map", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--phix-index",
        help="Optional Bowtie2 PhiX index prefix; omit when PhiX filtering is not approved.",
    )
    parser.add_argument("--min-unique", type=int, default=1)
    parser.add_argument("--include-multimapping", action="store_true")


def config_from_args(args: argparse.Namespace) -> PipelineConfig:
    return PipelineConfig(
        sample_id=args.sample_id,
        read1=args.read1.resolve(),
        read2=args.read2.resolve() if args.read2 else None,
        star_index=args.star_index.resolve(),
        gene_map=args.gene_map.resolve(),
        output_dir=args.output_dir.resolve(),
        threads=args.threads,
        phix_index=args.phix_index,
        minimum_unique_reads=args.min_unique,
        include_multimapping=args.include_multimapping,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="FASTQ-to-junction preprocessing for Splice Girl."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    check = subparsers.add_parser("check", help="Report required installed tools.")
    check.add_argument("--with-phix", action="store_true")
    add_pipeline_arguments(
        subparsers.add_parser("plan", help="Validate inputs and print commands only.")
    )
    add_pipeline_arguments(
        subparsers.add_parser("run", help="Execute the preprocessing pipeline.")
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.command == "check":
        statuses = tool_status(args.with_phix)
        for name, path in statuses.items():
            print(f"{name}: {path or 'missing'}")
        return 0 if all(statuses.values()) else 1

    config = config_from_args(args)
    try:
        config.validate()
        if args.command == "plan":
            print(command_plan(config))
        else:
            execute(config)
    except ValueError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
