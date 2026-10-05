#!/usr/bin/env python3
"""Validate the E-MTAB-15117 sample sheet and make a clean sample manifest."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

STUDY_ACCESSION = "E-MTAB-15117"
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = REPOSITORY_ROOT / "research" / "dataorg.csv"
DEFAULT_OUTPUT = REPOSITORY_ROOT / "data" / "metadata" / "E-MTAB-15117_samples.tsv"

REQUIRED_COLUMNS = {
    "Source Name",
    "Comment[ENA_SAMPLE]",
    "Comment[BioSD_SAMPLE]",
    "Characteristics[organism]",
    "Characteristics[age]",
    "Characteristics[sex]",
    "Characteristics[disease]",
    "Characteristics[individual]",
    "Characteristics[organism part]",
    "Comment[LIBRARY_LAYOUT]",
    "Comment[LIBRARY_SELECTION]",
    "Comment[LIBRARY_SOURCE]",
    "Comment[LIBRARY_STRATEGY]",
    "Comment[NOMINAL_LENGTH]",
    "Comment[ENA_EXPERIMENT]",
    "Comment[ENA_RUN]",
    "Comment[FASTQ_URI]",
}

OUTPUT_COLUMNS = [
    "study_accession",
    "sample_id",
    "condition",
    "original_disease",
    "tissue",
    "age_years",
    "sex",
    "individual_id",
    "ena_sample",
    "biosample",
    "ena_experiment",
    "ena_run",
    "library_layout",
    "library_selection",
    "library_source",
    "library_strategy",
    "nominal_length",
    "read1_url",
    "read2_url",
]

CONDITION_MAP = {
    "normal": "control",
    "Deep infiltration endometriosis (DIE)": "deep_infiltrating_endometriosis",
    # One source row contains this missing separator; preserve it in
    # original_disease while mapping it to the intended condition.
    "Deep infiltration endometriosis (DIE)endometriosis": "deep_infiltrating_endometriosis",
    "Ovarian endometriosis (OE)": "ovarian_endometriosis",
    "Peritoneal endometriosis (PE)": "peritoneal_endometriosis",
}


class ManifestError(ValueError):
    """Raised when the source sheet cannot be converted safely."""


def _required(row: dict[str, str], column: str, sample_id: str) -> str:
    value = (row.get(column) or "").strip()
    if not value:
        raise ManifestError(f"{sample_id}: missing {column}")
    return value


def _read_rows(input_path: Path) -> list[dict[str, str]]:
    try:
        with input_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = set(reader.fieldnames or [])
            missing = REQUIRED_COLUMNS - columns
            if missing:
                names = ", ".join(sorted(missing))
                raise ManifestError(f"source sheet is missing columns: {names}")
            rows = list(reader)
    except OSError as error:
        raise ManifestError(f"could not read {input_path}: {error}") from error

    if not rows:
        raise ManifestError(f"source sheet is empty: {input_path}")
    return rows


def _group_rows(rows: Iterable[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        sample_id = _required(row, "Source Name", "<unknown sample>")
        grouped[sample_id].append(row)
    return dict(grouped)


def _pair_urls(sample_id: str, rows: list[dict[str, str]]) -> tuple[str, str]:
    if len(rows) != 2:
        raise ManifestError(f"{sample_id}: expected two rows (R1 and R2), found {len(rows)}")

    urls: dict[str, str] = {}
    for row in rows:
        uri = _required(row, "Comment[FASTQ_URI]", sample_id)
        if uri.endswith("_1.fastq.gz"):
            read_name = "read1_url"
        elif uri.endswith("_2.fastq.gz"):
            read_name = "read2_url"
        else:
            raise ManifestError(f"{sample_id}: FASTQ URI is not an R1/R2 file: {uri}")
        if read_name in urls:
            raise ManifestError(f"{sample_id}: duplicate {read_name}")
        urls[read_name] = uri

    if set(urls) != {"read1_url", "read2_url"}:
        raise ManifestError(f"{sample_id}: both R1 and R2 FASTQ files are required")
    return urls["read1_url"], urls["read2_url"]


def _check_consistent(sample_id: str, rows: list[dict[str, str]]) -> None:
    fields = [
        "Comment[ENA_SAMPLE]",
        "Comment[BioSD_SAMPLE]",
        "Characteristics[organism]",
        "Characteristics[age]",
        "Characteristics[sex]",
        "Characteristics[disease]",
        "Characteristics[individual]",
        "Characteristics[organism part]",
        "Comment[LIBRARY_LAYOUT]",
        "Comment[LIBRARY_SELECTION]",
        "Comment[LIBRARY_SOURCE]",
        "Comment[LIBRARY_STRATEGY]",
        "Comment[NOMINAL_LENGTH]",
        "Comment[ENA_EXPERIMENT]",
        "Comment[ENA_RUN]",
    ]
    for field in fields:
        values = {_required(row, field, sample_id) for row in rows}
        if len(values) != 1:
            raise ManifestError(f"{sample_id}: R1/R2 disagree for {field}")


def _make_record(sample_id: str, rows: list[dict[str, str]]) -> dict[str, str]:
    _check_consistent(sample_id, rows)
    row = rows[0]
    original_disease = _required(row, "Characteristics[disease]", sample_id)
    try:
        condition = CONDITION_MAP[original_disease]
    except KeyError as error:
        raise ManifestError(
            f"{sample_id}: unrecognized disease label: {original_disease}"
        ) from error

    read1_url, read2_url = _pair_urls(sample_id, rows)
    age = _required(row, "Characteristics[age]", sample_id)
    try:
        age_years = str(int(age))
    except ValueError as error:
        raise ManifestError(f"{sample_id}: age is not an integer: {age}") from error

    layout = _required(row, "Comment[LIBRARY_LAYOUT]", sample_id)
    if layout != "PAIRED":
        raise ManifestError(f"{sample_id}: expected PAIRED layout, found {layout}")

    return {
        "study_accession": STUDY_ACCESSION,
        "sample_id": sample_id,
        "condition": condition,
        "original_disease": original_disease,
        "tissue": _required(row, "Characteristics[organism part]", sample_id),
        "age_years": age_years,
        "sex": _required(row, "Characteristics[sex]", sample_id),
        "individual_id": _required(row, "Characteristics[individual]", sample_id),
        "ena_sample": _required(row, "Comment[ENA_SAMPLE]", sample_id),
        "biosample": _required(row, "Comment[BioSD_SAMPLE]", sample_id),
        "ena_experiment": _required(row, "Comment[ENA_EXPERIMENT]", sample_id),
        "ena_run": _required(row, "Comment[ENA_RUN]", sample_id),
        "library_layout": layout,
        "library_selection": _required(row, "Comment[LIBRARY_SELECTION]", sample_id),
        "library_source": _required(row, "Comment[LIBRARY_SOURCE]", sample_id),
        "library_strategy": _required(row, "Comment[LIBRARY_STRATEGY]", sample_id),
        "nominal_length": _required(row, "Comment[NOMINAL_LENGTH]", sample_id),
        "read1_url": read1_url,
        "read2_url": read2_url,
    }


def build_manifest(input_path: Path) -> list[dict[str, str]]:
    rows = _read_rows(input_path)
    grouped = _group_rows(rows)
    records = [_make_record(sample_id, grouped[sample_id]) for sample_id in sorted(grouped)]

    runs = [record["ena_run"] for record in records]
    if len(runs) != len(set(runs)):
        raise ManifestError("a run accession is assigned to more than one sample")
    return records


def write_manifest(records: list[dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(records)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate dataorg.csv and write a one-row-per-sample manifest."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="validate the source sheet without writing a manifest",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        records = build_manifest(args.input)
        if not args.check_only:
            write_manifest(records, args.output)
        counts: dict[str, int] = defaultdict(int)
        for record in records:
            counts[record["condition"]] += 1
        print(f"validated {len(records)} samples from {args.input}")
        print("conditions: " + ", ".join(f"{key}={counts[key]}" for key in sorted(counts)))
        if not args.check_only:
            print(f"wrote {args.output}")
    except ManifestError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
