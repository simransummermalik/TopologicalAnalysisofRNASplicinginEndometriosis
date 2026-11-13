#!/usr/bin/env python3
"""Map, deduplicate, check, and combine Splice Girl candidate-gene tables.

These helpers support Priority 2 (tasks G3-G10 and E2-E10). They never decide
which evidence counts: people fill the evidence tables by following
``docs/gene_set_methods.md``. The code only maps identifiers, flags anything
uncertain for human review, builds membership files, checks schemas, and
computes the risk/environment overlap.

Subcommands:

    map      Map original gene labels to versioned GENCODE gene IDs.
    dedupe   Build one membership row per gene from mapped evidence rows.
    check    Schema-check any gene-set table.
    overlap  Build overlap.tsv and candidate_genes.tsv from both memberships.

Run ``python3 gene_sets/gene_sets.py <subcommand> --help`` for options.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, TextIO

# ---------------------------------------------------------------------------
# Schemas and vocabularies
# ---------------------------------------------------------------------------

MISSING = "NA"

RISK_EVIDENCE_COLUMNS = [
    "evidence_id",
    "source",
    "source_version",
    "source_record",
    "original_id",
    "variant_id",
    "p_value",
    "evidence_type",
    "evidence_level",
    "population",
    "gene_id",
    "gene_symbol",
    "mapping_status",
    "mapping_method",
    "notes",
]

ENVIRONMENT_EVIDENCE_COLUMNS = [
    "evidence_id",
    "source",
    "source_version",
    "source_record",
    "publications",
    "original_id",
    "exposure",
    "exposure_id",
    "dose_duration",
    "species",
    "tissue",
    "assay",
    "response_type",
    "evidence_type",
    "context",
    "gene_id",
    "gene_symbol",
    "mapping_status",
    "mapping_method",
    "notes",
]

RISK_MEMBERSHIP_COLUMNS = [
    "gene_id",
    "gene_symbol",
    "evidence_type",
    "source",
    "n_evidence",
    "evidence_ids",
    "notes",
]

ENVIRONMENT_MEMBERSHIP_COLUMNS = [
    "gene_id",
    "gene_symbol",
    "exposure",
    "exposure_id",
    "tissue",
    "context",
    "evidence_type",
    "source",
    "n_publications",
    "n_evidence",
    "evidence_ids",
    "notes",
]

OVERLAP_COLUMNS = ["gene_id", "gene_symbol", "risk_source", "environment_source"]

CANDIDATE_COLUMNS = [
    "gene_id",
    "gene_symbol",
    "is_risk",
    "is_environment",
    "is_overlap",
    "category",
]

SCHEMAS: dict[str, list[str]] = {
    "risk_evidence": RISK_EVIDENCE_COLUMNS,
    "environment_evidence": ENVIRONMENT_EVIDENCE_COLUMNS,
    "risk_membership": RISK_MEMBERSHIP_COLUMNS,
    "environment_membership": ENVIRONMENT_MEMBERSHIP_COLUMNS,
    "overlap": OVERLAP_COLUMNS,
    "candidates": CANDIDATE_COLUMNS,
}

MAPPING_STATUSES = {"pending", "mapped", "ambiguous", "unresolved", "excluded"}
RISK_EVIDENCE_TYPES = {"gwas_association"}
RISK_EVIDENCE_LEVELS = {"variant", "locus", "gene"}
ENVIRONMENT_EVIDENCE_TYPES = {"measured_response"}
RESPONSE_TYPES = {"expression", "protein", "activity"}
CATEGORIES = {"risk_only", "environment_only", "overlap"}

# Columns that may hold the explicit missing token "NA". Every other column
# must contain a real value. gene_id and gene_symbol are additionally allowed
# to be NA unless mapping_status is "mapped".
OPTIONAL_COLUMNS: dict[str, set[str]] = {
    # original_id may be NA when the source reports a locus with no gene; map
    # then marks the row unresolved so it is reviewed, never silently assigned.
    "risk_evidence": {"original_id", "variant_id", "p_value", "population", "notes"},
    "environment_evidence": {
        "dose_duration",
        "tissue",
        "assay",
        "context",
        "notes",
    },
    "risk_membership": {"notes"},
    "environment_membership": {"tissue", "context", "notes"},
    "overlap": set(),
    "candidates": set(),
}

GENE_ID_PATTERN = re.compile(r"^ENSG\d{11}\.\d+(_PAR_Y)?$")
UNVERSIONED_ENSG = re.compile(r"^ENSG\d{11}$")
LIST_SEPARATOR = "|"
# A source label naming more than one gene, e.g. "GREB1, KDR" or "WNT4/CDC42".
MULTI_GENE_SPLIT = re.compile(r"\s*[,;/]\s*")


# ---------------------------------------------------------------------------
# TSV helpers
# ---------------------------------------------------------------------------


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV file with a header row. Lines starting with '#' are skipped."""
    with open(path, newline="", encoding="utf-8") as handle:
        lines = (line for line in handle if not line.startswith("#"))
        reader = csv.DictReader(lines, delimiter="\t")
        rows = [dict(row) for row in reader]
        header = list(reader.fieldnames or [])
    return header, rows


def write_tsv(path: Path, columns: list[str], rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def join_unique(values: Iterable[str]) -> str:
    seen: list[str] = []
    for value in values:
        if value and value != MISSING and value not in seen:
            seen.append(value)
    return LIST_SEPARATOR.join(seen) if seen else MISSING


# ---------------------------------------------------------------------------
# Reference loading
# ---------------------------------------------------------------------------


@dataclass
class GeneReference:
    """Gene IDs and names from the GENCODE GTF plus HGNC symbol history."""

    # versioned gene_id -> gene_name
    gene_names: dict[str, str] = field(default_factory=dict)
    # unversioned ENSG -> versioned gene_ids (PAR_Y copies give more than one)
    by_unversioned: dict[str, list[str]] = field(default_factory=dict)
    # gene_name -> versioned gene_ids
    by_name: dict[str, list[str]] = field(default_factory=dict)
    # HGNC approved / previous / alias symbol (upper case) -> unversioned ENSGs
    hgnc_approved: dict[str, set[str]] = field(default_factory=dict)
    hgnc_previous: dict[str, set[str]] = field(default_factory=dict)
    hgnc_alias: dict[str, set[str]] = field(default_factory=dict)


def parse_gtf_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        key, _, value = item.strip().partition(" ")
        if key and key not in attributes:
            attributes[key] = value.strip().strip('"')
    return attributes


def load_gtf_genes(path: Path, reference: GeneReference) -> None:
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attributes = parse_gtf_attributes(fields[8])
            gene_id = attributes.get("gene_id")
            gene_name = attributes.get("gene_name", "")
            if not gene_id:
                continue
            reference.gene_names[gene_id] = gene_name
            unversioned = gene_id.split(".")[0]
            reference.by_unversioned.setdefault(unversioned, []).append(gene_id)
            if gene_name:
                reference.by_name.setdefault(gene_name, []).append(gene_id)


def _split_hgnc_list(value: str) -> list[str]:
    value = value.strip().strip('"')
    return [item.strip() for item in value.split("|") if item.strip()]


def load_hgnc(path: Path, reference: GeneReference) -> None:
    """Load an HGNC complete-set style TSV.

    Required columns: symbol, alias_symbol, prev_symbol, ensembl_gene_id.
    Multi-valued cells use '|' and may be wrapped in double quotes.
    """
    header, rows = read_tsv(path)
    required = {"symbol", "alias_symbol", "prev_symbol", "ensembl_gene_id"}
    missing = required - set(header)
    if missing:
        raise ValueError(f"{path}: HGNC file is missing columns {sorted(missing)}")
    for row in rows:
        ensembl = (row.get("ensembl_gene_id") or "").strip()
        if not UNVERSIONED_ENSG.match(ensembl):
            continue
        for symbol in _split_hgnc_list(row.get("symbol", "")):
            reference.hgnc_approved.setdefault(symbol.upper(), set()).add(ensembl)
        for symbol in _split_hgnc_list(row.get("prev_symbol", "")):
            reference.hgnc_previous.setdefault(symbol.upper(), set()).add(ensembl)
        for symbol in _split_hgnc_list(row.get("alias_symbol", "")):
            reference.hgnc_alias.setdefault(symbol.upper(), set()).add(ensembl)


def load_reference(gtf: Path, hgnc: Path | None) -> GeneReference:
    reference = GeneReference()
    load_gtf_genes(gtf, reference)
    if hgnc is not None:
        load_hgnc(hgnc, reference)
    return reference


# ---------------------------------------------------------------------------
# map
# ---------------------------------------------------------------------------


@dataclass
class MappingResult:
    status: str
    method: str
    candidates: list[str]

    def gene_id(self) -> str:
        return self.candidates[0] if self.status == "mapped" else MISSING


def _versioned(reference: GeneReference, unversioned: Iterable[str]) -> list[str]:
    found: list[str] = []
    for ensg in sorted(unversioned):
        found.extend(reference.by_unversioned.get(ensg, []))
    return found


def _primary(gene_ids: list[str]) -> list[str]:
    """Drop PAR_Y copies when the primary copy of the same gene is present."""
    primary = [gid for gid in gene_ids if not gid.endswith("_PAR_Y")]
    return primary or gene_ids


def map_label(label: str, reference: GeneReference) -> MappingResult:
    """Map one original label. Never guesses: anything not unique is flagged."""
    label = label.strip()
    if not label or label == MISSING:
        return MappingResult("unresolved", "empty_label", [])

    parts = [part for part in MULTI_GENE_SPLIT.split(label) if part]
    if len(parts) > 1:
        candidates: list[str] = []
        for part in parts:
            result = map_label(part, reference)
            candidates.extend(result.candidates)
        return MappingResult("ambiguous", "multi_gene_label", candidates)

    unversioned = label.split(".")[0]
    if UNVERSIONED_ENSG.match(unversioned):
        ids = _primary(reference.by_unversioned.get(unversioned, []))
        if len(ids) == 1:
            return MappingResult("mapped", "gencode_id", ids)
        if ids:
            return MappingResult("ambiguous", "gencode_id", ids)
        return MappingResult("unresolved", "gencode_id_not_in_release", [])

    ids = _primary(reference.by_name.get(label, []))
    if len(ids) == 1:
        return MappingResult("mapped", "gencode_symbol", ids)
    if len(ids) > 1:
        # GENCODE sometimes gives one name to several gene IDs. HGNC, the
        # naming authority, links an approved symbol to one Ensembl gene; use
        # that link only when it picks out exactly one of the candidates.
        approved = reference.hgnc_approved.get(label.upper(), set())
        linked = [gid for gid in ids if gid.split(".")[0] in approved]
        if len(approved) == 1 and len(linked) == 1:
            return MappingResult("mapped", "gencode_symbol_hgnc_link", linked)
        return MappingResult("ambiguous", "gencode_symbol_duplicate", ids)

    key = label.upper()
    for method, table in (
        ("hgnc_approved", reference.hgnc_approved),
        ("hgnc_previous", reference.hgnc_previous),
        ("hgnc_alias", reference.hgnc_alias),
    ):
        ensgs = table.get(key)
        if not ensgs:
            continue
        ids = _primary(_versioned(reference, ensgs))
        if len(ensgs) == 1 and len(ids) == 1:
            return MappingResult("mapped", method, ids)
        if len(ensgs) > 1:
            return MappingResult("ambiguous", method, ids)
        return MappingResult("unresolved", f"{method}_not_in_release", [])

    return MappingResult("unresolved", "no_match", [])


def map_evidence(
    rows: list[dict[str, str]], reference: GeneReference
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, int]]:
    """Map rows whose mapping_status is 'pending'. Other rows are left alone.

    Returns the updated rows, a mapping table (one row per distinct original
    label), and status counts.
    """
    mapping_rows: dict[str, dict[str, str]] = {}
    counts: dict[str, int] = defaultdict(int)
    updated: list[dict[str, str]] = []
    for row in rows:
        row = dict(row)
        if row.get("mapping_status", "pending").strip() in ("", "pending"):
            label = row.get("original_id", "")
            result = map_label(label, reference)
            row["mapping_status"] = result.status
            row["mapping_method"] = result.method
            row["gene_id"] = result.gene_id()
            row["gene_symbol"] = (
                reference.gene_names.get(result.gene_id(), MISSING) or MISSING
            )
            if result.status == "ambiguous":
                note = "candidates: " + list_candidates(result.candidates)
                existing = row.get("notes", "")
                row["notes"] = note if existing in ("", MISSING) else f"{existing}; {note}"
            if label not in mapping_rows:
                symbols = [reference.gene_names.get(gid, "") for gid in result.candidates]
                mapping_rows[label] = {
                    "original_id": label,
                    "mapping_status": result.status,
                    "mapping_method": result.method,
                    "n_candidates": str(len(result.candidates)),
                    "candidate_gene_ids": list_candidates(result.candidates),
                    "candidate_symbols": list_candidates(
                        list(dict.fromkeys(s for s in symbols if s))
                    ),
                }
        counts[row.get("mapping_status", "")] += 1
        updated.append(row)
    return updated, list(mapping_rows.values()), dict(counts)


MAX_LISTED_CANDIDATES = 20


def list_candidates(values: list[str]) -> str:
    """Join candidates for review, capped so multi-copy names stay readable."""
    if not values:
        return MISSING
    shown = LIST_SEPARATOR.join(values[:MAX_LISTED_CANDIDATES])
    hidden = len(values) - MAX_LISTED_CANDIDATES
    return f"{shown}|...(+{hidden} more)" if hidden > 0 else shown


MAPPING_TABLE_COLUMNS = [
    "original_id",
    "mapping_status",
    "mapping_method",
    "n_candidates",
    "candidate_gene_ids",
    "candidate_symbols",
]


def command_map(args: argparse.Namespace) -> int:
    columns = SCHEMAS[f"{args.kind}_evidence"]
    header, rows = read_tsv(args.evidence)
    if header != columns:
        print(f"error: {args.evidence} header does not match the {args.kind} evidence schema",
              file=sys.stderr)
        return 1
    reference = load_reference(args.gtf, args.hgnc)
    updated, mapping_table, counts = map_evidence(rows, reference)
    write_tsv(args.output, columns, updated)
    write_tsv(args.mapping_table, MAPPING_TABLE_COLUMNS, mapping_table)
    print(json.dumps({"rows": len(updated), "mapping_status": counts}, indent=2))
    return 0


# ---------------------------------------------------------------------------
# dedupe
# ---------------------------------------------------------------------------


def dedupe_risk(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], dict[str, int]]:
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    skipped = 0
    for row in rows:
        if row.get("mapping_status") != "mapped":
            skipped += 1
            continue
        groups[row["gene_id"]].append(row)

    membership = []
    for gene_id in sorted(groups):
        members = groups[gene_id]
        membership.append(
            {
                "gene_id": gene_id,
                "gene_symbol": members[0]["gene_symbol"],
                "evidence_type": join_unique(r["evidence_type"] for r in members),
                "source": join_unique(r["source"] for r in members),
                "n_evidence": str(len(members)),
                "evidence_ids": join_unique(r["evidence_id"] for r in members),
                "notes": MISSING,
            }
        )
    stats = {
        "evidence_rows": len(rows),
        "rows_not_mapped": skipped,
        "genes": len(membership),
    }
    return membership, stats


def _publications(row: dict[str, str]) -> set[str]:
    value = row.get("publications", "")
    return {p.strip() for p in value.split(LIST_SEPARATOR) if p.strip() and p.strip() != MISSING}


def dedupe_environment(
    rows: list[dict[str, str]], min_publications: int
) -> tuple[list[dict[str, str]], dict[str, int]]:
    """One membership row per (gene, exposure, tissue, context).

    The publication threshold is applied per gene-chemical pair (gene_id,
    exposure_id), counting distinct publications across all its rows.
    """
    skipped = 0
    mapped: list[dict[str, str]] = []
    for row in rows:
        if row.get("mapping_status") != "mapped":
            skipped += 1
        else:
            mapped.append(row)

    pair_pubs: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in mapped:
        pair_pubs[(row["gene_id"], row["exposure_id"])] |= _publications(row)

    kept = [
        row
        for row in mapped
        if len(pair_pubs[(row["gene_id"], row["exposure_id"])]) >= min_publications
    ]

    groups: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in kept:
        key = (row["gene_id"], row["exposure_id"], row["tissue"], row["context"])
        groups[key].append(row)

    membership = []
    for key in sorted(groups):
        members = groups[key]
        pubs: set[str] = set()
        for r in members:
            pubs |= _publications(r)
        membership.append(
            {
                "gene_id": key[0],
                "gene_symbol": members[0]["gene_symbol"],
                "exposure": members[0]["exposure"],
                "exposure_id": key[1],
                "tissue": key[2],
                "context": key[3],
                "evidence_type": join_unique(r["evidence_type"] for r in members),
                "source": join_unique(r["source"] for r in members),
                "n_publications": str(len(pubs)),
                "n_evidence": str(len(members)),
                "evidence_ids": join_unique(r["evidence_id"] for r in members),
                "notes": MISSING,
            }
        )
    stats = {
        "evidence_rows": len(rows),
        "rows_not_mapped": skipped,
        "rows_below_publication_threshold": len(mapped) - len(kept),
        "min_publications": min_publications,
        "membership_rows": len(membership),
        "genes": len({row["gene_id"] for row in membership}),
    }
    return membership, stats


def command_dedupe(args: argparse.Namespace) -> int:
    header, rows = read_tsv(args.evidence)
    if header != SCHEMAS[f"{args.kind}_evidence"]:
        print(f"error: {args.evidence} header does not match the {args.kind} evidence schema",
              file=sys.stderr)
        return 1
    if args.kind == "risk":
        membership, stats = dedupe_risk(rows)
        columns = RISK_MEMBERSHIP_COLUMNS
    else:
        membership, stats = dedupe_environment(rows, args.min_publications)
        columns = ENVIRONMENT_MEMBERSHIP_COLUMNS
    write_tsv(args.output, columns, membership)
    print(json.dumps(stats, indent=2))
    return 0


# ---------------------------------------------------------------------------
# check
# ---------------------------------------------------------------------------


def check_table(
    kind: str,
    header: list[str],
    rows: list[dict[str, str]],
    known_gene_ids: set[str] | None = None,
) -> list[str]:
    """Return a list of human-readable problems. Empty means the table passes."""
    expected = SCHEMAS[kind]
    if header != expected:
        missing = [c for c in expected if c not in header]
        extra = [c for c in header if c not in expected]
        problem = f"header does not match the {kind} schema"
        if missing:
            problem += f"; missing {missing}"
        if extra:
            problem += f"; unexpected {extra}"
        if not missing and not extra:
            problem += "; columns are out of order"
        return [problem]

    errors: list[str] = []
    optional = OPTIONAL_COLUMNS[kind]
    is_evidence = kind.endswith("_evidence")

    def err(line: int, message: str) -> None:
        errors.append(f"line {line}: {message}")

    seen: dict[tuple[str, ...], int] = {}
    for index, row in enumerate(rows):
        line = index + 2  # header is line 1
        status = row.get("mapping_status", "mapped")
        for column in expected:
            value = (row.get(column) or "").strip()
            if value == "":
                err(line, f"{column} is blank; use {MISSING} for a missing value")
            elif value == MISSING and column not in optional:
                if is_evidence and column in ("gene_id", "gene_symbol") and status != "mapped":
                    continue
                err(line, f"{column} is required and cannot be {MISSING}")

        gene_id = row.get("gene_id", "")
        symbol = row.get("gene_symbol", "")
        if not is_evidence or status == "mapped":
            if not GENE_ID_PATTERN.match(gene_id):
                err(line, f"gene_id {gene_id!r} is not a versioned GENCODE ID")
            elif known_gene_ids is not None and gene_id not in known_gene_ids:
                err(line, f"gene_id {gene_id} is not in the reference annotation")
        if symbol.split(".")[0].startswith("ENSG"):
            err(line, f"gene_symbol {symbol!r} looks like an ID; keep symbols separate")

        if is_evidence:
            if status not in MAPPING_STATUSES:
                err(line, f"mapping_status {status!r} not in {sorted(MAPPING_STATUSES)}")
            if status == "mapped" and gene_id == MISSING:
                err(line, "mapping_status is mapped but gene_id is NA")
            if status != "mapped" and gene_id not in (MISSING, ""):
                err(line, f"gene_id must be {MISSING} when mapping_status is {status}")

        if kind == "risk_evidence":
            if row["evidence_type"] not in RISK_EVIDENCE_TYPES:
                err(line, f"evidence_type {row['evidence_type']!r} not in {sorted(RISK_EVIDENCE_TYPES)}")
            if row["evidence_level"] not in RISK_EVIDENCE_LEVELS:
                err(line, f"evidence_level {row['evidence_level']!r} not in {sorted(RISK_EVIDENCE_LEVELS)}")
            p_value = row["p_value"]
            if p_value != MISSING:
                try:
                    if not 0 <= float(p_value) <= 1:
                        raise ValueError
                except ValueError:
                    err(line, f"p_value {p_value!r} is not a number between 0 and 1")
        if kind == "environment_evidence":
            if row["evidence_type"] not in ENVIRONMENT_EVIDENCE_TYPES:
                err(line, f"evidence_type {row['evidence_type']!r} not in {sorted(ENVIRONMENT_EVIDENCE_TYPES)}")
            if row["response_type"] not in RESPONSE_TYPES:
                err(line, f"response_type {row['response_type']!r} not in {sorted(RESPONSE_TYPES)}")
        if kind == "candidates":
            flags = {}
            for column in ("is_risk", "is_environment", "is_overlap"):
                if row[column] not in ("0", "1"):
                    err(line, f"{column} must be 0 or 1")
                flags[column] = row[column] == "1"
            if flags["is_overlap"] != (flags["is_risk"] and flags["is_environment"]):
                err(line, "is_overlap must equal is_risk AND is_environment")
            if not (flags["is_risk"] or flags["is_environment"]):
                err(line, "gene is in neither set")
            if row["category"] != _category(flags["is_risk"], flags["is_environment"]):
                err(line, f"category {row['category']!r} does not match the flags")

        key = _unique_key(kind, row)
        if key in seen:
            err(line, f"duplicate {'/'.join(_unique_columns(kind))} {key} (first on line {seen[key]})")
        else:
            seen[key] = line
    return errors


def _unique_columns(kind: str) -> list[str]:
    if kind.endswith("_evidence"):
        return ["evidence_id"]
    if kind == "environment_membership":
        return ["gene_id", "exposure_id", "tissue", "context"]
    return ["gene_id"]


def _unique_key(kind: str, row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row.get(column, "") for column in _unique_columns(kind))


def command_check(args: argparse.Namespace) -> int:
    header, rows = read_tsv(args.table)
    known = None
    if args.gtf is not None:
        reference = GeneReference()
        load_gtf_genes(args.gtf, reference)
        known = set(reference.gene_names)
    errors = check_table(args.kind, header, rows, known)
    if errors:
        for problem in errors:
            print(f"{args.table}: {problem}", file=sys.stderr)
        print(f"FAIL: {len(errors)} problem(s) in {len(rows)} row(s)", file=sys.stderr)
        return 1
    print(f"PASS: {args.table} ({args.kind}, {len(rows)} row(s))")
    return 0


# ---------------------------------------------------------------------------
# overlap
# ---------------------------------------------------------------------------


def _category(is_risk: bool, is_environment: bool) -> str:
    if is_risk and is_environment:
        return "overlap"
    return "risk_only" if is_risk else "environment_only"


def compute_overlap(
    risk: list[dict[str, str]], environment: list[dict[str, str]]
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, int]]:
    risk_by_gene = {row["gene_id"]: row for row in risk}
    env_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in environment:
        env_by_gene[row["gene_id"]].append(row)

    overlap_rows = []
    candidate_rows = []
    for gene_id in sorted(set(risk_by_gene) | set(env_by_gene)):
        in_risk = gene_id in risk_by_gene
        in_env = gene_id in env_by_gene
        symbol = (risk_by_gene[gene_id] if in_risk else env_by_gene[gene_id][0])["gene_symbol"]
        if in_risk and in_env:
            overlap_rows.append(
                {
                    "gene_id": gene_id,
                    "gene_symbol": symbol,
                    "risk_source": risk_by_gene[gene_id]["source"],
                    "environment_source": join_unique(
                        part
                        for r in env_by_gene[gene_id]
                        for part in r["source"].split(LIST_SEPARATOR)
                    ),
                }
            )
        candidate_rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": symbol,
                "is_risk": "1" if in_risk else "0",
                "is_environment": "1" if in_env else "0",
                "is_overlap": "1" if in_risk and in_env else "0",
                "category": _category(in_risk, in_env),
            }
        )
    counts = {category: 0 for category in ("risk_only", "environment_only", "overlap")}
    for row in candidate_rows:
        counts[row["category"]] += 1
    counts["union"] = len(candidate_rows)
    counts["risk_total"] = len(risk_by_gene)
    counts["environment_total"] = len(env_by_gene)
    return overlap_rows, candidate_rows, counts


def command_overlap(args: argparse.Namespace) -> int:
    problems = []
    loaded = {}
    for kind, path in (("risk_membership", args.risk), ("environment_membership", args.environment)):
        header, rows = read_tsv(path)
        errors = check_table(kind, header, rows)
        problems.extend(f"{path}: {e}" for e in errors)
        loaded[kind] = rows
    if problems:
        for problem in problems:
            print(problem, file=sys.stderr)
        print("FAIL: fix the membership files before computing the overlap", file=sys.stderr)
        return 1
    overlap_rows, candidate_rows, counts = compute_overlap(
        loaded["risk_membership"], loaded["environment_membership"]
    )
    write_tsv(args.overlap, OVERLAP_COLUMNS, overlap_rows)
    write_tsv(args.candidates, CANDIDATE_COLUMNS, candidate_rows)
    print(json.dumps(counts, indent=2))
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("map", help="map original labels to versioned GENCODE gene IDs")
    p.add_argument("--kind", choices=["risk", "environment"], required=True)
    p.add_argument("--evidence", type=Path, required=True, help="evidence TSV to map")
    p.add_argument("--gtf", type=Path, required=True, help="GENCODE v50 GTF (.gtf or .gtf.gz)")
    p.add_argument("--hgnc", type=Path, help="HGNC complete-set TSV for previous/alias symbols")
    p.add_argument("--output", type=Path, required=True, help="mapped evidence TSV")
    p.add_argument("--mapping-table", type=Path, required=True,
                   help="one row per original label with candidates, for review")
    p.set_defaults(func=command_map)

    p = sub.add_parser("dedupe", help="build one membership row per gene")
    p.add_argument("--kind", choices=["risk", "environment"], required=True)
    p.add_argument("--evidence", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--min-publications", type=int, default=1,
                   help="environment only: distinct publications required per gene-chemical pair")
    p.set_defaults(func=command_dedupe)

    p = sub.add_parser("check", help="schema-check a gene-set table")
    p.add_argument("--kind", choices=sorted(SCHEMAS), required=True)
    p.add_argument("--table", type=Path, required=True)
    p.add_argument("--gtf", type=Path, help="also require gene_ids to exist in this annotation")
    p.set_defaults(func=command_check)

    p = sub.add_parser("overlap", help="build overlap.tsv and candidate_genes.tsv")
    p.add_argument("--risk", type=Path, required=True, help="endometriosis_risk.tsv")
    p.add_argument("--environment", type=Path, required=True, help="environment_response.tsv")
    p.add_argument("--overlap", type=Path, required=True, help="output overlap.tsv")
    p.add_argument("--candidates", type=Path, required=True, help="output candidate_genes.tsv")
    p.set_defaults(func=command_overlap)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
