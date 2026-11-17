#!/usr/bin/env python3
"""Audit the environmental evidence and membership tables (E8).

Written separately from extract_ctd.py and review_ctd_context.py on purpose:
it re-reads the raw CTD file and the reference files and checks the saved
tables against them. Only the generic TSV reader and schema check are shared.

Checks:

    A  schema        both tables pass ``check`` against the GENCODE annotation
    B  provenance    every evidence row matches its raw CTD line exactly
    C  rederive      an independent pass over the raw file selects the same lines
    D  mapping       mapped gene IDs agree with CTD GeneID -> HGNC -> Ensembl
    E  exposure      exposure IDs/names are on the frozen list; context labels
                     agree with CTD's actions; exclusions carry a reason
    F  dedupe        membership is rebuilt from the evidence and must match
    G  sample        writes a fixed-seed sample of rows for checking against
                     the cited papers (done by a person)

Prints a JSON report; exits 1 if any automatic check fails.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gene_sets as gs  # noqa: E402

MAX_EXAMPLES = 5


class Audit:
    def __init__(self) -> None:
        self.results: dict[str, dict] = {}

    def record(self, name: str, problems: list[str], **facts) -> None:
        self.results[name] = {
            "pass": not problems,
            "problems": len(problems),
            "examples": problems[:MAX_EXAMPLES],
            **facts,
        }

    @property
    def passed(self) -> bool:
        return all(r["pass"] for r in self.results.values())


def raw_ctd_lines(path: Path, wanted: set[int]) -> dict[int, list[str]]:
    lines = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if number in wanted:
                lines[number] = line.rstrip("\n").split("\t")
    return lines


def rederive_lines(path: Path, ids: set[str]) -> set[int]:
    """Independent selection: listed chemical, human, the row's own gene measured."""
    kept = set()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[1] not in ids or f[6] != "Homo sapiens":
                continue
            kinds = {a.split("^")[1] for a in f[9].split("|") if "^" in a}
            if not kinds & {"expression", "activity"}:
                continue
            symbol = f[3]
            # The row's gene must directly follow "expression of"/"activity of",
            # allowing CTD's "expression of and affects the X of GENE" phrasing.
            for match in re.finditer(r"(?:expression|activity) of ((?:and [a-z ]+? of )*)(\S+)", f[8]):
                if match.group(2).rstrip("]") == symbol:
                    kept.add(number)
                    break
    return kept


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--membership", type=Path, required=True)
    parser.add_argument("--chemicals", type=Path, required=True)
    parser.add_argument("--ctd", type=Path, required=True, help="raw CTD_chem_gene_ixns.tsv.gz")
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--hgnc", type=Path, required=True)
    parser.add_argument("--min-publications", type=int, required=True)
    parser.add_argument("--sample-size", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20261117)
    parser.add_argument("--sample-output", type=Path, required=True)
    args = parser.parse_args(argv)

    audit = Audit()
    header, rows = gs.read_tsv(args.evidence)
    m_header, membership = gs.read_tsv(args.membership)
    _, chem_rows = gs.read_tsv(args.chemicals)
    chemicals = {r["ctd_chemical_id"]: r["ctd_name"] for r in chem_rows}
    reference = gs.GeneReference()
    gs.load_gtf_genes(args.gtf, reference)
    known = set(reference.gene_names)

    # A. schema
    problems = gs.check_table("environment_evidence", header, rows, known)
    problems += gs.check_table("environment_membership", m_header, membership, known)
    audit.record("A_schema", problems, evidence_rows=len(rows), membership_rows=len(membership))

    # B. provenance against the raw file
    line_of = {}
    problems = []
    for r in rows:
        m = re.fullmatch(r"CTD_chem_gene_ixns\.tsv\.gz line (\d+)", r["source_record"])
        if not m or r["source"] != "CTD":
            problems.append(f"{r['evidence_id']}: source/source_record not a CTD line")
            continue
        line_of[r["evidence_id"]] = int(m.group(1))
    versions = Counter(r["source_version"] for r in rows)
    if len(versions) != 1:
        problems.append(f"more than one source_version: {dict(versions)}")
    raw = raw_ctd_lines(args.ctd, set(line_of.values()))
    for r in rows:
        number = line_of.get(r["evidence_id"])
        if number is None:
            continue
        f = raw.get(number)
        if f is None or len(f) != 11:
            problems.append(f"{r['evidence_id']}: raw line {number} missing or malformed")
            continue
        expected = {
            "original_id": f[3],
            "exposure": f[0],
            "exposure_id": f"MESH:{f[1]}",
            "species": f[6],
            "publications": "|".join(f"PMID:{p}" for p in f[10].split("|") if p),
        }
        for column, value in expected.items():
            if r[column] != value:
                problems.append(f"{r['evidence_id']}: {column} {r[column]!r} != raw {value!r}")
        if f"interaction: {f[8]};" not in r["notes"] or f"CTD GeneID {f[4]};" not in r["notes"]:
            problems.append(f"{r['evidence_id']}: note does not quote raw interaction/GeneID")
        if not r["evidence_id"].startswith(f"CTD:{number}:"):
            problems.append(f"{r['evidence_id']}: evidence_id does not name line {number}")
    audit.record("B_provenance", problems, rows_traced=len(raw), source_version=next(iter(versions)))

    # C. independent re-derivation of which lines are in the table
    derived = rederive_lines(args.ctd, {c.removeprefix("MESH:") for c in chemicals})
    table_lines = set(line_of.values())
    problems = [f"line {n} rederived but not in table" for n in sorted(derived - table_lines)]
    problems += [f"line {n} in table but not rederived" for n in sorted(table_lines - derived)]
    audit.record("C_rederive", problems, rederived_lines=len(derived), table_lines=len(table_lines))

    # D. mapping cross-check through NCBI GeneID
    entrez_to_ensembl = {}
    with open(args.hgnc, newline="", encoding="utf-8") as handle:
        for h in csv.DictReader(handle, delimiter="\t"):
            if h.get("entrez_id") and h.get("ensembl_gene_id"):
                entrez_to_ensembl[h["entrez_id"]] = h["ensembl_gene_id"]
    problems, agree, unchecked = [], 0, []
    for r in rows:
        if r["mapping_status"] != "mapped":
            continue
        gene_id_ncbi = re.search(r"CTD GeneID (\S+?);", r["notes"]).group(1)
        ensembl = entrez_to_ensembl.get(gene_id_ncbi)
        if ensembl is None:
            unchecked.append(r["original_id"])
        elif r["gene_id"].split(".")[0] == ensembl:
            agree += 1
        else:
            problems.append(f"{r['evidence_id']}: {r['original_id']} -> {r['gene_id']}, NCBI says {ensembl}")
        if reference.gene_names.get(r["gene_id"]) != r["gene_symbol"]:
            problems.append(f"{r['evidence_id']}: gene_symbol {r['gene_symbol']} != GENCODE name")
    audit.record("D_mapping", problems, agree=agree, uncheckable=len(unchecked),
                 uncheckable_labels=sorted(set(unchecked)))

    # E. exposure, context, exclusions
    problems = []
    for r in rows:
        if chemicals.get(r["exposure_id"]) is None:
            problems.append(f"{r['evidence_id']}: exposure_id {r['exposure_id']} not on the frozen list")
        f = raw.get(line_of.get(r["evidence_id"], -1))
        if f is None:
            continue
        kinds = {a.split("^")[1] for a in f[9].split("|") if "^" in a}
        labels = set(r["context"].split("|"))
        if ("cotreatment" in kinds) != ("cotreatment" in labels):
            problems.append(f"{r['evidence_id']}: cotreatment label disagrees with CTD actions")
        has_reaction_label = bool(labels & {"dioxin_effect_modified", "dioxin_modifies_other"})
        if ("reaction" in kinds) != has_reaction_label:
            problems.append(f"{r['evidence_id']}: reaction label disagrees with CTD actions")
        if "cotreatment" in labels and r["mapping_status"] != "excluded":
            problems.append(f"{r['evidence_id']}: co-treatment row is not excluded (E6 rule)")
        if r["mapping_status"] == "excluded" and "excluded" not in r["notes"]:
            problems.append(f"{r['evidence_id']}: excluded without a recorded reason")
        if not r["tissue"].startswith("cited_paper:"):
            problems.append(f"{r['evidence_id']}: tissue {r['tissue']!r} lacks the paper-level tag")
    status = Counter(r["mapping_status"] for r in rows)
    audit.record("E_exposure_context", problems, mapping_status=dict(status))

    # F. rebuild membership independently and compare
    mapped = [r for r in rows if r["mapping_status"] == "mapped"]
    pair_pubs = defaultdict(set)
    for r in mapped:
        pair_pubs[(r["gene_id"], r["exposure_id"])].update(r["publications"].split("|"))
    expected = defaultdict(lambda: {"ids": set(), "pubs": set()})
    for r in mapped:
        if len(pair_pubs[(r["gene_id"], r["exposure_id"])]) >= args.min_publications:
            key = (r["gene_id"], r["exposure_id"], r["tissue"], r["context"])
            expected[key]["ids"].add(r["evidence_id"])
            expected[key]["pubs"].update(r["publications"].split("|"))
    problems = []
    seen_ids = Counter()
    actual_keys = set()
    for m in membership:
        key = (m["gene_id"], m["exposure_id"], m["tissue"], m["context"])
        actual_keys.add(key)
        ids = set(m["evidence_ids"].split("|"))
        seen_ids.update(ids)
        exp = expected.get(key)
        if exp is None:
            problems.append(f"membership row {key} not expected")
            continue
        if ids != exp["ids"]:
            problems.append(f"{key}: evidence_ids differ")
        if m["n_publications"] != str(len(exp["pubs"])) or m["n_evidence"] != str(len(ids)):
            problems.append(f"{key}: counts differ")
    problems += [f"expected membership row {k} missing" for k in sorted(set(expected) - actual_keys)]
    problems += [f"evidence {i} in {n} membership rows" for i, n in seen_ids.items() if n > 1]
    audit.record("F_dedupe", problems, membership_rows=len(membership),
                 genes=len({m["gene_id"] for m in membership}),
                 gene_chemical_pairs=len({(m["gene_id"], m["exposure_id"]) for m in membership}))

    # G. fixed-seed sample for a person to check against the papers
    rng = random.Random(args.seed)
    in_set = {i for m in membership for i in m["evidence_ids"].split("|")}
    pool = sorted(r["evidence_id"] for r in mapped if r["evidence_id"] in in_set)
    by_id = {r["evidence_id"]: r for r in rows}
    sample = [by_id[i] for i in rng.sample(pool, min(args.sample_size, len(pool)))]
    columns = ["evidence_id", "gene_symbol", "exposure", "context", "tissue", "publications",
               "interaction", "check_result", "check_note"]
    gs.write_tsv(args.sample_output, columns, [
        {**{c: r[c] for c in columns[:6]},
         "interaction": re.search(r"interaction: (.*?); CTD does not", r["notes"]).group(1),
         "check_result": "unchecked", "check_note": gs.MISSING}
        for r in sample
    ])
    audit.record("G_sample_written", [], seed=args.seed, size=len(sample), pool=len(pool),
                 path=str(args.sample_output))

    print(json.dumps({"pass": audit.passed, **audit.results}, indent=2))
    return 0 if audit.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
