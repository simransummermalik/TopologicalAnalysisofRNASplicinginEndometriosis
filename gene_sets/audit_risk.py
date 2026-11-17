#!/usr/bin/env python3
"""Audit the genetic-risk evidence and membership tables (G9-G10).

Written separately from extract_gwas.py and the review step: it re-reads the
cached GWAS Catalog responses, the Koller supplement, the review table, and
GENCODE, and checks the saved tables against them.

Checks:

    A  schema        both tables pass ``check``; membership unique by gene_id
    B  provenance    every Catalog row matches its cached association record;
                     every Supplementary Table 4 row matches its sheet row
    C  rederive      the extraction-time include/exclude outcome of every
                     Catalog row and the set of ST4 rows are recomputed from G1
    D  mapping       mapped IDs exist; symbols match GENCODE; automatic
                     mappings match their label; manual ones match the review
    E  review        no row is open; every exclusion has a reason; every
                     review decision was applied
    F  dedupe        membership is rebuilt from the evidence and must match
    G  evidence type only GWAS associations from the approved sources, all
                     p < 5e-8, no expression-based trait among included studies
    H  proximity     each mapped gene lies near its lead variant (GRCh38 for
                     Catalog rows; GRCh37 via GENCODE v19 for ST4 rows)
    I  sample        writes a fixed-seed sample of genes for a person to read

Prints a JSON report; exits 1 if any automatic check A-H fails.
"""

from __future__ import annotations

import argparse
import gzip
import json
import random
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import extract_gwas  # noqa: E402  (only its .xlsx reader)
import gene_sets as gs  # noqa: E402

THRESHOLD = 5e-8
MAX_EXAMPLES = 5
NEAR_BP = 1_000_000
APPROVED_SOURCES = {"GWAS Catalog", "Koller 2026 Supplementary Table 4"}
EXPRESSION_WORDS = re.compile(r"expression|eQTL|transcript|methylation|protein level", re.I)


def gene_spans(path: Path) -> dict[str, tuple[str, int, int]]:
    """Unversioned gene ID -> (chromosome without 'chr', start, end)."""
    spans = {}
    with gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene":
                gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
                if gid.endswith("_PAR_Y"):
                    continue
                spans[gid.split(".")[0]] = (f[0].removeprefix("chr"), int(f[3]), int(f[4]))
    return spans


def distance(span, chrom: str, pos: int):
    if span is None:
        return None
    c, s, e = span
    if c != chrom:
        return float("inf")
    return 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--membership", type=Path, required=True)
    parser.add_argument("--review", type=Path, required=True)
    parser.add_argument("--catalog-dir", type=Path, required=True)
    parser.add_argument("--koller-xlsx", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True, help="GENCODE v50 (GRCh38)")
    parser.add_argument("--gtf-v19", type=Path, required=True, help="GENCODE v19 (GRCh37)")
    parser.add_argument("--sample-size", type=int, default=15)
    parser.add_argument("--seed", type=int, default=20261117)
    parser.add_argument("--sample-output", type=Path, required=True)
    args = parser.parse_args(argv)

    results: dict[str, dict] = {}

    def record(name, problems, **facts):
        results[name] = {"pass": not problems, "problems": len(problems),
                         "examples": problems[:MAX_EXAMPLES], **facts}

    header, rows = gs.read_tsv(args.evidence)
    m_header, membership = gs.read_tsv(args.membership)
    reference = gs.GeneReference()
    gs.load_gtf_genes(args.gtf, reference)
    known = set(reference.gene_names)

    # A. schema and uniqueness
    problems = gs.check_table("risk_evidence", header, rows, known)
    problems += gs.check_table("risk_membership", m_header, membership, known)
    ids = Counter(m["gene_id"] for m in membership)
    problems += [f"gene_id {g} appears {n} times in membership" for g, n in ids.items() if n > 1]
    record("A_schema", problems, evidence_rows=len(rows), membership_rows=len(membership))

    # B. provenance
    catalog = {}
    for page in sorted(args.catalog_dir.glob("trait_p*.json")):
        for a in json.loads(page.read_text())["_embedded"]["associations"]:
            catalog[str(a["association_id"])] = a
    studies = {p.stem.removeprefix("study_"): json.loads(p.read_text())
               for p in args.catalog_dir.glob("study_*.json")}
    sheet = extract_gwas.read_xlsx_sheet(args.koller_xlsx, extract_gwas.ST4_SHEET)
    problems = []
    catalog_rows = [r for r in rows if r["source"] == "GWAS Catalog"]
    st4_rows = [r for r in rows if r["source"] != "GWAS Catalog"]
    for r in catalog_rows:
        aid = r["evidence_id"].removeprefix("GWASCAT:")
        a = catalog.get(aid)
        if a is None:
            problems.append(f"{r['evidence_id']}: not in the cached Catalog pages")
            continue
        expect = {
            "source_record": f"{a['accession_id']} association {aid}",
            "variant_id": "|".join(s["rs_id"] for s in a["snp_allele"]) or "NA",
            "original_id": ", ".join(a["mapped_genes"]) or "NA",
        }
        for column, value in expect.items():
            if r[column] != value:
                problems.append(f"{r['evidence_id']}: {column} {r[column]!r} != Catalog {value!r}")
        if float(r["p_value"]) != float(a["p_value"]):
            problems.append(f"{r['evidence_id']}: p_value differs from the Catalog")
    if len({r["evidence_id"] for r in catalog_rows}) != len(catalog):
        problems.append(f"{len(catalog)} cached associations but {len(catalog_rows)} Catalog rows")
    by_rs = {}
    for number, row in enumerate(sheet, start=1):
        if row and row[0].startswith(extract_gwas.ST4_SECTION_END):
            break
        if row and row[0].startswith("rs"):
            by_rs[row[0]] = (number, row)
    for r in st4_rows:
        rs = r["evidence_id"].removeprefix("KOLLER2026_ST4:").split("/")[0]
        number, row = by_rs.get(rs, (None, None))
        if row is None:
            problems.append(f"{r['evidence_id']}: lead SNP not in the first section of the sheet")
            continue
        if r["source_record"] != f"Supplementary Table 4 sheet row {number}":
            problems.append(f"{r['evidence_id']}: source_record does not give sheet row {number}")
        if r["original_id"] != row[5].strip() or float(r["p_value"]) != float(row[9]):
            problems.append(f"{r['evidence_id']}: gene label or p-value differs from the sheet")
    record("B_provenance", problems, catalog_rows=len(catalog_rows), st4_rows=len(st4_rows),
           st4_loci=len({r["variant_id"] for r in st4_rows}))

    # C. rederive the extraction-time outcome from G1 and the 2026-11-17 decisions
    problems = []
    for r in catalog_rows:
        a = catalog.get(r["evidence_id"].removeprefix("GWASCAT:"))
        if a is None:
            continue
        traits = {e["efo_trait"] for e in a["efo_traits"]}
        should_include = (a["p_value"] < THRESHOLD and traits == {"endometriosis"}
                          and not any("MTAG" in t for t in a["reported_trait"])
                          and not (a["pubmed_id"] == "42056605" and a["accession_id"] != "GCST90841381"))
        excluded_at_extraction = r["mapping_status"] == "excluded" and "manual review" not in r["notes"]
        if should_include == excluded_at_extraction:
            problems.append(f"{r['evidence_id']}: extraction outcome disagrees with G1")
    if set(by_rs) != {r["variant_id"] for r in st4_rows}:
        problems.append("Supplementary Table 4 lead SNPs differ from the sheet's first section")
    record("C_rederive", problems, st4_sheet_loci=len(by_rs))

    # D. mapping
    _, review_rows = gs.read_tsv(args.review)
    review = {r["evidence_id"]: r for r in review_rows}
    problems = []
    methods = Counter()
    for r in rows:
        if r["mapping_status"] != "mapped":
            continue
        methods[r["mapping_method"]] += 1
        if r["gene_symbol"] != gs.symbol_for(reference, r["gene_id"]):
            problems.append(f"{r['evidence_id']}: symbol {r['gene_symbol']} != GENCODE")
        if r["mapping_method"] == "gencode_symbol" and reference.gene_names[r["gene_id"]] != r["original_id"]:
            problems.append(f"{r['evidence_id']}: gencode_symbol mapping but label differs")
        if r["mapping_method"] == "manual":
            key = r["evidence_id"].split("/")[0]
            decision = review.get(key)
            if decision is None or r["gene_id"] not in decision["gene_id"].split("|"):
                problems.append(f"{r['evidence_id']}: manual gene not in the review table")
    record("D_mapping", problems, mapped_rows=sum(methods.values()), by_method=dict(methods))

    # E. review completeness
    problems = [f"{r['evidence_id']}: still {r['mapping_status']}" for r in rows
                if r["mapping_status"] in ("ambiguous", "unresolved", "pending")]
    problems += [f"{r['evidence_id']}: excluded without a reason" for r in rows
                 if r["mapping_status"] == "excluded" and "excluded" not in r["notes"]]
    row_ids = {r["evidence_id"].split("/")[0] for r in rows if "manual review" in r["notes"]}
    problems += [f"review decision {k} not applied" for k in review if k not in row_ids]
    record("E_review", problems, review_decisions=len(review),
           status=dict(Counter(r["mapping_status"] for r in rows)))

    # F. rebuild membership
    groups = defaultdict(set)
    for r in rows:
        if r["mapping_status"] == "mapped":
            groups[r["gene_id"]].add(r["evidence_id"])
    problems = []
    seen = set()
    for m in membership:
        seen.add(m["gene_id"])
        if set(m["evidence_ids"].split("|")) != groups.get(m["gene_id"], set()):
            problems.append(f"{m['gene_id']}: evidence_ids differ from the evidence table")
    problems += [f"{g}: mapped but missing from membership" for g in set(groups) - seen]
    record("F_dedupe", problems, genes=len(membership), mapped_rows=sum(len(v) for v in groups.values()))

    # G. evidence type: no expression-only or off-source rows among the included
    problems = []
    included = [r for r in rows if r["mapping_status"] == "mapped"]
    for r in included:
        if r["evidence_type"] != "gwas_association" or r["source"] not in APPROVED_SOURCES:
            problems.append(f"{r['evidence_id']}: not a GWAS association from an approved source")
        if float(r["p_value"]) >= THRESHOLD:
            problems.append(f"{r['evidence_id']}: p >= 5e-8")
    included_studies = {r["source_record"].split()[0] for r in included if r["source"] == "GWAS Catalog"}
    traits = {acc: studies[acc].get("disease_trait", "") for acc in included_studies}
    problems += [f"{acc}: trait {t!r} looks expression-based" for acc, t in traits.items()
                 if EXPRESSION_WORDS.search(t)]
    record("G_evidence_type", problems, included_study_traits=dict(Counter(traits.values())))

    # H. proximity of each mapped gene to its lead variant
    v50 = gene_spans(args.gtf)
    v19 = gene_spans(args.gtf_v19)
    problems, inside, near, unknown = [], 0, 0, []
    for r in included:
        uid = r["gene_id"].split(".")[0]
        if r["source"] == "GWAS Catalog":
            locs = re.search(r"location (\S+?);", r["notes"]).group(1).split("|")
            spans, positions = v50.get(uid), [loc.split(":") for loc in locs if ":" in loc]
        else:
            m = re.search(r"chr(\w+):(\d+)", r["notes"])
            spans, positions = v19.get(uid), [m.groups()]
        dists = [distance(spans, c, int(p)) for c, p in positions]
        dists = [d for d in dists if d is not None]
        if not dists:
            unknown.append(r["evidence_id"])
            continue
        best = min(dists)
        if best == 0:
            inside += 1
        elif best <= NEAR_BP:
            near += 1
        else:
            problems.append(f"{r['evidence_id']}: {r['gene_symbol']} is {best} bp from its lead variant")
    record("H_proximity", problems, inside_gene=inside, within_1Mb=near,
           no_coordinates=len(unknown), no_coordinates_ids=unknown[:10])

    # I. sample of genes for a person to read
    rng = random.Random(args.seed)
    sample = rng.sample(sorted(membership, key=lambda m: m["gene_id"]),
                        min(args.sample_size, len(membership)))
    gs.write_tsv(args.sample_output,
                 ["gene_id", "gene_symbol", "source", "n_evidence", "notes", "check_result", "check_note"],
                 [{**{k: m[k] for k in ("gene_id", "gene_symbol", "source", "n_evidence", "notes")},
                   "check_result": "unchecked", "check_note": gs.MISSING} for m in sample])
    results["I_sample_written"] = {"pass": True, "seed": args.seed, "size": len(sample)}

    passed = all(r["pass"] for r in results.values())
    print(json.dumps({"pass": passed, **results}, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
