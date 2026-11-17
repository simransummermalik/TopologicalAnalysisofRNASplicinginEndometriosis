#!/usr/bin/env python3
"""Draft genetic-risk evidence rows from the two approved sources (G4-G5).

Sources (docs/gene_set_methods.md, G2):

* GWAS Catalog, endometriosis trait. Every association is written. Rows that
  fail a G1 rule stay in the table as ``excluded`` with the reason, checked
  in this order:
    1. p >= 5e-8
    2. cross-trait (pleiotropy) analysis: an EFO trait other than endometriosis
    3. MTAG analysis
    4. a Koller et al. 2026 analysis other than the cross-ancestry combined
       definition (GCST90841381)
  Subtype analyses (for example stage 3 or 4) are kept, decision 2026-11-17.
  ``original_id`` is the Catalog's mapped genes, joined with ", " so that a
  two-gene locus is marked ambiguous by ``map``. Author-reported genes go in
  the notes for the G7 review.

* Koller et al. 2026 Supplementary Table 4, first section only (the 80
  cross-ancestry loci). ``original_id`` is the authors' "Nearest gene" exactly
  as printed, including ":"-joined lists and retired clone names.

Nothing is mapped or guessed here; included rows are ``pending`` for ``map``.

Usage:

    python3 gene_sets/extract_gwas.py \\
      --catalog-dir /path/to/gwas_cache --fetch \\
      --koller-xlsx /path/to/41588_2026_2582_MOESM4_ESM.xlsx \\
      --output data/gene_sets/risk_evidence.tsv
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gene_sets as gs  # noqa: E402

THRESHOLD = 5e-8
API = "https://www.ebi.ac.uk/gwas/rest/api/v2"
STATS_URL = "https://www.ebi.ac.uk/gwas/api/search/stats"
KOLLER_PMID = "42056605"
KOLLER_HEADLINE = "GCST90841381"
KOLLER_SOURCE = "Koller 2026 Supplementary Table 4"
KOLLER_VERSION = (
    "Nat Genet 58(5):1051-1061 (2026), doi:10.1038/s41588-026-02582-2; "
    "41588_2026_2582_MOESM4_ESM.xlsx"
)
KOLLER_POPULATION = (
    "cross-ancestry (AFR, AMR, EAS, EUR, MID, SAS) endometriosis, combined "
    "definition; 105,869 cases, 1,282,731 controls"
)
ST4_SHEET = "Supplementary Table 4"
ST4_SECTION_END = "Novel hits in other than cross-ancestry GWAS"

EXCLUDE_THRESHOLD = "excluded: p >= 5e-8 (G1 threshold)"
EXCLUDE_CROSS_TRAIT = (
    "excluded: cross-trait (pleiotropy) analysis with {traits}, not an "
    "endometriosis association alone (decision 2026-11-17)"
)
EXCLUDE_MTAG = (
    "excluded: MTAG analysis borrows power from correlated traits, not a plain "
    "GWAS or meta-analysis (decision 2026-11-17)"
)
EXCLUDE_KOLLER = (
    "excluded: Koller et al. 2026 analysis other than the cross-ancestry "
    "combined definition (decisions 2026-11-12 and 2026-11-17)"
)


# ---------------------------------------------------------------------------
# GWAS Catalog
# ---------------------------------------------------------------------------


def _get_json(url: str, path: Path, fetch: bool):
    if not path.exists():
        if not fetch:
            raise ValueError(f"{path} is missing; rerun with --fetch to download it")
        last = None
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=400) as response:
                    path.write_bytes(response.read())
                break
            except OSError as error:
                last = error
                time.sleep(5 * (attempt + 1))
        else:
            raise ValueError(f"could not download {url}: {last}")
        time.sleep(0.3)
    return json.loads(path.read_text(encoding="utf-8"))


def load_catalog(cache: Path, fetch: bool) -> tuple[list[dict], dict[str, dict], dict]:
    cache.mkdir(parents=True, exist_ok=True)
    stats = _get_json(STATS_URL, cache / "stats.json", fetch)
    associations, page = [], 0
    while True:
        data = _get_json(
            f"{API}/associations?efo_trait=endometriosis&size=200&page={page}",
            cache / f"trait_p{page}.json", fetch,
        )
        associations += data["_embedded"]["associations"]
        if page >= data["page"]["totalPages"] - 1:
            break
        page += 1
    if len({a["association_id"] for a in associations}) != len(associations):
        raise ValueError("GWAS Catalog pages returned a repeated association")
    studies = {
        acc: _get_json(f"{API}/studies/{acc}", cache / f"study_{acc}.json", fetch)
        for acc in sorted({a["accession_id"] for a in associations})
    }
    return associations, studies, stats


def exclusion_reason(assoc: dict) -> str | None:
    if assoc["p_value"] >= THRESHOLD:
        return EXCLUDE_THRESHOLD
    others = sorted({e["efo_trait"] for e in assoc["efo_traits"]} - {"endometriosis"})
    if others:
        return EXCLUDE_CROSS_TRAIT.format(traits=", ".join(others))
    if any("MTAG" in t for t in assoc["reported_trait"]):
        return EXCLUDE_MTAG
    if assoc["pubmed_id"] == KOLLER_PMID and assoc["accession_id"] != KOLLER_HEADLINE:
        return EXCLUDE_KOLLER
    return None


def author_reported_genes(assoc: dict, cache: Path, fetch: bool) -> list[str]:
    data = _get_json(f"{API}/associations/{assoc['association_id']}/loci",
                     cache / f"loci_{assoc['association_id']}.json", fetch)
    genes = []
    for locus in data.get("_embedded", {}).get("loci", []):
        for gene in locus.get("author_reported_genes") or []:
            name = gene.get("gene_name") if isinstance(gene, dict) else str(gene)
            if name and name not in genes:
                genes.append(name)
    return genes


def catalog_rows(associations, studies, stats, cache: Path, fetch: bool):
    version = (f"GWAS Catalog release {stats['date']} (Ensembl {stats['ensemblbuild']}, "
               f"dbSNP {stats['dbsnpbuild']})")
    rows, reasons = [], Counter()
    for assoc in sorted(associations, key=lambda a: a["association_id"]):
        study = studies[assoc["accession_id"]]
        reason = exclusion_reason(assoc)
        reasons[reason.split(":")[1].split("(")[0].strip() if reason else "included"] += 1
        mapped = assoc.get("mapped_genes") or []
        notes = [f"PMID {assoc['pubmed_id']} ({assoc.get('first_author', 'NA')})",
                 f"location {'|'.join(assoc.get('locations') or []) or 'NA'}",
                 f"Catalog mapped genes: {', '.join(mapped) or 'none'}"]
        if reason is None:
            reported = author_reported_genes(assoc, cache, fetch)
            notes.append(f"author-reported genes: {', '.join(reported) or 'none'}")
        else:
            notes.append(reason)
        rows.append({
            "evidence_id": f"GWASCAT:{assoc['association_id']}",
            "source": "GWAS Catalog",
            "source_version": version,
            "source_record": f"{assoc['accession_id']} association {assoc['association_id']}",
            "original_id": ", ".join(mapped) if mapped else gs.MISSING,
            "variant_id": "|".join(s["rs_id"] for s in assoc.get("snp_allele") or []) or gs.MISSING,
            "p_value": repr(float(assoc["p_value"])),
            "evidence_type": "gwas_association",
            "evidence_level": "variant",
            "population": f"{study.get('disease_trait', 'NA')}; {study.get('initial_sample_size', 'NA')}",
            "gene_id": gs.MISSING,
            "gene_symbol": gs.MISSING,
            "mapping_status": "pending" if reason is None else "excluded",
            "mapping_method": "not_run",
            "notes": "; ".join(notes),
        })
    return rows, reasons


# ---------------------------------------------------------------------------
# Koller Supplementary Table 4 (read with the standard library only)
# ---------------------------------------------------------------------------

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def read_xlsx_sheet(path: Path, sheet_name: str) -> list[list[str]]:
    with zipfile.ZipFile(path) as z:
        workbook = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        targets = {r.get("Id"): r.get("Target") for r in rels}
        target = None
        for sheet in workbook.iter(_M + "sheet"):
            if sheet.get("name") == sheet_name:
                target = targets[sheet.get(_R + "id")]
        if target is None:
            raise ValueError(f"{path}: no sheet named {sheet_name!r}")
        target = target.lstrip("/").removeprefix("xl/")
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            for item in ET.fromstring(z.read("xl/sharedStrings.xml")).iter(_M + "si"):
                shared.append("".join(t.text or "" for t in item.iter(_M + "t")))
        rows = []
        for row in ET.fromstring(z.read("xl/" + target)).iter(_M + "row"):
            cells = {}
            for cell in row.iter(_M + "c"):
                column = 0
                for ch in re.match(r"[A-Z]+", cell.get("r")).group(0):
                    column = column * 26 + ord(ch) - 64
                value = cell.find(_M + "v")
                if cell.get("t") == "inlineStr":
                    text = "".join(t.text or "" for t in cell.iter(_M + "t"))
                elif value is None:
                    text = ""
                elif cell.get("t") == "s":
                    text = shared[int(value.text)]
                else:
                    text = value.text or ""
                cells[column - 1] = text
            rows.append([cells.get(i, "") for i in range(max(cells) + 1)] if cells else [])
        return rows


def koller_rows(sheet: list[list[str]]):
    header_at = next((i for i, r in enumerate(sheet) if r and r[0] == "Lead SNP"), None)
    if header_at is None:
        raise ValueError("Supplementary Table 4: no 'Lead SNP' header row")
    header, sub = sheet[header_at], sheet[header_at + 2]
    expected = {0: "Lead SNP", 1: "Chromosome", 2: "Position", 5: "Nearest gene", 6: "Type"}
    for index, name in expected.items():
        if header[index].strip() != name:
            raise ValueError(f"Supplementary Table 4: column {index} is {header[index]!r}, expected {name!r}")
    if not header[7].startswith("Cross-ancestry") or sub[9] != "p-value":
        raise ValueError("Supplementary Table 4: cross-ancestry p-value is not in column 10")
    rows = []
    for number, row in enumerate(sheet[header_at + 3:], start=header_at + 4):
        if row and row[0].startswith(ST4_SECTION_END):
            break
        if not row or not row[0].startswith("rs"):
            continue
        p_value = float(row[9])
        rows.append({
            "evidence_id": f"KOLLER2026_ST4:{row[0]}",
            "source": KOLLER_SOURCE,
            "source_version": KOLLER_VERSION,
            "source_record": f"Supplementary Table 4 sheet row {number}",
            "original_id": row[5].strip() or gs.MISSING,
            "variant_id": row[0],
            "p_value": repr(p_value),
            "evidence_type": "gwas_association",
            "evidence_level": "locus",
            "population": KOLLER_POPULATION,
            "gene_id": gs.MISSING,
            "gene_symbol": gs.MISSING,
            "mapping_status": "pending" if p_value < THRESHOLD else "excluded",
            "mapping_method": "not_run",
            "notes": (f"chr{row[1]}:{row[2]}; {row[6].strip()}; nearest gene as reported by the authors"
                      + ("" if p_value < THRESHOLD else f"; {EXCLUDE_THRESHOLD}")),
        })
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--catalog-dir", type=Path, required=True,
                        help="cache folder for GWAS Catalog API responses")
    parser.add_argument("--fetch", action="store_true", help="download missing Catalog responses")
    parser.add_argument("--koller-xlsx", type=Path, required=True,
                        help="Koller et al. 2026 supplementary tables (MOESM4)")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        associations, studies, stats = load_catalog(args.catalog_dir, args.fetch)
        catalog, reasons = catalog_rows(associations, studies, stats, args.catalog_dir, args.fetch)
        koller = koller_rows(read_xlsx_sheet(args.koller_xlsx, ST4_SHEET))
        gs.write_tsv(args.output, gs.RISK_EVIDENCE_COLUMNS, catalog + koller)
    except (OSError, ValueError, KeyError, zipfile.BadZipFile, ET.ParseError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    print(json.dumps({
        "catalog_release": stats["date"],
        "catalog_associations": len(catalog),
        "catalog_by_outcome": dict(reasons),
        "catalog_included_studies": len({r["source_record"].split()[0] for r in catalog
                                         if r["mapping_status"] == "pending"}),
        "koller_st4_rows": len(koller),
        "koller_st4_included": sum(r["mapping_status"] == "pending" for r in koller),
        "evidence_rows": len(catalog) + len(koller),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
