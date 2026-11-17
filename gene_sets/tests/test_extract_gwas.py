"""Tests for extract_gwas.py with a made-up Catalog cache and a tiny .xlsx.

Nothing here is real evidence. Accessions, SNPs, and genes are fixtures.
"""

from __future__ import annotations

import io
import json
import sys
import tempfile
import unittest
import zipfile
from contextlib import redirect_stdout
from pathlib import Path

GENE_SETS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(GENE_SETS_DIR))

import extract_gwas as ex  # noqa: E402
import gene_sets as gs  # noqa: E402


def assoc(assoc_id, p, accession="GCST1", pmid="111", genes=("WNT4",), efo=("endometriosis",),
          trait="Endometriosis"):
    return {
        "association_id": assoc_id, "p_value": p, "accession_id": accession, "pubmed_id": pmid,
        "first_author": "Fixture A", "mapped_genes": list(genes), "locations": ["1:100"],
        "snp_allele": [{"rs_id": f"rs{assoc_id}", "effect_allele": "A"}],
        "efo_traits": [{"efo_trait": e} for e in efo], "reported_trait": [trait],
    }


ASSOCIATIONS = [
    assoc(1, 1e-10),                                              # included
    assoc(2, 1e-6),                                               # below threshold
    assoc(3, 1e-9, efo=("endometriosis", "asthma")),              # cross-trait
    assoc(4, 1e-9, trait="Endometriosis (MTAG)"),                 # MTAG
    assoc(5, 1e-9, accession="GCST_K_EUR", pmid=ex.KOLLER_PMID),  # other Koller analysis
    assoc(6, 1e-9, accession=ex.KOLLER_HEADLINE, pmid=ex.KOLLER_PMID, genes=("GREB1", "KDR")),
    assoc(7, 1e-9, genes=()),                                     # no mapped gene
]
STUDIES = {acc: {"disease_trait": "Endometriosis", "initial_sample_size": "10 cases, 20 controls"}
           for acc in ("GCST1", "GCST_K_EUR", ex.KOLLER_HEADLINE)}
STATS = {"date": "2026-10-29", "ensemblbuild": "116", "dbsnpbuild": "156"}


def write_cache(folder: Path) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "stats.json").write_text(json.dumps(STATS))
    (folder / "trait_p0.json").write_text(json.dumps(
        {"_embedded": {"associations": ASSOCIATIONS}, "page": {"totalPages": 1}}))
    for accession, study in STUDIES.items():
        (folder / f"study_{accession}.json").write_text(json.dumps(study))
    for a in ASSOCIATIONS:
        reported = [{"gene_name": "WNT4"}] if a["association_id"] == 1 else []
        (folder / f"loci_{a['association_id']}.json").write_text(json.dumps(
            {"_embedded": {"loci": [{"author_reported_genes": reported}]}}))


def inline_row(number: int, values: list[str]) -> str:
    cells = []
    for index, value in enumerate(values):
        if value == "":
            continue
        column = chr(65 + index) if index < 26 else "A" + chr(65 + index - 26)
        cells.append(f'<c r="{column}{number}" t="inlineStr"><is><t>{value}</t></is></c>')
    return f'<row r="{number}">{"".join(cells)}</row>'


def write_xlsx(path: Path) -> None:
    header = ["Lead SNP", "Chromosome", "Position", "Effect alelle", "Other allele",
              "Nearest gene", "Type", "Cross-ancestry (AFR, AMR, EAS, EUR, MID, SAS)"]
    sub = [""] * 7 + ["Beta", "SE", "p-value"]
    rows = [
        ["Supplemental Table 4. fixture"], [], header, [""] * 7 + ["Endometriosis combined"], sub,
        ["rs10", "1", "100", "A", "G", "WNT4", "Novel ", "0.1", "0.01", "4.5E-10"],
        ["rs11", "2", "200", "A", "G", "RP11-1A1.1:KDR", "Previously reported", "0.1", "0.01", "1E-9"],
        [ex.ST4_SECTION_END],
        ["rs12", "3", "300", "A", "G", "GREB1", "Clinical endometriosis, novel", "0.1", "0.01", "1E-5"],
    ]
    sheet = ('<?xml version="1.0"?><worksheet xmlns="http://schemas.openxmlformats.org/'
             'spreadsheetml/2006/main"><sheetData>'
             + "".join(inline_row(i + 1, r) for i, r in enumerate(rows) if r)
             + "</sheetData></worksheet>")
    workbook = ('<?xml version="1.0"?><workbook xmlns="http://schemas.openxmlformats.org/'
                'spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/'
                'officeDocument/2006/relationships"><sheets>'
                f'<sheet name="{ex.ST4_SHEET}" sheetId="1" r:id="rId1"/></sheets></workbook>')
    rels = ('<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/'
            'package/2006/relationships"><Relationship Id="rId1" Type="worksheet" '
            'Target="worksheets/sheet1.xml"/></Relationships>')
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("xl/workbook.xml", workbook)
        z.writestr("xl/_rels/workbook.xml.rels", rels)
        z.writestr("xl/worksheets/sheet1.xml", sheet)


class CatalogTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.cache = Path(self._tmp.name) / "cache"
        write_cache(self.cache)

    def tearDown(self):
        self._tmp.cleanup()

    def rows(self):
        associations, studies, stats = ex.load_catalog(self.cache, fetch=False)
        rows, reasons = ex.catalog_rows(associations, studies, stats, self.cache, fetch=False)
        return {r["evidence_id"]: r for r in rows}, reasons

    def test_exclusion_rules_in_order(self):
        rows, reasons = self.rows()
        status = {k: v["mapping_status"] for k, v in rows.items()}
        self.assertEqual(status, {
            "GWASCAT:1": "pending", "GWASCAT:2": "excluded", "GWASCAT:3": "excluded",
            "GWASCAT:4": "excluded", "GWASCAT:5": "excluded", "GWASCAT:6": "pending",
            "GWASCAT:7": "pending"})
        self.assertIn("p >= 5e-8", rows["GWASCAT:2"]["notes"])
        self.assertIn("pleiotropy) analysis with asthma", rows["GWASCAT:3"]["notes"])
        self.assertIn("MTAG", rows["GWASCAT:4"]["notes"])
        self.assertIn("Koller", rows["GWASCAT:5"]["notes"])
        self.assertEqual(reasons["included"], 3)

    def test_fields(self):
        rows, _ = self.rows()
        self.assertEqual(rows["GWASCAT:6"]["original_id"], "GREB1, KDR")
        self.assertEqual(rows["GWASCAT:7"]["original_id"], "NA")
        self.assertEqual(rows["GWASCAT:1"]["variant_id"], "rs1")
        self.assertIn("author-reported genes: WNT4", rows["GWASCAT:1"]["notes"])
        self.assertIn("release 2026-10-29", rows["GWASCAT:1"]["source_version"])
        self.assertEqual(gs.check_table("risk_evidence", gs.RISK_EVIDENCE_COLUMNS, list(rows.values())), [])

    def test_missing_cache_without_fetch_is_an_error(self):
        (self.cache / "stats.json").unlink()
        with self.assertRaises(ValueError):
            ex.load_catalog(self.cache, fetch=False)


class KollerTests(unittest.TestCase):
    def test_reads_first_section_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "supp.xlsx"
            write_xlsx(path)
            rows = ex.koller_rows(ex.read_xlsx_sheet(path, ex.ST4_SHEET))
        self.assertEqual([r["variant_id"] for r in rows], ["rs10", "rs11"])
        self.assertEqual(rows[1]["original_id"], "RP11-1A1.1:KDR")
        self.assertEqual(rows[0]["p_value"], repr(4.5e-10))
        self.assertEqual(rows[0]["evidence_level"], "locus")
        self.assertIn("Novel", rows[0]["notes"])
        self.assertEqual(gs.check_table("risk_evidence", gs.RISK_EVIDENCE_COLUMNS, rows), [])

    def test_wrong_layout_is_refused(self):
        sheet = [["Lead SNP", "Chromosome", "Position", "", "", "Gene", "Type", "Cross-ancestry"],
                 [], ["", "", "", "", "", "", "", "Beta", "SE", "p-value"]]
        with self.assertRaises(ValueError):
            ex.koller_rows(sheet)


class CommandLineTests(unittest.TestCase):
    def test_main(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            write_cache(tmp / "cache")
            write_xlsx(tmp / "supp.xlsx")
            out = io.StringIO()
            with redirect_stdout(out):
                code = ex.main(["--catalog-dir", str(tmp / "cache"), "--koller-xlsx", str(tmp / "supp.xlsx"),
                                "--output", str(tmp / "risk.tsv")])
            self.assertEqual(code, 0)
            counts = json.loads(out.getvalue())
            self.assertEqual((counts["catalog_associations"], counts["koller_st4_rows"]), (7, 2))
            header, rows = gs.read_tsv(tmp / "risk.tsv")
            self.assertEqual(header, gs.RISK_EVIDENCE_COLUMNS)
            self.assertEqual(len(rows), 9)


if __name__ == "__main__":
    unittest.main()
