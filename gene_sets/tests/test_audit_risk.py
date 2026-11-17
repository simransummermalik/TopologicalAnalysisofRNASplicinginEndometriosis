"""Tests for audit_risk.py: a clean fixture pipeline passes, tampering fails."""

from __future__ import annotations

import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TESTS_DIR.parent))
sys.path.insert(0, str(TESTS_DIR))

import audit_risk  # noqa: E402
import extract_gwas as ex  # noqa: E402
import gene_sets as gs  # noqa: E402
from test_extract_gwas import write_cache, write_xlsx  # noqa: E402
from test_gene_sets import GTF, HGNC  # noqa: E402

REVIEW = (
    "evidence_id\tdecision\tgene_id\treason\n"
    "GWASCAT:6\texcluded\tNA\tCatalog copy of a Koller locus\n"
    "GWASCAT:7\texcluded\tNA\tno gene\n"
    "KOLLER2026_ST4:rs11\tmapped\tENSG00000000004.2\tclone name translated\n"
)


class AuditRiskTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        tmp = Path(self._tmp.name)
        self.p = {name: tmp / name for name in (
            "cache", "supp.xlsx", "annotation.gtf", "hgnc.tsv", "review.tsv",
            "evidence.tsv", "membership.tsv", "sample.tsv")}
        write_cache(self.p["cache"])
        write_xlsx(self.p["supp.xlsx"])
        self.p["annotation.gtf"].write_text(GTF)
        self.p["hgnc.tsv"].write_text(HGNC)
        self.p["review.tsv"].write_text(REVIEW)

        associations, studies, stats = ex.load_catalog(self.p["cache"], fetch=False)
        rows, _ = ex.catalog_rows(associations, studies, stats, self.p["cache"], fetch=False)
        rows += ex.koller_rows(ex.read_xlsx_sheet(self.p["supp.xlsx"], ex.ST4_SHEET))
        reference = gs.load_reference(self.p["annotation.gtf"], self.p["hgnc.tsv"])
        rows, _, _ = gs.map_evidence(rows, reference)
        key, decisions = gs.read_mapping_review(self.p["review.tsv"], reference)
        rows, _ = gs.apply_mapping_review(rows, decisions, reference, key)
        self.rows = rows
        gs.write_tsv(self.p["evidence.tsv"], gs.RISK_EVIDENCE_COLUMNS, rows)
        membership, _ = gs.dedupe_risk(rows)
        gs.write_tsv(self.p["membership.tsv"], gs.RISK_MEMBERSHIP_COLUMNS, membership)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def run_audit(self) -> tuple[int, dict]:
        p = self.p
        out = io.StringIO()
        with redirect_stdout(out):
            code = audit_risk.main([
                "--evidence", str(p["evidence.tsv"]), "--membership", str(p["membership.tsv"]),
                "--review", str(p["review.tsv"]), "--catalog-dir", str(p["cache"]),
                "--koller-xlsx", str(p["supp.xlsx"]), "--gtf", str(p["annotation.gtf"]),
                "--gtf-v19", str(p["annotation.gtf"]), "--sample-size", "2",
                "--sample-output", str(p["sample.tsv"]),
            ])
        return code, json.loads(out.getvalue())

    def test_clean_pipeline_passes(self) -> None:
        code, report = self.run_audit()
        self.assertEqual(code, 0, report)
        self.assertEqual(report["A_schema"]["membership_rows"], 2)  # WNT4 and KDR
        self.assertEqual(report["H_proximity"]["no_coordinates"], 0)

    def test_tampered_p_value_fails_provenance(self) -> None:
        row = next(r for r in self.rows if r["evidence_id"] == "GWASCAT:1")
        row["p_value"] = "1e-12"
        gs.write_tsv(self.p["evidence.tsv"], gs.RISK_EVIDENCE_COLUMNS, self.rows)
        code, report = self.run_audit()
        self.assertEqual(code, 1)
        self.assertFalse(report["B_provenance"]["pass"])

    def test_wrong_rule_outcome_fails_rederive(self) -> None:
        row = next(r for r in self.rows if r["evidence_id"] == "GWASCAT:4")  # MTAG row
        row.update(mapping_status="mapped", gene_id="ENSG00000000001.5", gene_symbol="WNT4",
                   mapping_method="gencode_symbol")
        gs.write_tsv(self.p["evidence.tsv"], gs.RISK_EVIDENCE_COLUMNS, self.rows)
        code, report = self.run_audit()
        self.assertEqual(code, 1)
        self.assertFalse(report["C_rederive"]["pass"])


if __name__ == "__main__":
    unittest.main()
