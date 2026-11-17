"""Tests for audit_environment.py: a clean fixture pipeline passes, tampering fails."""

from __future__ import annotations

import gzip
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

import audit_environment as audit  # noqa: E402
import extract_ctd as ex  # noqa: E402
import gene_sets as gs  # noqa: E402
import review_ctd_context as rv  # noqa: E402
from test_extract_ctd import CTD  # noqa: E402
from test_gene_sets import GTF, HGNC  # noqa: E402

CHEMICALS = (
    "congener\tclass\tctd_name\tctd_chemical_id\n"
    "2,3,7,8-TCDD\tPCDD\tTCDD\tMESH:D013749\n"
    "PCB 126\tnon-ortho PCB\tPCB 126\tMESH:C023035\n"
)


class AuditTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        tmp = Path(self._tmp.name)
        self.paths = {name: tmp / name for name in (
            "CTD_chem_gene_ixns.tsv.gz", "chemicals.tsv", "annotation.gtf", "hgnc.tsv",
            "evidence.tsv", "membership.tsv", "sample.tsv")}
        with gzip.open(self.paths["CTD_chem_gene_ixns.tsv.gz"], "wt", encoding="utf-8") as handle:
            handle.write(CTD)
        self.paths["chemicals.tsv"].write_text(CHEMICALS)
        self.paths["annotation.gtf"].write_text(GTF)
        self.paths["hgnc.tsv"].write_text(HGNC)

        chemicals = ex.read_chemicals(self.paths["chemicals.tsv"])
        rows, _ = ex.extract(io.StringIO(CTD), chemicals, "CTD_chem_gene_ixns.tsv.gz")
        reference = gs.load_reference(self.paths["annotation.gtf"], self.paths["hgnc.tsv"])
        rows, _, _ = gs.map_evidence(rows, reference)
        pmids = {p.removeprefix("PMID:") for r in rows for p in r["publications"].split("|")}
        rows, _ = rv.apply_review(rows, {p: "other_or_unstated" for p in pmids})
        self.rows = rows
        gs.write_tsv(self.paths["evidence.tsv"], gs.ENVIRONMENT_EVIDENCE_COLUMNS, rows)
        membership, _ = gs.dedupe_environment(rows, 1)
        gs.write_tsv(self.paths["membership.tsv"], gs.ENVIRONMENT_MEMBERSHIP_COLUMNS, membership)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def run_audit(self) -> tuple[int, dict]:
        p = self.paths
        out = io.StringIO()
        with redirect_stdout(out):
            code = audit.main([
                "--evidence", str(p["evidence.tsv"]), "--membership", str(p["membership.tsv"]),
                "--chemicals", str(p["chemicals.tsv"]), "--ctd", str(p["CTD_chem_gene_ixns.tsv.gz"]),
                "--gtf", str(p["annotation.gtf"]), "--hgnc", str(p["hgnc.tsv"]),
                "--min-publications", "1", "--sample-size", "3", "--sample-output", str(p["sample.tsv"]),
            ])
        return code, json.loads(out.getvalue())

    def test_clean_pipeline_passes(self) -> None:
        code, report = self.run_audit()
        self.assertEqual(code, 0, report)
        self.assertEqual(report["C_rederive"]["rederived_lines"], 7)
        self.assertEqual(len(gs.read_tsv(self.paths["sample.tsv"])[1]), 3)

    def test_tampered_exposure_fails_provenance(self) -> None:
        self.rows[0]["exposure"] = "Bisphenol A"
        gs.write_tsv(self.paths["evidence.tsv"], gs.ENVIRONMENT_EVIDENCE_COLUMNS, self.rows)
        code, report = self.run_audit()
        self.assertEqual(code, 1)
        self.assertFalse(report["B_provenance"]["pass"])

    def test_dropped_membership_row_fails_dedupe(self) -> None:
        header, membership = gs.read_tsv(self.paths["membership.tsv"])
        gs.write_tsv(self.paths["membership.tsv"], header, membership[1:])
        code, report = self.run_audit()
        self.assertEqual(code, 1)
        self.assertFalse(report["F_dedupe"]["pass"])


if __name__ == "__main__":
    unittest.main()
