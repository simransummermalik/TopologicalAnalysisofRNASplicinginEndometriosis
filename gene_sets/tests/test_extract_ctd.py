"""Tests for extract_ctd.py using a small, made-up CTD report.

Nothing here is real evidence. Chemicals, genes, and PubMed IDs are fixtures.
"""

from __future__ import annotations

import gzip
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

GENE_SETS_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = GENE_SETS_DIR.parent
sys.path.insert(0, str(GENE_SETS_DIR))

import extract_ctd as ex  # noqa: E402
import gene_sets as gs  # noqa: E402

HUMAN = "Homo sapiens\t9606"
RAT = "Rattus norvegicus\t10116"


def row(chem, cid, sym, gid, forms, organism, text, actions, pmids):
    return "\t".join([chem, cid, "", sym, gid, forms, organism, text, actions, pmids])


CTD = "\n".join(
    [
        "# The Comparative Toxicogenomics Database (CTD) - fixture",
        "# Report created: Thu Oct 29 12:43:06 EDT 2026",
        "#",
        "# Fields:",
        "# " + "\t".join(ex.CTD_FIELDS),
        "#",
        # line 7: kept, direct mRNA
        row("TCDD", "D013749", "WNT4", "1", "mRNA", HUMAN,
            "TCDD results in increased expression of WNT4 mRNA", "increases^expression", "11|12"),
        # line 8: chemical not on the list
        row("Bisphenol A", "C006780", "WNT4", "1", "mRNA", HUMAN,
            "Bisphenol A results in increased expression of WNT4 mRNA", "increases^expression", "13"),
        # line 9: wrong species
        row("TCDD", "D013749", "WNT4", "1", "mRNA", RAT,
            "TCDD results in increased expression of WNT4 mRNA", "increases^expression", "14"),
        # line 10: binding only
        row("TCDD", "D013749", "KDR", "4", "protein", HUMAN,
            "TCDD binds to KDR protein", "affects^binding", "15"),
        # line 11: the row's gene only takes part; GREB1 is the measured gene
        row("TCDD", "D013749", "AHR", "9", "protein", HUMAN,
            "AHR protein affects the reaction [TCDD results in increased expression of GREB1 mRNA]",
            "affects^reaction|increases^expression", "16"),
        # line 12: kept, modulated reaction, protein
        row("PCB 126", "C023035", "GREB1", "3", "protein", HUMAN,
            "Fulvestrant inhibits the reaction [PCB 126 results in increased expression of GREB1 protein]",
            "decreases^reaction|increases^expression", "17"),
        # line 13: kept, co-treatment, activity
        row("TCDD", "D013749", "KDR", "4", "protein", HUMAN,
            "[TCDD co-treated with Estradiol] results in decreased activity of KDR protein",
            "affects^cotreatment|decreases^activity", "18"),
        # line 14: kept, form not stated
        row("TCDD", "D013749", "IL15", "8", "", HUMAN,
            "TCDD results in increased expression of IL15", "increases^expression", "19"),
        # line 15: kept, mRNA and protein in one interaction -> two evidence rows
        row("TCDD", "D013749", "BCL2L11", "2", "mRNA|protein", HUMAN,
            "TCDD results in increased expression of BCL2L11 mRNA and results in increased expression of BCL2L11 protein",
            "increases^expression", "20"),
        # line 16: kept, CTD's joined phrasing
        row("TCDD", "D013749", "KDR", "4", "protein", HUMAN,
            "TCDD affects the expression of and affects the secretion of KDR protein",
            "affects^expression|affects^secretion", "21"),
        # line 17: secretion only, not a measured response under the rule
        row("TCDD", "D013749", "WNT4", "1", "protein", HUMAN,
            "TCDD results in increased secretion of WNT4 protein", "increases^secretion", "22"),
    ]
) + "\n"

CHEMICALS = {"D013749": "2,3,7,8-TCDD", "C023035": "PCB 126"}


class ExtractTests(unittest.TestCase):
    def setUp(self):
        self.rows, self.counts = ex.extract(io.StringIO(CTD), CHEMICALS, "fixture.tsv")
        self.by_id = {r["evidence_id"]: r for r in self.rows}

    def test_filter_counts_in_order(self):
        self.assertEqual(
            self.counts["rows_after_each_filter"],
            {"input": 11, "chemical": 10, "species": 9, "measured": 7, "target": 6},
        )
        self.assertEqual(self.counts["per_chemical"]["MESH:C023035 (PCB 126)"]["target"], 1)
        self.assertEqual(self.counts["source_version"], "Thu Oct 29 12:43:06 EDT 2026")

    def test_direct_row_fields(self):
        r = self.by_id["CTD:7:expression"]
        self.assertEqual(r["original_id"], "WNT4")
        self.assertEqual(r["exposure_id"], "MESH:D013749")
        self.assertEqual(r["publications"], "PMID:11|PMID:12")
        self.assertEqual(r["context"], "direct")
        self.assertEqual(r["mapping_status"], "pending")
        self.assertEqual(r["source_record"], "fixture.tsv line 7")
        self.assertEqual((r["tissue"], r["dose_duration"], r["assay"]), ("NA", "NA", "NA"))
        self.assertIn("interaction: TCDD results", r["notes"])

    def test_participant_gene_is_dropped(self):
        self.assertFalse(any(r["original_id"] == "AHR" for r in self.rows))

    def test_contexts_and_response_types(self):
        self.assertEqual(self.by_id["CTD:12:protein"]["context"], "modulated_reaction")
        self.assertEqual(self.by_id["CTD:13:activity"]["context"], "cotreatment")
        self.assertEqual(self.by_id["CTD:16:protein"]["original_id"], "KDR")
        self.assertIn("does not state mRNA or protein", self.by_id["CTD:14:expression"]["notes"])
        self.assertNotIn("does not state", self.by_id["CTD:7:expression"]["notes"])

    def test_one_interaction_with_two_forms_gives_two_rows(self):
        self.assertIn("CTD:15:expression", self.by_id)
        self.assertIn("CTD:15:protein", self.by_id)

    def test_output_passes_check(self):
        self.assertEqual(
            gs.check_table("environment_evidence", gs.ENVIRONMENT_EVIDENCE_COLUMNS, self.rows), []
        )

    def test_missing_report_line_is_an_error(self):
        body = "\n".join(line for line in CTD.splitlines() if "Report created" not in line)
        with self.assertRaises(ValueError):
            ex.extract(io.StringIO(body), CHEMICALS, "fixture.tsv")

    def test_wrong_field_count_is_an_error(self):
        with self.assertRaises(ValueError):
            ex.extract(io.StringIO(CTD + "too\tfew\n"), CHEMICALS, "fixture.tsv")


class ChemicalListTests(unittest.TestCase):
    def test_project_list_has_21_unique_ids(self):
        chemicals = ex.read_chemicals(REPO_ROOT / "data/gene_sets/dioxin_like_chemicals.tsv")
        self.assertEqual(len(chemicals), 21)
        self.assertEqual(chemicals["D013749"], "2,3,7,8-TCDD")


class CommandLineTests(unittest.TestCase):
    def test_main_writes_table_from_gzip(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            ctd = tmp / "ctd.tsv.gz"
            with gzip.open(ctd, "wt", encoding="utf-8") as handle:
                handle.write(CTD)
            chems = tmp / "chems.tsv"
            chems.write_text(
                "congener\tclass\tctd_name\tctd_chemical_id\n"
                "2,3,7,8-TCDD\tPCDD\tTCDD\tMESH:D013749\n",
                encoding="utf-8",
            )
            out = tmp / "env.tsv"
            stdout = io.StringIO()
            with redirect_stdout(stdout):
                code = ex.main(["--ctd", str(ctd), "--chemicals", str(chems), "--output", str(out)])
            self.assertEqual(code, 0)
            self.assertEqual(json.loads(stdout.getvalue())["rows_after_each_filter"]["target"], 5)
            header, rows = gs.read_tsv(out)
            self.assertEqual(header, gs.ENVIRONMENT_EVIDENCE_COLUMNS)
            self.assertEqual(len(rows), 6)


if __name__ == "__main__":
    unittest.main()
