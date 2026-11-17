"""Tests for review_ctd_context.py using made-up rows and PubMed records."""

from __future__ import annotations

import io
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

GENE_SETS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(GENE_SETS_DIR))

import gene_sets as gs  # noqa: E402
import review_ctd_context as rv  # noqa: E402

PUBMED_XML = """<?xml version="1.0"?>
<PubmedArticleSet>
 <PubmedArticle><MedlineCitation><PMID>101</PMID>
  <Article><ArticleTitle>TCDD and <i>CYP1A1</i> in liver cells</ArticleTitle></Article>
  <MeshHeadingList><MeshHeading><DescriptorName>Hep G2 Cells</DescriptorName></MeshHeading></MeshHeadingList>
 </MedlineCitation></PubmedArticle>
 <PubmedArticle><MedlineCitation><PMID>102</PMID>
  <Article><ArticleTitle>Dioxin and glycodelin</ArticleTitle></Article>
  <MeshHeadingList><MeshHeading><DescriptorName>Endometrium</DescriptorName></MeshHeading></MeshHeadingList>
 </MedlineCitation></PubmedArticle>
 <PubmedArticle><MedlineCitation><PMID>103</PMID>
  <Article><ArticleTitle>TCDD in human granulosa cells</ArticleTitle></Article>
 </MedlineCitation></PubmedArticle>
</PubmedArticleSet>
"""


def evidence_row(evidence_id, pmids, context, status="mapped"):
    mapped = status == "mapped"
    return {
        "evidence_id": evidence_id, "source": "CTD", "source_version": "fixture",
        "source_record": "fixture line 1", "publications": pmids, "original_id": "WNT4",
        "exposure": "TCDD", "exposure_id": "MESH:D013749", "dose_duration": "NA",
        "species": "Homo sapiens", "tissue": "NA", "assay": "NA",
        "response_type": "expression", "evidence_type": "measured_response",
        "context": context, "gene_id": "ENSG00000000001.5" if mapped else "NA",
        "gene_symbol": "WNT4" if mapped else "NA", "mapping_status": status,
        "mapping_method": "gencode_symbol", "notes": "CTD GeneID 1",
    }


class ClassifyTests(unittest.TestCase):
    def test_parse_and_classify(self):
        records = rv.parse_pubmed_xml(PUBMED_XML)
        self.assertEqual(records["101"][1], "TCDD and CYP1A1 in liver cells")
        self.assertEqual(rv.classify_publication(*records["101"]), ("other_or_unstated", "NA"))
        self.assertEqual(rv.classify_publication(*records["102"]),
                         ("endometrium_uterus", "MeSH:Endometrium"))
        self.assertEqual(rv.classify_publication(*records["103"]),
                         ("other_female_reproductive", "title:granulosa"))

    def test_broad_terms_do_not_count(self):
        self.assertEqual(rv.classify_publication(["Pregnancy", "Oocytes"], "hERG in oocytes"),
                         ("other_or_unstated", "NA"))


class ApplyTests(unittest.TestCase):
    TIERS = {"101": "other_or_unstated", "102": "endometrium_uterus"}

    def test_tissue_takes_most_relevant_paper(self):
        rows, _ = rv.apply_review([evidence_row("a", "PMID:101|PMID:102", "direct")], self.TIERS)
        self.assertEqual(rows[0]["tissue"], "cited_paper:endometrium_uterus")
        self.assertEqual(rows[0]["mapping_status"], "mapped")

    def test_cotreatment_excluded_with_reason_and_rerun_is_stable(self):
        rows = [
            evidence_row("a", "PMID:101", "cotreatment"),
            evidence_row("b", "PMID:101", "cotreatment|dioxin_effect_modified"),
            evidence_row("c", "PMID:101", "dioxin_modifies_other"),
        ]
        once, counts = rv.apply_review(rows, self.TIERS)
        self.assertEqual([r["mapping_status"] for r in once], ["excluded", "excluded", "mapped"])
        self.assertIn("was mapped to ENSG00000000001.5 WNT4", once[0]["notes"])
        self.assertEqual(once[0]["gene_id"], "NA")
        self.assertEqual(counts["newly_excluded_cotreatment"], 2)
        twice, counts = rv.apply_review(once, self.TIERS)
        self.assertEqual(twice, once)
        self.assertEqual(counts["newly_excluded_cotreatment"], 0)
        self.assertEqual(gs.check_table("environment_evidence", gs.ENVIRONMENT_EVIDENCE_COLUMNS, twice), [])

    def test_untagged_pmid_is_an_error(self):
        with self.assertRaises(ValueError):
            rv.apply_review([evidence_row("a", "PMID:999", "direct")], self.TIERS)


class CommandLineTests(unittest.TestCase):
    def test_tag_then_apply(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            xml_dir = tmp / "pubmed"
            xml_dir.mkdir()
            (xml_dir / "pubmed_000.xml").write_text(PUBMED_XML, encoding="utf-8")
            evidence = tmp / "env.tsv"
            gs.write_tsv(evidence, gs.ENVIRONMENT_EVIDENCE_COLUMNS,
                         [evidence_row("a", "PMID:102", "direct"),
                          evidence_row("b", "PMID:103", "cotreatment")])
            pubs = tmp / "pubs.tsv"
            with redirect_stdout(io.StringIO()):
                self.assertEqual(rv.main(["tag-publications", "--evidence", str(evidence),
                                          "--pubmed-dir", str(xml_dir), "--output", str(pubs)]), 0)
                self.assertEqual(rv.main(["apply", "--evidence", str(evidence),
                                          "--publications", str(pubs), "--output", str(evidence)]), 0)
            _, rows = gs.read_tsv(evidence)
            self.assertEqual([r["tissue"] for r in rows],
                             ["cited_paper:endometrium_uterus", "cited_paper:other_female_reproductive"])
            self.assertEqual(rows[1]["mapping_status"], "excluded")


if __name__ == "__main__":
    unittest.main()
