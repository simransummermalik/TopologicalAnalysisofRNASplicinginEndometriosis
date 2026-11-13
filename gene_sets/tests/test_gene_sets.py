"""Tests for gene_sets.py using small, made-up genes and evidence.

Nothing here is real evidence. Gene IDs, symbols, and sources are fixtures.
"""

from __future__ import annotations

import io
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

GENE_SETS_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = GENE_SETS_DIR.parent
sys.path.insert(0, str(GENE_SETS_DIR))

import gene_sets as gs  # noqa: E402

GTF = """\
##description: tiny fixture
chr1\tFIXTURE\tgene\t100\t900\t.\t+\t.\tgene_id "ENSG00000000001.5"; gene_type "protein_coding"; gene_name "WNT4";
chr1\tFIXTURE\texon\t100\t200\t.\t+\t.\tgene_id "ENSG00000000001.5"; gene_name "WNT4";
chr2\tFIXTURE\tgene\t100\t900\t.\t-\t.\tgene_id "ENSG00000000002.3"; gene_type "protein_coding"; gene_name "BCL2L11";
chr2\tFIXTURE\tgene\t1000\t1900\t.\t+\t.\tgene_id "ENSG00000000003.1"; gene_type "protein_coding"; gene_name "GREB1";
chr2\tFIXTURE\tgene\t2000\t2900\t.\t+\t.\tgene_id "ENSG00000000004.2"; gene_type "protein_coding"; gene_name "KDR";
chrX\tFIXTURE\tgene\t100\t900\t.\t+\t.\tgene_id "ENSG00000000005.8"; gene_type "protein_coding"; gene_name "PARGENE";
chrY\tFIXTURE\tgene\t100\t900\t.\t+\t.\tgene_id "ENSG00000000005.8_PAR_Y"; gene_type "protein_coding"; gene_name "PARGENE";
chr3\tFIXTURE\tgene\t100\t900\t.\t+\t.\tgene_id "ENSG00000000006.1"; gene_type "lncRNA"; gene_name "DUPNAME";
chr4\tFIXTURE\tgene\t100\t900\t.\t+\t.\tgene_id "ENSG00000000007.1"; gene_type "lncRNA"; gene_name "DUPNAME";
chr5\tFIXTURE\tgene\t100\t900\t.\t+\t.\tgene_id "ENSG00000000008.4"; gene_type "protein_coding"; gene_name "IL15";
chr6\tFIXTURE\tgene\t100\t900\t.\t+\t.\tgene_id "ENSG00000000009.1"; gene_type "protein_coding"; gene_name "SHAREDA";
chr7\tFIXTURE\tgene\t100\t900\t.\t+\t.\tgene_id "ENSG00000000010.1"; gene_type "protein_coding"; gene_name "SHAREDB";
chr8\tFIXTURE\tgene\t100\t900\t.\t+\t.\tgene_id "ENSG00000000011.9"; gene_type "protein_coding"; gene_name "LINKED";
chr8\tFIXTURE\tgene\t1000\t1900\t.\t+\t.\tgene_id "ENSG00000000012.1"; gene_type "protein_coding"; gene_name "LINKED";
""" + "".join(
    f'chr9\tFIXTURE\tgene\t{i}\t{i + 50}\t.\t+\t.\tgene_id "ENSG{900000 + i:011d}.1"; '
    f'gene_type "misc_RNA"; gene_name "Y_RNA";\n'
    for i in range(25)
)

HGNC = (
    "hgnc_id\tsymbol\talias_symbol\tprev_symbol\tensembl_gene_id\n"
    "HGNC:1\tBCL2L11\tBIM\t\tENSG00000000002\n"
    "HGNC:2\tIL15\t\"IL-15\"\t\tENSG00000000008\n"
    "HGNC:3\tSHAREDA\tCOMMONALIAS\t\tENSG00000000009\n"
    "HGNC:4\tSHAREDB\tCOMMONALIAS\tOLDB\tENSG00000000010\n"
    "HGNC:5\tWITHDRAWNLIKE\t\t\tENSG00000099999\n"
    "HGNC:6\tLINKED\t\t\tENSG00000000011\n"
)


def risk_row(evidence_id: str, original_id: str, **overrides: str) -> dict[str, str]:
    row = {
        "evidence_id": evidence_id,
        "source": "FIXTURE_GWAS",
        "source_version": "fixture-1",
        "source_record": f"table1:{evidence_id}",
        "original_id": original_id,
        "variant_id": "rs0000001",
        "p_value": "1e-9",
        "evidence_type": "gwas_association",
        "evidence_level": "locus",
        "population": "fixture",
        "gene_id": "NA",
        "gene_symbol": "NA",
        "mapping_status": "pending",
        "mapping_method": "NA",
        "notes": "NA",
    }
    row.update(overrides)
    return row


def env_row(evidence_id: str, gene_id: str, symbol: str, pubs: str, **overrides: str) -> dict[str, str]:
    row = {
        "evidence_id": evidence_id,
        "source": "FIXTURE_CTD",
        "source_version": "fixture-1",
        "source_record": f"row:{evidence_id}",
        "publications": pubs,
        "original_id": symbol,
        "exposure": "Tetrachlorodibenzodioxin",
        "exposure_id": "MESH:D013749",
        "dose_duration": "NA",
        "species": "Homo sapiens",
        "tissue": "NA",
        "assay": "NA",
        "response_type": "expression",
        "evidence_type": "measured_response",
        "context": "NA",
        "gene_id": gene_id,
        "gene_symbol": symbol,
        "mapping_status": "mapped",
        "mapping_method": "gencode_symbol",
        "notes": "NA",
    }
    row.update(overrides)
    return row


class GeneSetTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self._temp = tempfile.TemporaryDirectory()
        self.temp = Path(self._temp.name)
        self.gtf = self.temp / "annotation.gtf"
        self.gtf.write_text(GTF)
        self.hgnc = self.temp / "hgnc.tsv"
        self.hgnc.write_text(HGNC)
        self.reference = gs.load_reference(self.gtf, self.hgnc)

    def tearDown(self) -> None:
        self._temp.cleanup()

    def run_cli(self, *argv: str) -> tuple[int, str, str]:
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = gs.main(list(argv))
        return code, out.getvalue(), err.getvalue()


class MapLabelTests(GeneSetTestCase):
    def assertMapped(self, label: str, gene_id: str, method: str) -> None:
        result = gs.map_label(label, self.reference)
        self.assertEqual((result.status, result.method, result.gene_id()), ("mapped", method, gene_id))

    def test_exact_gencode_symbol(self) -> None:
        self.assertMapped("WNT4", "ENSG00000000001.5", "gencode_symbol")

    def test_unversioned_and_versioned_ensembl_ids(self) -> None:
        self.assertMapped("ENSG00000000001", "ENSG00000000001.5", "gencode_id")
        self.assertMapped("ENSG00000000001.2", "ENSG00000000001.5", "gencode_id")

    def test_alias_maps_through_hgnc(self) -> None:
        self.assertMapped("BIM", "ENSG00000000002.3", "hgnc_alias")
        self.assertMapped("IL-15", "ENSG00000000008.4", "hgnc_alias")

    def test_previous_symbol(self) -> None:
        self.assertMapped("OLDB", "ENSG00000000010.1", "hgnc_previous")

    def test_par_y_copy_does_not_make_symbol_ambiguous(self) -> None:
        self.assertMapped("PARGENE", "ENSG00000000005.8", "gencode_symbol")

    def test_duplicate_gencode_name_is_ambiguous(self) -> None:
        result = gs.map_label("DUPNAME", self.reference)
        self.assertEqual(result.status, "ambiguous")
        self.assertEqual(sorted(result.candidates), ["ENSG00000000006.1", "ENSG00000000007.1"])

    def test_duplicate_gencode_name_resolved_by_hgnc_link(self) -> None:
        self.assertMapped("LINKED", "ENSG00000000011.9", "gencode_symbol_hgnc_link")

    def test_multi_copy_name_lists_capped_candidates(self) -> None:
        result = gs.map_label("Y_RNA", self.reference)
        self.assertEqual(result.status, "ambiguous")
        self.assertEqual(len(result.candidates), 25)
        listed = gs.list_candidates(result.candidates)
        self.assertTrue(listed.endswith("|...(+5 more)"))
        self.assertEqual(listed.count("ENSG"), 20)

    def test_alias_shared_by_two_genes_is_ambiguous(self) -> None:
        result = gs.map_label("COMMONALIAS", self.reference)
        self.assertEqual(result.status, "ambiguous")
        self.assertEqual(len(result.candidates), 2)

    def test_multi_gene_label_is_ambiguous_not_split(self) -> None:
        result = gs.map_label("GREB1, KDR", self.reference)
        self.assertEqual((result.status, result.method), ("ambiguous", "multi_gene_label"))
        self.assertEqual(result.candidates, ["ENSG00000000003.1", "ENSG00000000004.2"])
        self.assertEqual(result.gene_id(), "NA")

    def test_unknown_and_missing_labels_are_unresolved(self) -> None:
        self.assertEqual(gs.map_label("NOTAGENE", self.reference).status, "unresolved")
        self.assertEqual(gs.map_label("NA", self.reference).status, "unresolved")
        self.assertEqual(gs.map_label("WITHDRAWNLIKE", self.reference).method,
                         "hgnc_approved_not_in_release")
        self.assertEqual(gs.map_label("ENSG00000012345", self.reference).status, "unresolved")


class MapEvidenceTests(GeneSetTestCase):
    def test_only_pending_rows_change_and_mapping_table_is_written(self) -> None:
        manual = risk_row("R3", "HLA class II", gene_id="NA", gene_symbol="NA",
                          mapping_status="excluded", mapping_method="manual",
                          notes="gene family; excluded in G7")
        rows = [risk_row("R1", "WNT4"), risk_row("R2", "GREB1, KDR"), manual, risk_row("R4", "WNT4")]
        evidence = self.temp / "risk_evidence.tsv"
        gs.write_tsv(evidence, gs.RISK_EVIDENCE_COLUMNS, rows)

        code, out, _ = self.run_cli(
            "map", "--kind", "risk", "--evidence", str(evidence), "--gtf", str(self.gtf),
            "--hgnc", str(self.hgnc), "--output", str(self.temp / "mapped.tsv"),
            "--mapping-table", str(self.temp / "mapping.tsv"),
        )
        self.assertEqual(code, 0)
        self.assertIn('"mapped": 2', out)

        _, mapped = gs.read_tsv(self.temp / "mapped.tsv")
        self.assertEqual(mapped[0]["gene_id"], "ENSG00000000001.5")
        self.assertEqual(mapped[0]["gene_symbol"], "WNT4")
        self.assertEqual(mapped[1]["mapping_status"], "ambiguous")
        self.assertIn("ENSG00000000003.1|ENSG00000000004.2", mapped[1]["notes"])
        self.assertEqual(mapped[2], manual)
        self.assertEqual(mapped[3]["original_id"], "WNT4")

        _, table = gs.read_tsv(self.temp / "mapping.tsv")
        self.assertEqual([r["original_id"] for r in table], ["WNT4", "GREB1, KDR"])
        self.assertEqual(table[1]["candidate_symbols"], "GREB1|KDR")

        code, _, _ = self.run_cli("check", "--kind", "risk_evidence",
                                  "--table", str(self.temp / "mapped.tsv"), "--gtf", str(self.gtf))
        self.assertEqual(code, 0)

    def test_wrong_header_is_refused(self) -> None:
        evidence = self.temp / "bad.tsv"
        evidence.write_text("gene\tsource\nWNT4\tX\n")
        code, _, err = self.run_cli(
            "map", "--kind", "risk", "--evidence", str(evidence), "--gtf", str(self.gtf),
            "--output", str(self.temp / "o.tsv"), "--mapping-table", str(self.temp / "m.tsv"),
        )
        self.assertEqual(code, 1)
        self.assertIn("header does not match", err)
        self.assertFalse((self.temp / "o.tsv").exists())


class DedupeTests(GeneSetTestCase):
    def test_risk_membership_keeps_every_evidence_id(self) -> None:
        rows = [
            risk_row("R1", "WNT4", gene_id="ENSG00000000001.5", gene_symbol="WNT4", mapping_status="mapped"),
            risk_row("R2", "WNT4", gene_id="ENSG00000000001.5", gene_symbol="WNT4",
                     mapping_status="mapped", source="FIXTURE_META"),
            risk_row("R3", "GREB1, KDR", mapping_status="ambiguous"),
        ]
        membership, stats = gs.dedupe_risk(rows)
        self.assertEqual(stats, {"evidence_rows": 3, "rows_not_mapped": 1, "genes": 1})
        self.assertEqual(membership[0]["evidence_ids"], "R1|R2")
        self.assertEqual(membership[0]["source"], "FIXTURE_GWAS|FIXTURE_META")
        self.assertEqual(gs.check_table("risk_membership", gs.RISK_MEMBERSHIP_COLUMNS, membership), [])

    def test_environment_publication_threshold_is_per_gene_chemical_pair(self) -> None:
        rows = [
            # WNT4: two rows, one publication each, different -> 2 distinct pubs.
            env_row("E1", "ENSG00000000001.5", "WNT4", "PMID:1"),
            env_row("E2", "ENSG00000000001.5", "WNT4", "PMID:2", tissue="endometrium"),
            # IL15: two rows citing the same publication -> 1 distinct pub.
            env_row("E3", "ENSG00000000008.4", "IL15", "PMID:3"),
            env_row("E4", "ENSG00000000008.4", "IL15", "PMID:3"),
            # KDR: one row, two publications.
            env_row("E5", "ENSG00000000004.2", "KDR", "PMID:4|PMID:5"),
        ]
        membership, stats = gs.dedupe_environment(rows, min_publications=2)
        self.assertEqual(stats["rows_below_publication_threshold"], 2)
        self.assertEqual(stats["genes"], 2)
        # WNT4 keeps two context rows (tissue NA and endometrium).
        wnt4 = [r for r in membership if r["gene_symbol"] == "WNT4"]
        self.assertEqual(sorted(r["tissue"] for r in wnt4), ["NA", "endometrium"])
        self.assertEqual(gs.check_table("environment_membership",
                                        gs.ENVIRONMENT_MEMBERSHIP_COLUMNS, membership), [])

        everything, _ = gs.dedupe_environment(rows, min_publications=1)
        self.assertEqual(len({r["gene_id"] for r in everything}), 3)


class CheckTests(GeneSetTestCase):
    def problems(self, kind: str, rows: list[dict[str, str]]) -> list[str]:
        return gs.check_table(kind, gs.SCHEMAS[kind], rows)

    def test_valid_rows_pass(self) -> None:
        rows = [risk_row("R1", "WNT4", gene_id="ENSG00000000001.5", gene_symbol="WNT4",
                         mapping_status="mapped", mapping_method="gencode_symbol"),
                risk_row("R2", "NOTAGENE", mapping_status="unresolved", mapping_method="no_match")]
        self.assertEqual(self.problems("risk_evidence", rows), [])
        # A locus the source reports with no gene is valid evidence to record.
        no_gene = risk_row("R3", "NA", mapping_status="unresolved", mapping_method="empty_label")
        self.assertEqual(self.problems("risk_evidence", [no_gene]), [])
        self.assertTrue(self.problems("environment_evidence", [
            env_row("E1", "NA", "NA", "PMID:1", original_id="NA", mapping_status="unresolved")]))
        self.assertEqual(self.problems("environment_evidence",
                                       [env_row("E1", "ENSG00000000001.5", "WNT4", "PMID:1")]), [])

    def test_header_only_template_passes(self) -> None:
        for name, kind in (("risk_evidence.tsv", "risk_evidence"),
                           ("environment_evidence.tsv", "environment_evidence")):
            code, _, err = self.run_cli("check", "--kind", kind,
                                        "--table", str(REPO_ROOT / "data" / "gene_sets" / name))
            self.assertEqual(code, 0, err)

    def test_detects_each_kind_of_problem(self) -> None:
        cases = {
            "blank source": risk_row("R1", "X", source=""),
            "NA in required column": risk_row("R1", "X", evidence_type="NA"),
            "unversioned gene id": risk_row("R1", "X", gene_id="ENSG00000000001",
                                            gene_symbol="WNT4", mapping_status="mapped"),
            "symbol holds an ID": risk_row("R1", "X", gene_id="ENSG00000000001.5",
                                           gene_symbol="ENSG00000000001.5", mapping_status="mapped"),
            "unmapped row with an ID": risk_row("R1", "X", gene_id="ENSG00000000001.5",
                                                mapping_status="ambiguous"),
            "bad status": risk_row("R1", "X", mapping_status="guessed"),
            "expression evidence in risk": risk_row("R1", "X", evidence_type="expression"),
            "bad p-value": risk_row("R1", "X", p_value="significant"),
        }
        for label, row in cases.items():
            with self.subTest(label):
                self.assertTrue(self.problems("risk_evidence", [row]), label)

    def test_duplicates_and_flags(self) -> None:
        member = {"gene_id": "ENSG00000000001.5", "gene_symbol": "WNT4", "evidence_type": "gwas_association",
                  "source": "S", "n_evidence": "1", "evidence_ids": "R1", "notes": "NA"}
        dup = self.problems("risk_membership", [member, dict(member)])
        self.assertTrue(any("duplicate" in p for p in dup))

        bad_flags = {"gene_id": "ENSG00000000001.5", "gene_symbol": "WNT4", "is_risk": "1",
                     "is_environment": "0", "is_overlap": "1", "category": "overlap"}
        self.assertTrue(self.problems("candidates", [bad_flags]))

    def test_gene_id_must_exist_in_annotation_when_gtf_given(self) -> None:
        table = self.temp / "risk.tsv"
        gs.write_tsv(table, gs.RISK_MEMBERSHIP_COLUMNS, [
            {"gene_id": "ENSG00000099999.1", "gene_symbol": "MISSING", "evidence_type": "gwas_association",
             "source": "S", "n_evidence": "1", "evidence_ids": "R1", "notes": "NA"}])
        code, _, _ = self.run_cli("check", "--kind", "risk_membership", "--table", str(table))
        self.assertEqual(code, 0)
        code, _, err = self.run_cli("check", "--kind", "risk_membership", "--table", str(table),
                                    "--gtf", str(self.gtf))
        self.assertEqual(code, 1)
        self.assertIn("not in the reference annotation", err)


class OverlapTests(GeneSetTestCase):
    def test_categories_and_counts(self) -> None:
        def risk(gene_id: str, symbol: str) -> dict[str, str]:
            return {"gene_id": gene_id, "gene_symbol": symbol, "evidence_type": "gwas_association",
                    "source": "FIXTURE_GWAS", "n_evidence": "1", "evidence_ids": "R", "notes": "NA"}

        risk_rows = [risk("ENSG00000000001.5", "WNT4"), risk("ENSG00000000003.1", "GREB1")]
        env_rows, _ = gs.dedupe_environment([
            env_row("E1", "ENSG00000000001.5", "WNT4", "PMID:1"),
            env_row("E2", "ENSG00000000001.5", "WNT4", "PMID:2", tissue="liver", source="OTHER"),
            env_row("E3", "ENSG00000000008.4", "IL15", "PMID:3"),
        ], min_publications=1)

        risk_path, env_path = self.temp / "risk.tsv", self.temp / "env.tsv"
        gs.write_tsv(risk_path, gs.RISK_MEMBERSHIP_COLUMNS, risk_rows)
        gs.write_tsv(env_path, gs.ENVIRONMENT_MEMBERSHIP_COLUMNS, env_rows)
        overlap_path, cand_path = self.temp / "overlap.tsv", self.temp / "candidate_genes.tsv"

        code, out, err = self.run_cli("overlap", "--risk", str(risk_path), "--environment", str(env_path),
                                      "--overlap", str(overlap_path), "--candidates", str(cand_path))
        self.assertEqual(code, 0, err)
        for expected in ('"risk_only": 1', '"environment_only": 1', '"overlap": 1', '"union": 3'):
            self.assertIn(expected, out)

        _, overlap = gs.read_tsv(overlap_path)
        self.assertEqual(overlap, [{"gene_id": "ENSG00000000001.5", "gene_symbol": "WNT4",
                                    "risk_source": "FIXTURE_GWAS",
                                    "environment_source": "FIXTURE_CTD|OTHER"}])
        header, candidates = gs.read_tsv(cand_path)
        self.assertEqual(gs.check_table("candidates", header, candidates), [])
        self.assertEqual([c["category"] for c in candidates],
                         ["overlap", "risk_only", "environment_only"])

    def test_refuses_invalid_membership(self) -> None:
        risk_path, env_path = self.temp / "risk.tsv", self.temp / "env.tsv"
        risk_path.write_text("gene_id\tgene_symbol\nENSG00000000001.5\tWNT4\n")
        gs.write_tsv(env_path, gs.ENVIRONMENT_MEMBERSHIP_COLUMNS, [])
        code, _, err = self.run_cli("overlap", "--risk", str(risk_path), "--environment", str(env_path),
                                    "--overlap", str(self.temp / "o.tsv"),
                                    "--candidates", str(self.temp / "c.tsv"))
        self.assertEqual(code, 1)
        self.assertIn("header does not match", err)
        self.assertFalse((self.temp / "c.tsv").exists())


if __name__ == "__main__":
    unittest.main()
