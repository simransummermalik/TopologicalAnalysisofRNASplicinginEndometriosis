#!/usr/bin/env python3
"""Draft environment-evidence rows from CTD chemical-gene interactions (E3-E4).

Applies the frozen E3 filters from ``docs/gene_set_methods.md`` in order and
writes one ``pending`` row per measured response, ready for ``gene_sets.py
map``. It does not apply the publication rule (``dedupe`` does that) and it
does not decide which context counts: every kept row is labeled so the E6
reviewers can decide.

Filters, in order:

1. chemical  - ChemicalID is one of the fixed dioxin-like IDs (exact IDs only,
               never names).
2. species   - Organism is Homo sapiens.
3. measured  - InteractionActions include an expression or activity change.
4. target    - The row's own gene is the one whose expression or activity was
               measured. CTD also files an interaction under genes that only
               take part in it, e.g. AHR in "AHR affects the reaction [TCDD
               increases expression of CYP1A1]"; those rows are dropped.

Kept rows get a ``context`` label:

    direct                  the chemical alone changes the gene
    cotreatment             the chemical is part of a co-treatment (mixture)
    dioxin_effect_modified  another agent alters the chemical's own effect
    dioxin_modifies_other   the chemical alters another agent's effect

Labels combine with "|" (for example "cotreatment|dioxin_effect_modified").

CTD does not record tissue, dose, duration, or assay, so those are NA and the
note points to the interaction text and cited publications.

Usage:

    python3 gene_sets/extract_ctd.py \\
      --ctd CTD_chem_gene_ixns.tsv.gz \\
      --chemicals data/gene_sets/dioxin_like_chemicals.tsv \\
      --output data/gene_sets/environment_evidence.tsv
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gene_sets as gs  # noqa: E402

CTD_FIELDS = [
    "ChemicalName",
    "ChemicalID",
    "CasRN",
    "GeneSymbol",
    "GeneID",
    "GeneForms",
    "Organism",
    "OrganismID",
    "Interaction",
    "InteractionActions",
    "PubMedIDs",
]
SPECIES = "Homo sapiens"
MEASURED_ACTION = re.compile(r"\^(expression|activity)$")
FORMS = ("mRNA", "protein", "promoter", "enhancer", "gene")
REPORT_CREATED = re.compile(r"^#\s*Report created:\s*(.+?)\s*$")
CONTEXT_NOTE = (
    "CTD does not record tissue, dose, duration, or assay; "
    "see the interaction text and cited publications"
)


def read_chemicals(path: Path) -> dict[str, str]:
    """Return {CTD ID without 'MESH:': congener} from the frozen list."""
    header, rows = gs.read_tsv(path)
    if header[:4] != ["congener", "class", "ctd_name", "ctd_chemical_id"]:
        raise ValueError(f"{path}: unexpected header {header}")
    chemicals = {}
    for row in rows:
        ctd_id = row["ctd_chemical_id"]
        if not ctd_id.startswith("MESH:"):
            raise ValueError(f"{path}: chemical ID {ctd_id!r} lacks the MESH: prefix")
        chemicals[ctd_id.removeprefix("MESH:")] = row["congener"]
    return chemicals


def measured_responses(symbol: str, interaction: str) -> list[tuple[str, bool]]:
    """Response types measured for ``symbol`` itself, with a form-not-stated flag.

    Matches "expression of SYMBOL mRNA", "activity of SYMBOL protein", and
    CTD's joined phrasing "expression of and affects the secretion of SYMBOL
    protein". Returns [] when the gene only takes part in the interaction.
    """
    pattern = re.compile(
        r"\b(expression|activity) of (?:and [a-z ]+? of )*"
        + re.escape(symbol)
        + r"(?: (" + "|".join(FORMS) + r"))?(?=[\s\]]|$)"
    )
    found: dict[str, bool] = {}
    for kind, form in pattern.findall(interaction):
        if kind == "activity":
            response = "activity"
        elif form == "protein":
            response = "protein"
        else:
            response = "expression"
        no_form = form == ""
        found[response] = found.get(response, True) and no_form
    return sorted(found.items())


def classify_context(actions: list[str], chemical: str, interaction: str) -> str:
    """Label how the chemical relates to the measured change.

    A reaction is split by where the chemical sits: outside the brackets it
    alters another agent's effect ("TCDD inhibits [estradiol increases X]");
    inside, its own effect is altered by another agent ("AHR affects [TCDD
    increases X]").
    """
    kinds = {action.split("^", 1)[-1] for action in actions}
    labels = []
    if "cotreatment" in kinds:
        labels.append("cotreatment")
    if "reaction" in kinds:
        outer = re.sub(r"\[.*\]", "[...]", interaction)
        labels.append("dioxin_modifies_other" if chemical in outer else "dioxin_effect_modified")
    return gs.LIST_SEPARATOR.join(labels) if labels else "direct"


def extract(lines, chemicals: dict[str, str], source_name: str):
    """Return (evidence rows, counts after each filter and per chemical)."""
    counts: dict[str, object] = {}
    steps = ["input", "chemical", "species", "measured", "target"]
    step_counts = Counter()
    per_chemical: dict[str, Counter] = defaultdict(Counter)
    by_context = Counter()
    by_response = Counter()
    source_version = None
    rows = []

    for line_number, raw in enumerate(lines, start=1):
        line = raw.rstrip("\n")
        if line.startswith("#"):
            match = REPORT_CREATED.match(line)
            if match:
                source_version = match.group(1)
            continue
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != len(CTD_FIELDS):
            raise ValueError(
                f"{source_name} line {line_number}: expected {len(CTD_FIELDS)} "
                f"fields, found {len(fields)}"
            )
        record = dict(zip(CTD_FIELDS, fields))
        step_counts["input"] += 1

        chemical_id = record["ChemicalID"]
        if chemical_id not in chemicals:
            continue
        step_counts["chemical"] += 1
        per_chemical[chemical_id]["chemical"] += 1

        if record["Organism"] != SPECIES:
            continue
        step_counts["species"] += 1
        per_chemical[chemical_id]["species"] += 1

        actions = [a for a in record["InteractionActions"].split("|") if a]
        if not any(MEASURED_ACTION.search(a) for a in actions):
            continue
        step_counts["measured"] += 1
        per_chemical[chemical_id]["measured"] += 1

        responses = measured_responses(record["GeneSymbol"], record["Interaction"])
        if not responses:
            continue
        step_counts["target"] += 1
        per_chemical[chemical_id]["target"] += 1

        context = classify_context(actions, record["ChemicalName"], record["Interaction"])
        by_context[context] += 1
        pmids = [p for p in record["PubMedIDs"].split("|") if p]
        publications = gs.LIST_SEPARATOR.join(f"PMID:{p}" for p in pmids) or gs.MISSING
        for response, no_form in responses:
            by_response[response] += 1
            notes = [
                f"CTD GeneID {record['GeneID'] or gs.MISSING}",
                f"interaction: {record['Interaction']}",
                CONTEXT_NOTE,
            ]
            if no_form:
                notes.append("CTD does not state mRNA or protein")
            rows.append(
                {
                    "evidence_id": f"CTD:{line_number}:{response}",
                    "source": "CTD",
                    "source_version": f"CTD_chem_gene_ixns report {source_version or gs.MISSING}",
                    "source_record": f"{source_name} line {line_number}",
                    "publications": publications,
                    "original_id": record["GeneSymbol"],
                    "exposure": record["ChemicalName"],
                    "exposure_id": f"MESH:{chemical_id}",
                    "dose_duration": gs.MISSING,
                    "species": SPECIES,
                    "tissue": gs.MISSING,
                    "assay": gs.MISSING,
                    "response_type": response,
                    "evidence_type": "measured_response",
                    "context": context,
                    "gene_id": gs.MISSING,
                    "gene_symbol": gs.MISSING,
                    "mapping_status": "pending",
                    "mapping_method": "not_run",
                    "notes": "; ".join(notes),
                }
            )

    if source_version is None:
        raise ValueError(f"{source_name}: no 'Report created' line; is this a CTD report?")
    counts["source_version"] = source_version
    counts["rows_after_each_filter"] = {step: step_counts[step] for step in steps}
    counts["per_chemical"] = {
        f"MESH:{cid} ({name})": {step: per_chemical[cid][step] for step in steps[1:]}
        for cid, name in chemicals.items()
    }
    counts["evidence_rows"] = len(rows)
    counts["kept_interactions_by_context"] = dict(sorted(by_context.items()))
    counts["evidence_rows_by_response_type"] = dict(sorted(by_response.items()))
    return rows, counts


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--ctd", type=Path, required=True, help="CTD_chem_gene_ixns.tsv(.gz)")
    parser.add_argument("--chemicals", type=Path, required=True, help="frozen chemical list TSV")
    parser.add_argument("--output", type=Path, required=True, help="environment_evidence.tsv to write")
    args = parser.parse_args(argv)
    try:
        chemicals = read_chemicals(args.chemicals)
        with gs.open_text(args.ctd) as handle:
            rows, counts = extract(handle, chemicals, args.ctd.name)
        gs.write_tsv(args.output, gs.ENVIRONMENT_EVIDENCE_COLUMNS, rows)
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    print(json.dumps(counts, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
