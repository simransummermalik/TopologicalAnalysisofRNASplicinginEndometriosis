#!/usr/bin/env python3
"""Apply the recorded E6 context decisions to the CTD environment evidence.

The decisions themselves are in ``docs/gene_set_methods.md`` (E6 context
review, 2026-11-17). This script only applies them, so a rerun after a fresh
extraction gives the same table.

Subcommands:

    tag-publications  Label each cited publication by tissue from its PubMed
                      MeSH headings and title. Reads saved PubMed XML
                      (``efetch`` output); fetches it first with --fetch.
    apply             Write the paper-level tissue label into ``tissue`` and
                      exclude co-treatment rows with a recorded reason.

Tissue labels describe the cited paper, not the individual measurement: a
paper indexed under Endometrium may also report liver cells. That is why the
value is prefixed ``cited_paper:``.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gene_sets as gs  # noqa: E402

# Ordered from most to least relevant to the study tissue (human endometrium).
TISSUE_TIERS = ["endometrium_uterus", "other_female_reproductive", "other_or_unstated"]
TIER_MESH = {
    "endometrium_uterus": {
        "Endometrium", "Endometriosis", "Endometrial Neoplasms", "Endometrial Hyperplasia",
        "Uterus", "Uterine Neoplasms", "Uterine Diseases", "Decidua", "Myometrium",
        "Leiomyoma", "Adenomyosis", "Menstrual Cycle",
    },
    # "Pregnancy" and "Oocytes" are left out on purpose: in this corpus they
    # tag fetal-tissue, eye, allergy, and Xenopus-oocyte ion-channel studies.
    "other_female_reproductive": {
        "Ovary", "Ovarian Neoplasms", "Granulosa Cells", "Theca Cells",
        "Corpus Luteum", "Luteal Cells", "Fallopian Tubes", "Cervix Uteri", "Vagina",
        "Placenta", "Trophoblasts", "Chorion", "Amnion",
    },
}
TIER_TITLE = {
    "endometrium_uterus": re.compile(r"endometri|uter|ishikawa|decidu|myometri", re.I),
    "other_female_reproductive": re.compile(
        r"ovar|granulosa|theca|luteal|placent|trophoblast|chorio|cervi", re.I
    ),
}
# Papers whose MeSH tag points to the wrong human tissue, checked against the
# abstract during the E6 review (2026-11-17).
OVERRIDES = {
    "16327978": ("other_female_reproductive",
                 "BeWo choriocarcinoma (placental) cells; MeSH Uterine Neoplasms is the tumor site"),
    "16489053": ("other_or_unstated",
                 "human work in MCF-7 breast cancer cells; MeSH Uterus is from the mouse experiment"),
}
TISSUE_PREFIX = "cited_paper:"
PUBLICATION_COLUMNS = ["pmid", "tissue_tier", "matched", "title"]
E6_EXCLUSION = (
    "E6 excluded 2026-11-17: co-treatment, so the response cannot be credited to "
    "the dioxin-like compound alone"
)
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def classify_publication(mesh: list[str], title: str) -> tuple[str, str]:
    """Return (tier, the MeSH headings or title words that decided it)."""
    for tier in TISSUE_TIERS[:-1]:
        hits = sorted(set(mesh) & TIER_MESH[tier])
        title_hits = sorted({m.lower() for m in TIER_TITLE[tier].findall(title)})
        if hits or title_hits:
            matched = [f"MeSH:{h}" for h in hits] + [f"title:{t}" for t in title_hits]
            return tier, gs.LIST_SEPARATOR.join(matched)
    return TISSUE_TIERS[-1], gs.MISSING


def parse_pubmed_xml(text: str) -> dict[str, tuple[list[str], str]]:
    root = ET.fromstring(text)
    records = {}
    for article in root.iter("PubmedArticle"):
        pmid = article.findtext("MedlineCitation/PMID")
        mesh = [d.text or "" for d in article.iter("DescriptorName")]
        title_element = article.find("MedlineCitation/Article/ArticleTitle")
        title = "".join(title_element.itertext()) if title_element is not None else ""
        records[pmid] = (mesh, " ".join(title.split()))
    return records


def cited_pmids(rows: list[dict[str, str]]) -> list[str]:
    pmids = {
        p.removeprefix("PMID:")
        for r in rows
        for p in r["publications"].split(gs.LIST_SEPARATOR)
        if p.startswith("PMID:")
    }
    return sorted(pmids, key=int)


def fetch_pubmed(pmids: list[str], cache_dir: Path, batch: int = 200) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    for start in range(0, len(pmids), batch):
        out = cache_dir / f"pubmed_{start // batch:03d}.xml"
        if out.exists():
            continue
        query = urllib.parse.urlencode(
            {"db": "pubmed", "retmode": "xml", "tool": "gene_set_review",
             "id": ",".join(pmids[start:start + batch])}
        )
        with urllib.request.urlopen(f"{EUTILS}?{query}", timeout=120) as response:
            out.write_bytes(response.read())
        time.sleep(1)  # NCBI asks for at most 3 requests per second without a key


def command_tag(args: argparse.Namespace) -> int:
    _, rows = gs.read_tsv(args.evidence)
    pmids = cited_pmids(rows)
    if args.fetch:
        fetch_pubmed(pmids, args.pubmed_dir)
    records: dict[str, tuple[list[str], str]] = {}
    for path in sorted(args.pubmed_dir.glob("*.xml")):
        records.update(parse_pubmed_xml(path.read_text(encoding="utf-8")))
    missing = [p for p in pmids if p not in records]
    if missing:
        raise ValueError(f"{len(missing)} cited PMIDs have no PubMed record, e.g. {missing[:5]}")
    out_rows = []
    for pmid in pmids:
        mesh, title = records[pmid]
        if pmid in OVERRIDES:
            tier, reason = OVERRIDES[pmid]
            matched = f"override: {reason}"
        else:
            tier, matched = classify_publication(mesh, title)
        out_rows.append({"pmid": pmid, "tissue_tier": tier, "matched": matched,
                         "title": title or gs.MISSING})
    gs.write_tsv(args.output, PUBLICATION_COLUMNS, out_rows)
    print(json.dumps({"publications": len(out_rows),
                      "tissue_tier": dict(Counter(r["tissue_tier"] for r in out_rows))}, indent=2))
    return 0


def apply_review(
    rows: list[dict[str, str]], tiers: dict[str, str]
) -> tuple[list[dict[str, str]], dict[str, object]]:
    """Return reviewed rows and counts. Safe to rerun on its own output."""
    counts: Counter = Counter()
    tissue_counts: Counter = Counter()
    out = []
    for row in rows:
        row = dict(row)
        pmids = [p.removeprefix("PMID:") for p in row["publications"].split(gs.LIST_SEPARATOR)]
        try:
            row_tiers = [tiers[p] for p in pmids]
        except KeyError as error:
            raise ValueError(f"{row['evidence_id']}: PMID {error.args[0]} is not tagged") from None
        tier = min(row_tiers, key=TISSUE_TIERS.index)
        row["tissue"] = TISSUE_PREFIX + tier
        tissue_counts[tier] += 1

        if "cotreatment" in row["context"].split(gs.LIST_SEPARATOR) and E6_EXCLUSION not in row["notes"]:
            if row["mapping_status"] == "mapped":
                note = f"{E6_EXCLUSION} (was mapped to {row['gene_id']} {row['gene_symbol']})"
            else:
                note = f"{E6_EXCLUSION} (was {row['mapping_status']})"
            row["notes"] = f"{row['notes']}; {note}"
            row["mapping_status"] = "excluded"
            row["gene_id"] = gs.MISSING
            row["gene_symbol"] = gs.MISSING
            counts["newly_excluded_cotreatment"] += 1
        counts[row["mapping_status"]] += 1
        out.append(row)
    return out, {"mapping_status": {k: v for k, v in counts.items() if k != "newly_excluded_cotreatment"},
                 "newly_excluded_cotreatment": counts["newly_excluded_cotreatment"],
                 "tissue_tier": dict(tissue_counts)}


def command_apply(args: argparse.Namespace) -> int:
    header, rows = gs.read_tsv(args.evidence)
    if header != gs.ENVIRONMENT_EVIDENCE_COLUMNS:
        raise ValueError(f"{args.evidence}: not an environment_evidence table")
    _, pubs = gs.read_tsv(args.publications)
    tiers = {r["pmid"]: r["tissue_tier"] for r in pubs}
    reviewed, counts = apply_review(rows, tiers)
    gs.write_tsv(args.output, gs.ENVIRONMENT_EVIDENCE_COLUMNS, reviewed)
    print(json.dumps(counts, indent=2))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    tag = sub.add_parser("tag-publications", help="label cited papers by tissue")
    tag.add_argument("--evidence", type=Path, required=True)
    tag.add_argument("--pubmed-dir", type=Path, required=True, help="folder of PubMed efetch XML")
    tag.add_argument("--fetch", action="store_true", help="download missing PubMed XML first")
    tag.add_argument("--output", type=Path, required=True)
    tag.set_defaults(func=command_tag)
    apply = sub.add_parser("apply", help="apply the E6 decisions")
    apply.add_argument("--evidence", type=Path, required=True)
    apply.add_argument("--publications", type=Path, required=True)
    apply.add_argument("--output", type=Path, required=True)
    apply.set_defaults(func=command_apply)
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (OSError, ValueError, ET.ParseError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
