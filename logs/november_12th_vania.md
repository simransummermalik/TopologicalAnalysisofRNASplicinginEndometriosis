# November 12, 2026 — Vania

**Status:** Pushed

**Commit:** `0441b28` — `Priority 2 work`

## Contribution

I set up Priority 2 (candidate-gene evidence, tasks G1–G10 and E1–E11).
I wrote the evidence rules, approved and verified the sources, wrote the
dioxin-like chemical list, and built and tested the tooling that turns
evidence tables into gene sets. No real evidence rows were extracted yet.
Choosing which evidence counts is still a human task under the new rules.

### Rules and decisions (G1, G2, E1, E3)

New file `docs/gene_set_methods.md`, with dated entries in `docs/decisions.md`:

- **Join key.** Versioned GENCODE v50 `gene_id`, the same as the junction
  tables. The NCBI GRCh38.p14 records in `human-reference-genome/` are marked
  as not canonical.
- **Genetic risk (G1).** Endometriosis associations at P < 5 × 10⁻⁸, using the
  gene the source reports or maps. Loci without a single gene are marked
  ambiguous, never guessed. Expression-only evidence never counts as risk.
- **Risk sources (G2).** The GWAS Catalog, plus Koller et al., *Nature
  Genetics* 58(5):1051–1061 (2026). Koller replaces Rahmioglu et al. 2023
  (42 loci) as the largest meta-analysis. Only its cross-ancestry
  combined-definition analysis is used: 105,869 cases, 80 loci, from
  Supplementary Table 4.
- **Exposure scope (E1).** Dioxin-like aryl-hydrocarbon-receptor (AhR)
  ligands.
- **Environmental source (E3).** CTD curated chemical–gene interactions, with
  these limits:
  - *Homo sapiens* only
  - measured expression, protein, or activity responses only
  - a provisional rule of at least 2 publications per gene–chemical pair,
    to be confirmed after the E10 counts
- **Literature genes.** The expression-only genes from the group research PDF
  (for example SULT1E1, NNMT, and BCL2L11 reported as "BIM") go into a
  separate literature-context list. They are never counted as evidence.

### Sources checked against the original papers

- **WHO 2005 TEF list.** Checked against Table 1 of Van den Berg et al.,
  *Toxicol Sci* 93(2):223–241 (2006), doi:10.1093/toxsci/kfl055. It lists
  exactly 29 congeners: 7 PCDDs, 10 PCDFs, and 12 dioxin-like PCBs.
- **CTD chemical IDs.** Matched by name and synonym in the CTD chemical
  vocabulary (report dated 2026-10-29; checksum recorded). 21 of the 29
  congeners have CTD IDs; TCDD is `MESH:D013749`. The 8 without IDs are
  listed. Class and mixture terms are explicitly excluded, because
  `Polychlorinated Dibenzodioxins` lists TCDD as a synonym.
- **Koller 2026.** Volume and pages are from Crossref. The threshold, primary
  analysis, locus table, and Zenodo summary statistics are from the
  open-access accepted manuscript (https://hdl.handle.net/2445/229292).
- **2022 WHO update** (DeVito et al. 2024). The publisher blocks full-text
  download, so its congener table was not checked. The definition stays
  anchored to the verified 2005 table, recorded as a decision.
- **HGNC export.** Recorded in `docs/data_methods.md` with the download date
  and MD5/SHA-256, in the same format as the GENCODE files.

### Tooling (G3, G6–G8, G10, E2, E5, E7–E10)

- `gene_sets/gene_sets.py` uses the Python standard library only, like
  `preprocessing/`. It has four subcommands:
  - `map`: original labels → versioned GENCODE IDs. It tries the GTF ID, then
    the exact GTF name, then the HGNC approved, previous, and alias symbols.
    It only touches `pending` rows, so manual review decisions survive
    reruns. Multi-gene labels, shared names, and shared aliases are marked
    `ambiguous`, with the candidates listed. Nothing is guessed.
  - `dedupe`: one membership row per gene (risk) or per gene, chemical,
    tissue, and context (environment). Every evidence ID is kept. The
    publication threshold is an option.
  - `check`: schema validation for all six gene-set tables. It checks the
    header, blank cells, the ID format and that IDs exist in the GTF,
    symbol/ID separation, allowed values, duplicates, and consistent overlap
    flags.
  - `overlap`: writes `overlap.tsv` and `candidate_genes.tsv` with
    risk-only, environment-only, and overlap categories and counts.
- `data/gene_sets/` has header-only evidence templates and a README covering
  columns, allowed values, and the commands to run.
- `gene_sets/tests/test_gene_sets.py` has 22 unit tests on made-up genes.
  CI now runs them.

### Verification on real data

The tools were run end to end on these files:

- the GENCODE v50 GTF, whose MD5 matches `docs/data_methods.md`
- the HGNC export
- 535 GWAS Catalog endometriosis associations at P < 5 × 10⁻⁸
- 8,358 human CTD interactions for the 21 chemicals

These were converted by a throwaway script, not the real extraction. Results:

- Every output passed `check`, including the check that each ID exists in the
  annotation.
- Environmental mappings were independently cross-checked: all 8,331
  checkable rows match HGNC's link from CTD's NCBI gene ID to Ensembl, with
  0 disagreements.
- The real data exposed three problems, all fixed before the commit:
  - `check` rejected GWAS loci that have no gene, which the rules say to
    record
  - 20 names are duplicated in GENCODE v50 (for example PDE4C); these are now
    resolved by HGNC's link only when it picks exactly one ID
  - the candidate list for multi-copy names such as `Y_RNA` (759 IDs) is now
    capped
- Correctly left for people to decide: COX2 (PTGS2 or MT-CO2), DOC1, and
  alternate-haplotype genes such as HLA-DRB4.
- CI passed on the pushed commit.

Early numbers from the test data (not results):

- About a third of the risk rows list two flanking genes and will need G7
  review.
- TCDD makes up 90% of the environmental rows.
- With the 2-publication rule the environmental set has 1,706 genes, about
  17 times the 102-gene risk set. Without the rule it has 4,417 genes.

## Current scope

The setup work for G1–G3, G6, G8, G10, E1–E3, E5, and E7–E10 is in place.
Still to do:

- a cold read of the rules by a teammate (G1 step 7)
- inspecting the gene columns of Koller's Supplementary Table 4
- the extraction itself (G4–G5, E3–E4)
- human review of ambiguous mappings and context (G7, E6)
- the independent audits (G9, E8)
- the final gene sets
- the interpretation-boundary paragraph (E11)
- the literature-context list
