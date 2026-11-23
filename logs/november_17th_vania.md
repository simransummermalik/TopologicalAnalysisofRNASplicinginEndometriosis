# November 17, 2026 — Vania

**Status:** Pushed

**Commits:**

- `3e3a0e4` — Draft CTD environmental evidence (E3–E4)
- `548b76a` — Complete E6 context review of CTD environmental evidence
- `1389da1` — E5–E8 work
- `8a825a7` — Environmental interpretation boundary (E11)
- `584466b` — Draft genetic-risk evidence (G4–G6)
- `9b660b8` — Resolving genetic risk rows (G7)
- `28c9987` — Genetic-risk set work (G8–G10)
- `bfabb0f` — Wrapping up priority 2 work (E9–E10)
- `261a0c4` — Finalized the environmental publication rule
- `7486bb0`, `68d2992` — Spot-check sheet and task-board link
- `f08fc0c` — Priority 4 (S5, S7, S8)

## Contribution

I extracted and reviewed the candidate-gene evidence (Priority 2) and fixed
the statistical inference rules (Priority 4). Both were done before any real
splicing result existed.

### Environmental-response evidence (E3–E11)

- Extracted CTD curated chemical–gene interactions for TCDD and the fixed list
  of dioxin-like compounds, *Homo sapiens* only, measured responses only
  (`gene_sets/extract_ctd.py`).
- Context review (E6): excluded 1,028 co-treatment rows, because a mixture
  response cannot be credited to the dioxin alone; kept and labeled rows where
  another agent modifies the dioxin's effect or the reverse
  (`gene_sets/review_ctd_context.py`).
- Tagged each row with the cited paper's tissue from PubMed MeSH headings and
  titles. Tissue is never a reason to exclude; 98% of kept rows come from
  non-reproductive tissue, and the tag keeps that limit visible.
- Mapping review and audit sample (E5, E7, E8; `gene_sets/audit_environment.py`).
- Final publication rule, replacing the provisional one: a gene–chemical pair
  needs at least 2 publications, and at least 1 must not be a large screen
  (a publication behind 100 or more evidence rows). Under the provisional rule
  64% of qualifying pairs rested only on genome-wide screens. Result:
  565 environmental genes.
- Wrote the interpretation boundary for environmental genes (E11) in the gene
  set methods and presentation notes: membership does not show exposure.

### Genetic-risk evidence (G4–G10)

- Extracted GWAS Catalog endometriosis associations and the Koller et al. 2026
  Supplementary Table 4 loci (`gene_sets/extract_gwas.py`).
- Excluded cross-trait (75 rows) and MTAG (92 rows) analyses, and Koller
  analyses other than the cross-ancestry combined definition (139 rows).
- G7 ambiguity rules: multi-gene loci keep every author-listed gene, flagged;
  retired clone names are translated through GENCODE v19 IDs; flanking-gene
  and no-gene rows are excluded rather than guessed.
- Built the deduplicated risk set and its audit sample (G8–G10,
  `gene_sets/audit_risk.py`): 103 risk genes.

### Candidate set (E9–E10)

- Combined both sets into `data/gene_sets/candidate_genes.tsv`: 661 candidate
  genes, 7 in both sets (`overlap.tsv`).
- Wrote `docs/spot_checks.md`, a worksheet for a person to check sampled rows
  against their sources. It still needs a reviewer.

### Statistics rules (S5, S7, S8)

- Decided the inference rules: every test two-sided; a gene is testable only with
  all 4 cases and 4 controls; the single confirmatory test is one gene-set
  test on the 661 candidates at α = 0.05; risk, environment and overlap set
  tests and gene-level BH are exploratory. Recorded in `docs/decisions.md` and
  `docs/statistics_spec.md` (v0.5).
- Implemented the gene-set permutation test in Rust
  (`set_permutation_test` in `src/permutation.rs`, `test_gene_set` in
  `src/statistics.rs`): the mean over genes of the absolute case-minus-control
  cycle-fraction effect, with each of the 70 relabelings applied to every gene
  at once. Added `set_test.tsv`, `set_scores.tsv`, a report section and tests.

## Current scope

The evidence tables and statistics code are complete, but the spot checks and
the G1 cold read still need a person. No real case/control result was
produced.
