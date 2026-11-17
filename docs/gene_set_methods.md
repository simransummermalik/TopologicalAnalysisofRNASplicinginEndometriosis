# Candidate-gene set methods

This file holds the written rules for Priority 2 (tasks G1–G10 and E1–E11).
A reader should be able to classify a new evidence row using only these rules.
Dated decisions are mirrored in [`docs/decisions.md`](decisions.md). Do not
rewrite a rule silently; add a dated revision instead.

Status: **draft rules recorded 2026-11-12; awaiting independent read-through
(G1 step 7).**

## Stable identifiers

- Join key: the versioned GENCODE v50 `gene_id` (for example
  `ENSG00000123456.7`), the same namespace the junction tables use
  ([`docs/data_contract.md`](data_contract.md)).
- Mapping resource: the GENCODE v50 comprehensive GTF recorded in
  [`docs/data_methods.md`](data_methods.md), plus a dated HGNC export for
  previous and alias symbols.
- The original identifier from the source is always kept in its own column.
  The readable symbol is kept separately and is never used as a join key.

## Genetic-susceptibility rules (G1)

**Definition.** For this project, a gene shows genetic susceptibility to
endometriosis when an inherited DNA variant associated with endometriosis
reaches genome-wide significance and the source links that variant to the gene.

**Qualifying evidence.**
- Endometriosis association at p < 5×10⁻⁸ in a genome-wide association study
  or meta-analysis.
- The gene is the one the source reports or maps for that association.

**Does not qualify as inherited risk.**
- Expression differences between cases and controls (bulk, single-cell, or
  lesion versus endometrium).
- Associations below genome-wide significance.
- Genes inferred by the extractor from proximity alone when the source names
  no gene.

**Recording level.** Each evidence row records the variant or locus exactly as
the source gives it, its p-value, the reported/mapped gene label, and
`evidence_level` = `variant`, `locus`, or `gene`.

**Ambiguous evidence.** A locus whose source lists several genes, or no gene,
is recorded with `mapping_status = ambiguous` and resolved or excluded in G7
with a written reason. It is never silently assigned to the nearest gene.

**Exclusions.** Rows failing the threshold, non-endometriosis traits, and
non-human studies are excluded with a reason in the notes column.

## Approved risk sources (G2)

| Source | Access location | Version / release | Accessed | Status |
|---|---|---|---|---|
| GWAS Catalog, endometriosis trait | https://www.ebi.ac.uk/gwas/ | Record release date at extraction | — | Approved, not yet extracted |
| Largest endometriosis GWAS meta-analysis: Koller, Cormand et al., *Nat Genet* 2026, doi:10.1038/s41588-026-02582-2 | Supplementary Table 4 (cross-ancestry combined definition only) | *Nat Genet* 58(5):1051–1061, 2026; record supplement version at extraction | — | Approved 2026-11-12; full-text checks done (below) |

**Meta-analysis candidates (checked 2026-11-12).**

- *Verified:* Rahmioglu N, …, Zondervan KT. The genetic basis of endometriosis
  and comorbidity with other pain and inflammatory conditions. *Nature
  Genetics* 55(3):423–436 (2023). doi:10.1038/s41588-023-01323-z. 60,674 cases
  and 701,926 controls of European and East Asian ancestry; 42 genome-wide
  significant loci (49 signals).
- *Verified 2026-11-12:* Koller D, He J, Lokhammer S, …, Cormand B, …,
  Polimanti R. Multi-ancestry genome-wide association and integrated
  multi-omics analyses of endometriosis and its clinical manifestations.
  *Nature Genetics* 58(5):1051–1061 (2026). doi:10.1038/s41588-026-02582-2.
  Bibliographic details from Crossref. Study details from the open-access
  accepted manuscript in the University of Barcelona repository
  (https://hdl.handle.net/2445/229292).

**Decision (2026-11-12):** use Koller et al. 2026, recorded in
[`docs/decisions.md`](decisions.md). Rahmioglu et al. 2023 is kept here only as
the superseded comparison and is not an extraction source.

**Full-text checks (done 2026-11-12, from the accepted manuscript):**

1. **Threshold.** The paper uses P < 5 × 10⁻⁸ ("based on Bonferroni multiple
   testing correction"), the same as the G1 rule. Conditional (COJO)
   secondary signals use the same threshold on the joint-analysis P.
2. **Analysis to extract.** The 80 loci, 37 of them new, come from the
   *cross-ancestry endometriosis combined definition* meta-analysis: 8
   cohorts, 105,869 cases, and 1,282,731 female controls, with cases from
   self-report and electronic health records. Extract this analysis only.
   Excluded: ancestry-specific analyses (for example 63 loci in EUR
   combined), EUR clinical, EUR self-reported, adenomyosis, and
   endometriosis without adenomyosis. Loci found only in those analyses do
   not enter the risk set.
3. **Where the loci are.** Main-text Fig. 2 and Supplementary Table 4, which
   lists genome-wide significant loci across definitions. The supplement is
   not part of the accepted manuscript, so its gene columns are still to be
   inspected at G4. The GWAS Catalog has curated this paper (first author
   "Koller D", PubMed 42056605, 219 associations at P < 5 × 10⁻⁸ across its
   analyses), which gives a second, cross-checkable copy of the loci with
   mapped genes.
4. **Summary statistics.** Zenodo, doi:10.5281/zenodo.18983492.

A source is described as used only after its evidence rows are saved.

## Environmental-response scope (E1)

**Exposure.** Dioxin-like aryl-hydrocarbon-receptor (AhR) ligands: TCDD and
other dioxin-like compounds. This is a chemical class, not a single chemical.

**Relevant evidence.** A measured gene-level response (expression, protein, or
activity change) to a dioxin-like compound in an experimental or observational
study.

**Outside scope.** Other endocrine disruptors (for example BPA or phthalates),
general estrogen response without a dioxin-like exposure, and
disease-versus-control expression findings.

**Context.** Species, tissue, assay, dose, and duration are recorded for every
row; mismatches are labeled and reviewed in E6.

Scope approved: 2026-11-12, recorded by Vania Zeledon. Changing this scope
requires a new dated decision.

## Approved environmental source (E3)

| Source | Access location | Version / release | Accessed | Status |
|---|---|---|---|---|
| CTD curated chemical–gene interactions | https://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz | Report created 2026-10-29 | 2026-11-17 | Extracted as a draft; awaiting E6 review and E8 audit (see the extraction record below) |

**Filters, applied in this order and recorded with row counts after each:**

1. **Chemicals.** TCDD plus a fixed list of dioxin-like compounds, written
   below with CTD chemical identifiers *before* any genes are extracted. The
   list is not changed after seeing results without a new dated decision.
2. **Species.** *Homo sapiens* only. Rodent or other species may enter later
   only as a labeled indirect tier through a new dated decision.
3. **Interaction type.** Measured responses only: changes in expression,
   protein, or activity. Binding-only or co-mention interactions are excluded.
4. **Support (provisional).** At least 2 distinct supporting publications per
   gene–chemical pair. Confirm or revise after the E10 category counts.

### Dioxin-like chemical list

Written 2026-11-12, before any gene extraction.

**Inclusion rule.** A chemical is dioxin-like for this project when it is one
of the 29 congeners assigned a toxic equivalency factor (TEF) by the WHO
scheme: 7 polychlorinated dibenzo-p-dioxins (PCDDs), 10 polychlorinated
dibenzofurans (PCDFs), and 12 dioxin-like polychlorinated biphenyls (PCBs).
These are the compounds the WHO treats as acting through the AhR in the same
way as TCDD.

Reference (verified 2026-11-12 against the paper's Table 1, "Summary of WHO
1998 and WHO 2005 TEF Values"): Van den Berg M, Birnbaum LS, Denison M, et al.
The 2005 World Health Organization reevaluation of human and mammalian toxic
equivalency factors for dioxins and dioxin-like compounds. *Toxicological
Sciences* 93(2):223–241 (2006). doi:10.1093/toxsci/kfl055. Table 1 lists
exactly the 29 congeners used here: 7 PCDDs, 10 PCDFs, 4 non-ortho PCBs (77,
81, 126, 169), and 8 mono-ortho PCBs (105, 114, 118, 123, 156, 157, 167, 189).

A 2022 WHO reevaluation also exists (DeVito et al., *Regulatory Toxicology and
Pharmacology* 146, 2024; doi:10.1016/j.yrtph.2023.105525). It revised TEF
values. Its abstract says the panel reevaluated "the 2005 WHO TEFs for
chlorinated dioxin-like compounds", assigning new values "for almost all
congeners" and keeping the 2005 values for mono-ortho PCBs. That suggests the
same set of compounds, but the congener table itself was not checked: the
publisher blocks full-text download through PMC, Europe PMC, and NCBI
E-utilities (checked 2026-11-12).

**Resolution:** the project's definition is anchored to the verified 2005
table. This does not depend on the 2022 values, because the project uses the
list of congeners and not the TEFs. If someone with journal access finds a
congener added or removed in 2022, record it as a new dated decision.

**ID lookup.** CTD chemical IDs were matched by name and synonym in the CTD
chemical vocabulary:

| Field | Recorded value |
|---|---|
| File | `CTD_chemicals.tsv.gz` from https://ctdbase.org/reports/CTD_chemicals.tsv.gz |
| CTD report date | 2026-10-29 |
| Accessed | 2026-11-12 |
| SHA-256 | `3939d41b3fb9e39408b5f06f1060d60eb535c924de7fef259288626fc6298b4c` |
| Rows | 179,669 chemicals |

**Included chemicals (21 with CTD IDs).**

| Congener | Class | CTD name | CTD chemical ID |
|---|---|---|---|
| 2,3,7,8-TCDD | PCDD | Tetrachlorodibenzodioxin | MESH:D013749 |
| 1,2,3,7,8-PeCDD | PCDD | 1,2,3,7,8-pentachlorodibenzo-p-dioxin | MESH:C051452 |
| 1,2,3,4,7,8-HxCDD | PCDD | 1,2,3,4,7,8-hexachlorodibenzodioxin | MESH:C078219 |
| 1,2,3,6,7,8-HxCDD | PCDD | 1,2,3,6,7,8-hexachlorodibenzodioxin | MESH:C046838 |
| 1,2,3,4,6,7,8-HpCDD | PCDD | 1,2,3,4,6,7,8-heptachlorodibenzodioxin | MESH:C046839 |
| OCDD | PCDD | octachlorodibenzo-4-dioxin | MESH:C018117 |
| 2,3,7,8-TCDF | PCDF | 2,3,7,8-tetrachlorodibenzofuran | MESH:C014211 |
| 1,2,3,7,8-PeCDF | PCDF | 1,2,3,7,8-pentachlorodibenzofuran | MESH:C057584 |
| 2,3,4,7,8-PeCDF | PCDF | 2,3,4,7,8-pentachlorodibenzofuran | MESH:C038890 |
| 1,2,3,4,7,8-HxCDF | PCDF | 1,2,3,4,7,8-hexachlorodibenzofuran | MESH:C051412 |
| 1,2,3,4,6,7,8-HpCDF | PCDF | 1,2,3,4,6,7,8-heptachlorodibenzofuran | MESH:C473650 |
| OCDF | PCDF | octachlorodibenzofuran | MESH:C008694 |
| PCB 77 | non-ortho PCB | 3,4,3',4'-tetrachlorobiphenyl | MESH:C028451 |
| PCB 126 | non-ortho PCB | 3,4,5,3',4'-pentachlorobiphenyl | MESH:C023035 |
| PCB 169 | non-ortho PCB | 3,4,5,3',4',5'-hexachlorobiphenyl | MESH:C027309 |
| PCB 105 | mono-ortho PCB | 2,3,3',4,4'-pentachlorobiphenyl | MESH:C067940 |
| PCB 114 | mono-ortho PCB | 2,3,4,4'5-pentachlorobiphenyl *(CTD's spelling)* | MESH:C035834 |
| PCB 118 | mono-ortho PCB | 2,3',4,4',5-pentachlorobiphenyl | MESH:C070055 |
| PCB 156 | mono-ortho PCB | 2,3,3',4,4',5-hexachlorobiphenyl | MESH:C087667 |
| PCB 167 | mono-ortho PCB | 2,4,5,3',4',5'-hexachlorobiphenyl | MESH:C038372 |
| PCB 189 | mono-ortho PCB | 2,3,4,5,3',4',5'-heptachlorobiphenyl | MESH:C034559 |

PCB 114, 167, and 189 matched through CTD synonyms (`PCB 114`, `PCB 167`,
`PCB 189`) because CTD lists them under a different locant order.

**WHO congeners with no CTD entry (8).** These stay in the definition but
cannot contribute rows: 1,2,3,7,8,9-HxCDD; 1,2,3,6,7,8-HxCDF;
1,2,3,7,8,9-HxCDF; 2,3,4,6,7,8-HxCDF; 1,2,3,4,7,8,9-HpCDF; PCB 81; PCB 123;
PCB 157. Re-check them if a newer CTD release is used.

**Explicitly excluded CTD terms.** Mixture or class terms such as
`Polychlorinated Dibenzodioxins` (MESH:D000072317),
`Dibenzofurans, Polychlorinated` (MESH:D000072338), and
`Polychlorinated Biphenyls` also cover congeners that are not dioxin-like.
`Polychlorinated Dibenzodioxins` even lists 2,3,7,8-TCDD as a synonym, so
extraction must filter on the exact IDs above, never on names. Non-WHO
congeners (for example 1,2,3,4-TCDD) and substituted or labeled derivatives
are also excluded.

**Expected skew.** Human interaction data in CTD are likely to be dominated by
TCDD. Report per-chemical row counts at E3 so that skew is visible.

Source and filters approved: 2026-11-12, recorded by Vania Zeledon.

## CTD extraction record (E3–E4)

Drafted 2026-11-17 by Vania Zeledon. The rows are a machine-applied draft of
the frozen filters. They are not reviewed evidence until E6 and E8 are done.

| Field | Recorded value |
|---|---|
| File | `CTD_chem_gene_ixns.tsv.gz` from https://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz |
| CTD report date | Thu Oct 29 12:43:06 EDT 2026 (same release as the chemical vocabulary above) |
| Accessed | 2026-11-17 |
| Size | 43,514,023 bytes compressed; 3,190,383 interaction rows |
| SHA-256 | `032f1b6fc674ca8614b91efeefe40e2cab8315639c1b1552a5df4df98d88f8c9` |
| MD5 | `cd371d942b0214f2d395a59dd7db2ca0` |
| Chemical list | [`data/gene_sets/dioxin_like_chemicals.tsv`](../data/gene_sets/dioxin_like_chemicals.tsv), a copy of the 21-row table above |
| Script | [`gene_sets/extract_ctd.py`](../gene_sets/extract_ctd.py), then `gene_sets.py map` |
| Output | [`data/gene_sets/environment_evidence.tsv`](../data/gene_sets/environment_evidence.tsv) (7,958 rows) and `environment_mapping.tsv` |

The file is not kept in Git. To rebuild, download the same release, check the
SHA-256, and run the commands in
[`data/gene_sets/README.md`](../data/gene_sets/README.md).

### How the filters were applied

1. **Chemicals.** Exact match of `ChemicalID` to the 21 IDs. Names are never
   used.
2. **Species.** `Organism` is exactly `Homo sapiens`.
3. **Interaction type.** At least one `InteractionActions` entry is an
   `expression` or `activity` change (`increases`, `decreases`, or
   `affects`). Binding, secretion, phosphorylation, localization, abundance,
   and other action types do not count on their own.
4. **Measured gene (implementation check).** The row's own gene must be the
   one whose expression or activity was measured, read from the interaction
   text (`expression of GENE mRNA`, `activity of GENE protein`). CTD also
   files an interaction under every gene that takes part in it, so without
   this check AHR would be counted for "AHR affects the reaction [TCDD
   increases expression of CYP1A1]". This step removed 384 rows.

Each kept row records:

- **`response_type`.** `expression` for mRNA or an unstated form, `protein`
  for protein expression, `activity` for activity. Eight interactions do not
  say mRNA or protein; they are recorded as `expression` with a note.
- **`context`.** One of:
  - `direct`: the chemical alone changes the gene.
  - `cotreatment`: the chemical is part of a mixture or co-exposure.
  - `modulated_reaction`: another agent alters the chemical's effect, or the
    chemical alters another agent's effect, for example "TCDD inhibits
    [estradiol increases BRCA1]".
  - `cotreatment|modulated_reaction`: both.
- **Missing fields.** CTD records no tissue, dose, duration, or assay, so
  `tissue`, `dose_duration`, and `assay` are `NA`. The note on each row gives
  CTD's interaction text and NCBI GeneID; the cited PubMed IDs hold the
  details.

**Choices for E6 to confirm.** These choices are not covered by the frozen
wording:

- the measured-gene check
- keeping `cotreatment` and `modulated_reaction` rows, labeled rather than
  dropped
- reading form-less interactions as `expression`

### Counts

Rows remaining after each filter, in order:

| Step | Rows |
|---|---|
| All CTD interactions | 3,190,383 |
| 1. Dioxin-like chemical IDs | 65,053 |
| 2. *Homo sapiens* | 9,331 |
| 3. Expression or activity change | 8,342 |
| 4. Row's gene is the measured gene | 7,958 |

Rows kept per chemical:

| Chemical | After step 1 | After step 2 | After step 4 |
|---|---|---|---|
| 2,3,7,8-TCDD | 55,910 | 8,368 | 7,215 |
| PCB 126 | 3,997 | 619 | 515 |
| PCB 77 | 727 | 105 | 73 |
| PCB 118 | 1,497 | 86 | 54 |
| PCB 156 | 56 | 25 | 24 |
| 2,3,4,7,8-PeCDF | 145 | 42 | 21 |
| 2,3,7,8-TCDF | 2,463 | 36 | 12 |
| 1,2,3,4,7,8-HxCDF | 27 | 8 | 8 |
| 1,2,3,4,6,7,8-HpCDF | 23 | 8 | 8 |
| 1,2,3,7,8-PeCDD | 30 | 7 | 6 |
| PCB 105 | 46 | 9 | 6 |
| PCB 169 | 72 | 6 | 5 |
| 1,2,3,6,7,8-HxCDD | 18 | 4 | 4 |
| 1,2,3,4,6,7,8-HpCDD | 15 | 4 | 4 |
| PCB 114 | 4 | 4 | 3 |
| 1,2,3,4,7,8-HxCDD, 1,2,3,7,8-PeCDF, PCB 167, PCB 189 | 10, 5, 5, 3 | 0 | 0 |
| OCDD, OCDF | 0 | 0 | 0 |

As expected, TCDD dominates: 91% of kept rows.

Kept rows by context:

| Context | Rows |
|---|---|
| `direct` | 5,548 |
| `modulated_reaction` | 1,382 |
| `cotreatment` | 964 |
| `cotreatment\|modulated_reaction` | 64 |

Kept rows by response type:

| Response type | Rows |
|---|---|
| `expression` | 7,020 |
| `protein` | 508 |
| `activity` | 430 |

### Gene mapping

`map` results:

- **7,943 rows mapped.**
- **2 ambiguous.** `COX2` (CTD GeneID 4513, which is MT-CO2, not PTGS2) and
  `DOC1`.
- **13 unresolved.** `HLA-DRB4` and `KIR2DL5A` are on alternate haplotypes
  only. `ERV3-2`, `HGC6.3`, `LINC00312`, `PGBD3`, `RNA28SN4`, `RUNX1-IT1`,
  `SOS`, `TRK-TTT3-5`, and `WNT5` have no match.

These 15 rows need E7/G7-style review, using the CTD GeneID in each note.

**Independent check.** Each mapped row was also checked through a separate
route: CTD's NCBI GeneID, then HGNC `entrez_id`, then Ensembl ID. Result:
7,931 agree, 0 disagree, 12 cannot be checked (no HGNC Entrez–Ensembl
link).

### Evidence for the provisional publication rule

These counts cover the 7,943 mapped rows from 628 distinct publications.
Each cell is the number of genes with at least *n* distinct publications for
some gene–chemical pair:

| Rows used | n ≥ 1 | n ≥ 2 | n ≥ 3 | n ≥ 5 |
|---|---|---|---|---|
| All contexts | 4,393 | 1,693 | 711 | 146 |
| `direct` only | 4,175 | 1,628 | 667 | 133 |

**Large screening papers drive the 2-publication count.** 17 publications
contribute at least 100 rows each. They supply 6,735 of the 11,230
publication–row links; the largest, PMID 20106945, has 1,361 rows. Of the
1,753 gene–chemical pairs with at least 2 publications, 1,119 are supported
only by these large screens. Only 248 genes have at least 2 supporting
publications with fewer than 100 rows each.

So at this size, "2 publications" mostly means "appears in two genome-wide
screens", not "independently replicated". The team should weigh this when the
provisional rule is confirmed or revised. The draft genetic-risk set from the
November 12 test run had 102 genes.

## Literature-context list

Expression-only genes from the group research notes (for example SULT1E1,
NNMT, BCL2L11 reported as "BIM", CXCL10, IL15, and the HLA class II family) are
kept in a separate, labeled list for interpretation only. They are never
counted as risk or environment evidence, and they do not enter the overlap
calculation.

## Interpretation boundary (E11)

*To be written.* It must state that environmental evidence only prioritizes
genes, and that without patient-level exposure data no exposure causation can
be claimed.
