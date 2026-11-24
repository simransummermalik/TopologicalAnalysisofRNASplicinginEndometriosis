# Interpretation plan: how results will be reported

**Status:** approved without changes by Vania Zeledon on 2026-11-24 (drafted
2026-11-23), before the eight-sample results were produced. Recorded in
[`decisions.md`](decisions.md). Nothing in it has been
computed on real data. The only case/control output seen so far is the
exploratory one-case-versus-one-control comparison of 2026-11-20 (per-sample
graphs, cycle fraction only; see the 2026-11-23 correction in
[`decisions.md`](decisions.md)). The junction-usage baseline below has never
been computed on real samples.

**Amendment (section 6):** drafted 2026-11-24, before any eight-sample result
existed, after checking this plan against the project proposal
([`hodgetheorynew.pdf`](../assets/hodgetheorynew.pdf)). It adds three
exploratory analyses and two claim limits; it does not change the
confirmatory test or the four outcomes. Approved without changes by Vania
Zeledon on 2026-11-24, before any eight-sample result existed.

`scripts/interpret_results.py` applies these rules mechanically and writes the
outcome to `outcome.json`. The written results must use the matching wording
below.

## 1. What can be called a finding

Only the **primary candidate-set test** (`set_test.tsv` of the primary
scenario, 2026-11-17 inference rules) can support the word "significant", and
only at p < 0.05. Everything else is exploratory and is reported with its raw
p-value and no significance claim:

- the risk, environment and overlap set tests;
- gene-level p-values and BH-adjusted values;
- the six sensitivity scenarios;
- the baseline comparison in section 2.

With 4 cases and 4 controls the smallest possible p-value is 2/70 ≈ 0.029.
Among hundreds of genes tested at that resolution, some raw p-values below
0.05 are expected by chance alone, so a gene's raw p-value is never described
as evidence on its own.

## 2. Does the Hodge summary add anything? The baseline

Task S10 asks whether the Hodge decomposition adds information over a
simpler summary. The comparison must be fair: same genes, same samples, same
shared edges, same read counts, same test. Only the per-sample summary
changes.

- **Hodge summary (primary):** cycle fraction, the share of a gene's
  edge-flow energy in the cycle space. It depends on how the edges connect.
- **Baseline:** junction-usage entropy, the normalized Shannon entropy of the
  read shares over the gene's shared edges:
  `H = −Σ pₑ log pₑ / log E`, where `pₑ` is edge e's share of the gene's
  reads in that sample and `E` the number of shared edges. 0 means all reads
  on one edge; 1 means reads spread evenly. It uses no topology: rearranging
  how the edges connect would not change it.

The existing count and topology baselines (node count, edge count, cycle rank,
total reads) are kept as descriptive output. With one shared graph per gene,
the first three are identical in every sample, and total reads mostly
measures depth, so none of them can answer this question.

Both summaries go through the identical exact set test (mean over genes of
the absolute case-minus-control difference, 70 relabelings applied to all
genes at once). The script checks that its recomputation of the Hodge test
reproduces the Rust primary result exactly, and stops if not. It also reports,
descriptively:

- the Spearman correlation across tested genes between the absolute
  cycle-fraction effect and the absolute entropy effect (a high value means
  the two summaries rank genes alike);
- the Spearman correlation between the absolute cycle-fraction effect and
  cycle rank (a high value would suggest effects mainly track graph size).

## 3. Pre-declared outcomes and the wording to use

| `outcome` | When | Say | Do not say |
|---|---|---|---|
| `null` | primary p ≥ 0.05 | "No difference in cycle fraction between endometriosis and control samples was detected across the candidate genes." Report the p-value, the number of genes tested, and that with 4 versus 4 only large, consistent differences are detectable. | "No difference exists"; "splicing topology is unaffected"; any gene-level "hit" presented as a finding. |
| `difference_also_in_baseline` | primary p < 0.05 and baseline p < 0.05 | "Cycle fraction differed between groups across the candidate genes, and so did junction-usage entropy; the Hodge summary is not shown to add information beyond junction usage." | "Hodge analysis revealed…" as if it were unique to the method. |
| `difference_not_in_baseline` | primary p < 0.05, baseline p ≥ 0.05 | "Cycle fraction differed between groups across the candidate genes, while the junction-usage entropy baseline did not reach p < 0.05." State that the comparison is exploratory. | "Hodge decomposition outperforms simpler methods" (one baseline, one dataset, n = 4 per group). |
| `not_tested` | the primary test could not run | "The primary test could not be run" and why (for example a sample excluded, or no gene testable). | Any case/control conclusion. |

A null result is a valid result and is reported with the same prominence as
a positive one (task S10.6). The method's complexity is never a reason to
claim novelty.

**Sensitivity.** If the primary outcome changes in any of the six scenarios
(for example p crosses 0.05 or the signed effect changes sign), the result is
described as depending on the threshold choice, naming the scenarios.

**Unusual samples.** Before interpreting a primary p < 0.05, check
`set_scores.tsv` and `baseline_comparisons.tsv`. If one sample drives the set
(for example, it alone has the extreme score), say so.

## 4. Which genes are discussed

To avoid picking genes after seeing which look interesting:

- The genes discussed by name are the **five tested genes with the largest
  absolute cycle-fraction effect**, ties broken by smaller raw p-value, then
  gene ID. They are described as illustrative, not as findings.
- Any other gene discussed must be named here before results are produced.
  None are named.
- Each discussed gene is shown with its row in `gene_interpretation.tsv`
  (category, evidence sources and IDs, eligibility, effect, raw p) and its
  per-sample edge counts in `gene_junctions.tsv`, with a graph figure for one
  case and one control.

## 5. Claims and their limits

Every claim in the results, slides or website must carry its matching limit.

| If the text says… | It must also say… |
|---|---|
| a gene's graph has cycles or a high cycle fraction | Cycles are alternative paths through annotated exons; they do not show circular RNA. Back-splice junctions were not detected by this pipeline. |
| cycle fraction differs between groups | This is an association in 4 versus 4 samples of eutopic endometrium in the proliferative phase from one dataset (GSE315857); it does not show the difference causes or results from endometriosis. |
| a gene is an environmental-response gene | Its membership comes from curated dioxin-like compound responses, mostly in non-reproductive tissue (CTD); no exposure was measured in these subjects. |
| a gene is a risk gene | Membership comes from GWAS association reported for the gene; the associated variant may act through another gene. |
| overlap genes behave differently | The overlap set has 7 genes; it is not automatically the most important group. |
| any count-based result | Only uniquely mapped reads at annotated junctions are counted; ambiguous and unannotated junctions are excluded, and sequencing depth differs between samples (49.5–89.7 million pairs; corrected 2026-11-24 from 80.0). |
| any category comparison | Category tests are exploratory, and the categories differ greatly in size (96 risk-only, 558 environment-only, 7 overlap). |
| the project identified genes with large changes in splicing topology (proposal success criterion 7) | Genes were ranked by effect size; with 4 versus 4 samples no single gene can be shown to differ, so the ranked genes are illustrations, not individual findings. (Amendment A4.) |
| the proposal's back-splice example or synthetic back-splice test | The back-splice case was shown on synthetic graphs only; in the real data every cycle comes from alternative paths between annotated exons (for example exon skipping), because back-splice junctions were not detected. (Amendment A4.) |

## 6. Amendment of 2026-11-24: coverage of the project proposal

Checked against the project proposal, sections 1–5 above left three of its
commitments without a pre-declared analysis. All three are **exploratory**:
they report raw p-values or correlations, never "significant", and they do
not change the primary test or which outcome in section 3 applies.

### A1. The overlap part of the hypothesis

The proposal's hypothesis (section 24) says differences "may be especially
pronounced among genes associated with both genetic susceptibility and
environmental response". The overlap set test in section 1 asks whether
overlap genes differ between groups, not whether they differ *more* than
other candidates. The added comparison (`overlap_comparison.tsv`):

- **Statistic:** mean absolute cycle-fraction effect of tested overlap genes
  minus that of the other tested candidates (risk-only and
  environment-only together).
- **Test:** one-sided, because the hypothesis names a direction. Category
  labels are permuted across tested genes: every group of tested genes the
  size of the overlap group is a possible "overlap" group. The test is exact
  when there are at most 100,000 such groups, otherwise 100,000 random groups
  drawn with the fixed seed 20261124, with p = (extreme + 1) / (draws + 1).
- **Minimum:** at least 3 tested overlap genes, otherwise `not_assessable`.

| `verdict` | When | Say |
|---|---|---|
| `not_assessable` | fewer than 3 overlap genes tested | "Too few overlap genes could be tested to assess whether their splicing differs more than other candidates'." |
| `consistent_with_hypothesis` | overlap mean above both the risk-only and environment-only means, and p < 0.05 | "In an exploratory comparison, overlap genes had larger cycle-fraction effects than the other candidates." Add that the overlap group is small and that genes are not independent. |
| `not_supported` | otherwise | "Overlap genes did not show larger cycle-fraction effects than the other candidates." |

Do not say that overlap genes "are more affected by endometriosis": the
comparison is between effect sizes, and if the primary test is null those
effects are not themselves shown to be differences.

### A2. Which junctions carry a gene's change

The proposal's program "identifies the junctions contributing most strongly"
(section 19, step 10). For each discussed gene (section 4),
`junction_contributions.tsv` lists every shared edge with its case and
control mean read count, normalized signal, signed cycle component and
**cycle-energy share** (the edge's squared cycle component divided by the
sample's signal energy). A sample's shares add up to its cycle fraction, so
the case-minus-control share differences add up exactly to the gene's
effect; the script checks this and stops if not. Edges are ranked within the
gene by the absolute share difference.

Say "the change in GENE's cycle fraction is concentrated on the junction(s)
…", not that a junction is differently used: no junction-level test is run.
A shared edge with no reads in a sample can still carry a cycle component
(see [shared graphs](shared_graphs.md#interpretation-changes-to-explain)),
so read the counts beside the shares.

### A3. Expression level

The hypothesis says differences occur "not only in expression level". A gene's
cycle-fraction change could simply follow a change in how strongly it is
expressed. For every gene in the method comparison, the script computes
**gene junction expression**: log2 of the gene's reads on its shared graph per
million assigned junction reads in that sample (the whole input junction
table, all genes). It reports:

- the case-minus-control difference per gene (`gene_method_comparison.tsv`);
- the same exact set test as the primary test (`method_comparison.tsv`,
  labeled descriptive);
- the Spearman correlation across tested genes between the absolute
  cycle-fraction effect and the absolute expression effect.

Required wording, recorded in `outcome.json` under `expression_check`:

- if the primary outcome is a difference and the expression set test has
  p < 0.05: "Gene junction expression also differed across the candidate
  genes; the cycle-fraction difference is not shown to be separate from an
  expression difference."
- if the Spearman correlation is 0.5 or higher, whatever the outcome:
  "Genes with larger cycle-fraction effects also tended to have larger
  expression effects."

This is a junction-read measure, not a full expression analysis; say so if it
is reported.

### A4. Two claim limits

Two rows were added to the table in section 5: one for the proposal's
success criterion 7 (genes are ranked, not individually shown to differ) and
one for back-splicing (shown on synthetic graphs only).

### Amendment sign-off

- **Approved by:** sign-off by Vania Zeledon approves this amendment:
  - [x] Vania Zeledon (2026-11-24)
- **Date:** 2026-11-24
- **Changes from this draft, if any:** none

## Sign-off

- **Approved by:** sign-off by Vania Zeledon approves the plan (sections 1–5):
  - [x] Vania Zeledon (2026-11-24)
- **Date:** 2026-11-24
- **Changes from this proposal, if any:** none
