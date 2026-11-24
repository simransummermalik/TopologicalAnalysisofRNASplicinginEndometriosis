# Results draft

Updated November 24, 2026. **Part 1 is a skeleton for the eight-sample
results**, which do not exist yet. Every number is a blank marked
**[…]** with the output file it comes from. Where the result decides the
wording, every pre-declared alternative is written out; keep the one that
applies and delete the rest. The wording is copied from the
[interpretation plan](interpretation_plan.md) (approved 2026-11-24); if the
two ever differ, the plan wins. **Part 2** reports the completed technical
pilot and stays as written.

File paths below are relative to `analysis/cohort_primary/results/`;
`interpretation/` is the output folder of `scripts/interpret_results.py`.

---

# Part 1. Eight-sample results (to fill in)

## 1. Samples and alignment

Sources: each run's `alignment_summary.txt` and `run_manifest.json`
(strand check); raw read pairs from the MultiQC tables in
`data/galaxy_data/october_20/`.

| Sample | Group | Raw read pairs (M) | Overall alignment | Read-1 antisense | Strand check |
|---|---|---|---|---|---|
| SRR36710383 | case | 80.0 | […] | […] | […] |
| SRR36710384 | case | 75.4 | […] | […] | […] |
| SRR36710385 | case | 89.7 | […] | […] | […] |
| SRR36710386 | case | 54.7 | […] | […] | […] |
| SRR36710387 | control | 49.5 | […] | […] | […] |
| SRR36710388 | control | 58.2 | […] | […] | […] |
| SRR36710389 | control | 72.5 | […] | […] | […] |
| SRR36710390 | control | 65.1 | […] | […] | […] |

Keep one:

- "All eight samples passed the strand check (read 1 antisense in […]–[…]%
  of eligible spliced alignments) and entered the analysis."
- "[Sample] was held by the strand check ([…]% read-1 agreement, below the
  95% rule) and excluded. The primary test requires all four cases and four
  controls, so it could not be run." → outcome `not_tested` (section 3).

## 2. Genes that could be tested

Source: `interpretation/category_summary.tsv`.

| | Candidates | With shared edges | Eligible (cycle, all 8 covered) | Tested |
|---|---|---|---|---|
| Risk-only | 96 | […] | […] | […] |
| Environment-only | 558 | […] | […] | […] |
| Overlap | 7 | […] | […] | […] |
| All | 661 | […] | […] | […] |

"Of the 661 candidate genes, […] had a shared graph with at least one cycle
that was covered in all eight samples and were tested. The others lacked
shared edges ([…]), had no cycle ([…]) or lacked coverage in at least one
sample ([…]); a gene without a cycle has a cycle fraction of zero by
definition." Counts from `analysis/primary/gene_eligibility.tsv`.

## 3. Primary result: candidate-set test

Sources: `interpretation/set_tests.tsv` (row `confirmatory`) and
`interpretation/outcome.json`.

"Across the […] tested candidate genes, the mean absolute case-minus-control
difference in cycle fraction was […] (mean signed difference […]); exact
permutation p = […] ([…] of 70 relabelings at least as extreme)."

Then the sentence for the outcome in `outcome.json` (keep one):

- **`null`:** "No difference in cycle fraction between endometriosis and
  control samples was detected across the candidate genes." Add the p-value,
  the number of genes tested, and: "With four samples per group, only large,
  consistent differences are detectable."
- **`difference_also_in_baseline`:** "Cycle fraction differed between groups
  across the candidate genes, and so did junction-usage entropy; the Hodge
  summary is not shown to add information beyond junction usage."
- **`difference_not_in_baseline`:** "Cycle fraction differed between groups
  across the candidate genes, while the junction-usage entropy baseline did
  not reach p < 0.05." Add: "This comparison is exploratory."
- **`not_tested`:** "The primary test could not be run because […]." No
  case/control conclusion follows.

Do not write: "no difference exists", "splicing topology is unaffected",
"Hodge analysis revealed …" as if unique to the method, or "Hodge
decomposition outperforms simpler methods".

If the outcome is a difference, check `analysis/primary/statistics/set_scores.tsv`
and `baseline_comparisons.tsv` first; if one sample alone drives the result,
say so: "[…]".

## 4. Does the Hodge summary add anything?

Source: `interpretation/method_comparison.tsv`, `outcome.json`.

| Summary (same genes, samples and test) | Genes | Set statistic | p (exploratory except primary) |
|---|---|---|---|
| Cycle fraction (primary) | […] | […] | […] |
| Junction-usage entropy (baseline) | […] | […] | […] |
| Gene junction expression, log2 (descriptive) | […] | […] | […] |

"Across tested genes, the Spearman correlation between absolute
cycle-fraction and absolute entropy effects was […]; between absolute
cycle-fraction effect and cycle rank, […]; and between absolute
cycle-fraction and absolute expression effects, […]."

Add each sentence listed under `expression_check.required_statements` in
`outcome.json`, word for word. If it lists none, write: "Cycle-fraction
effects did not closely follow gene junction expression (Spearman […])."
Note that this is a junction-read measure, not a full expression analysis.

## 5. Exploratory category tests

Source: `interpretation/set_tests.tsv` (rows `exploratory`). Raw p-values
only; none is described as significant.

| Set | Genes tested | Mean absolute effect | Raw p |
|---|---|---|---|
| Risk | […] | […] | […] |
| Environment | […] | […] | […] |
| Overlap | […] | […] | […] |

"These category tests are exploratory, and the categories differ greatly in
size (96 risk-only, 558 environment-only, 7 overlap)."

## 6. Overlap genes and the hypothesis

Source: `interpretation/overlap_comparison.tsv`. Keep one:

- **`not_assessable`:** "Too few overlap genes could be tested ([…]) to
  assess whether their splicing differs more than other candidates'."
- **`consistent_with_hypothesis`:** "In an exploratory comparison, overlap
  genes had larger cycle-fraction effects than the other candidates (mean
  absolute effect […] versus […]; one-sided permutation p = […])." Add that
  the overlap group is small and that genes are not independent.
- **`not_supported`:** "Overlap genes did not show larger cycle-fraction
  effects than the other candidates (mean absolute effect […] versus […];
  p = […])."

Do not write that overlap genes "are more affected by endometriosis".

## 7. Sensitivity to the thresholds

Source: `analysis/sensitivity.tsv` and each scenario's `statistics/set_test.tsv`.

| Scenario | Changed setting | Genes tested | Primary p | Signed effect |
|---|---|---|---|---|
| Primary | — | […] | […] | […] |
| Edge reads 3 | edge ≥ 3 reads | […] | […] | […] |
| Edge reads 10 | edge ≥ 10 reads | […] | […] | […] |
| Samples per edge 2 | ≥ 2 of 8 samples | […] | […] | […] |
| Samples per edge 6 | ≥ 6 of 8 samples | […] | […] | […] |
| Coverage 20 | ≥ 20 reads per sample | […] | […] | […] |
| Coverage 100 | ≥ 100 reads per sample | […] | […] | […] |

Keep one:

- "The primary outcome was the same in all six sensitivity scenarios."
- "The result depended on the threshold choice: in [scenarios], p crossed
  0.05 / the signed effect changed sign."

## 8. Genes discussed

Source: `interpretation/discussed_genes.tsv`, `junction_contributions.tsv`,
`gene_junctions.tsv`. Only these five genes are named, chosen by rule (largest
absolute effect; ties by raw p, then gene ID). They are illustrations, not
findings.

| Rank | Gene | Category | Case mean | Control mean | Difference | Raw p |
|---|---|---|---|---|---|---|
| 1 | […] | […] | […] | […] | […] | […] |
| 2 | […] | […] | […] | […] | […] | […] |
| 3 | […] | […] | […] | […] | […] | […] |
| 4 | […] | […] | […] | […] | […] | […] |
| 5 | […] | […] | […] | […] | […] | […] |

For each gene, one or two sentences: "The change in [GENE]'s cycle fraction
is concentrated on the junction(s) [from → to] ([…] of the […] difference)."
Give the mean read counts on those junctions beside the shares; no
junction-level test was run. Show one case and one control graph figure for
rank 1.

"With four samples per group, no single gene can be shown to differ; these
genes are ranked by effect size and shown as examples."

## 9. Limits that go with these results

Each claim made above carries its limit from the plan's claims table:

- Cycles are alternative paths through annotated exons; they do not show
  circular RNA, and back-splice junctions were not detected.
- Any difference is an association in 4 versus 4 samples of proliferative
  eutopic endometrium from one dataset; it does not show cause.
- Environmental-response genes come from curated responses to dioxin-like
  compounds, mostly in non-reproductive tissue; no exposure was measured in
  these subjects.
- Risk genes come from reported GWAS associations; the variant may act
  through another gene.
- Only uniquely mapped reads at annotated junctions were counted, and depth
  differs between samples (49.5–89.7 million pairs).

---

# Part 2. Completed technical pilot

This section reports work done on November 17 with **earlier, exploratory
settings**: unstranded alignment and one graph node per exon interval. The
cohort analysis uses RF alignment and the shared-splice-site node rule, so
none of these numbers belong to Part 1.

## Software and evidence inputs

Splice Girl implements input validation, graph construction, within-gene
normalization, graph Hodge decomposition, descriptive comparisons, explicit
exact permutation tests, candidate-set tests, BH correction, baseline exports,
support/coverage sensitivity runs and offline reports. The candidate table
contains 661 versioned gene IDs: 96 risk-only, 558 environment-only and seven
overlap genes. Those labels come from evidence curation and do not describe
the sequenced donors' exposure histories. Independent evidence sign-offs
remain listed in [spot_checks.md](spot_checks.md).

## Human-reference technical pilot

The first one million synchronized cleaned pairs from control SRR36710387
were aligned against the verified GENCODE v50 human primary assembly. HISAT2
reported 97.23% overall alignment. The annotation converter exported 32,575
junction rows. A predeclared coverage/size rule selected CDH11; at minimum
support three its graph retained six nodes, four edges and 103 read-junction
observations. Independent BAM recounting matched the six base edges and 106
observations, and the Rust projection agreed with an independent SVD to
approximately 4.6e-11 per exported component. See the
[versioned pilot bundle](../data/pilots/SRR36710387_prefix1m_nov17/README.md).

The selected graph had two connected components and cycle rank zero.
Its cycle fraction was effectively zero at support thresholds 1, 3, 5 and 10.
This is an expected mathematical result for this graph and a successful
technical check. It does not demonstrate a nonzero cycle component in real
data or establish absence of alternative splicing or circular RNA.

## Assignment quality and execution scale

Of 587,089 assignment-stage observations, 152,907 were uniquely assigned,
432,599 were ambiguous, and 1,583 were unannotated. The detailed audit found
that 409,613 of the ambiguous observations identified a single gene but
multiple full exon pairs. This loss led to the shared-splice-site node rule
adopted on 2026-11-18, under which 98.5% of annotated candidate junctions map
uniquely ([review](data_methods.md#node-rule-review-2026-11-18)). The
[audit figures](../data/pilots/SRR36710387_prefix1m_nov17/audit_figures/README.md)
show these diagnostics without presenting them as biological differences.

Across all 243 candidate genes observed in this prefix, the Rust analysis and
table export processed 1,388 junction rows in a median of 11.2 ms over five
runs, with a maximum measured child-process peak of about 4.0 MiB. All groups
passed internal math checks on every run. The largest observed graph had 66
nodes. This measures analysis/export on one prefix, not alignment, report
rendering or full-cohort runtime. The
[benchmark record](../data/pilots/SRR36710387_prefix1m_nov17/candidate_benchmark/benchmark.json)
retains per-run measurements and hashes. The 418 unobserved candidates were
not assigned zero-valued results.

No case/control effect, disease-associated gene, significant gene set, or
advantage over simpler baselines is claimed from this one-control pilot.
