# Splice Girl statistics specification

**Spec version:** 0.1 (2026-10-29) — covers S1 (biological unit), S4 (primary
effect), and S5 (permutation design) for the fixed eight-sample GSE315857
cohort recorded in [`data/metadata/samples.tsv`](../data/metadata/samples.tsv).

S2 (freezing eligible samples/genes and support rules), S3 (validating real
per-sample input), and S6–S12 (implementation validation, FDR, baselines,
sensitivity, and export) are deferred. They need real per-sample
`cycle_fraction` values to work against, and every sample's `junction_path`
is still `pending_D6` — there is no real per-sample output yet to validate,
baseline, or export. Writing those sections now would mean guessing at data
that doesn't exist.

## S1 — Biological unit and permutation restrictions

**Source:** [math specification](math_spec.md), [ROADMAP.md](../ROADMAP.md)
Phase 8, and [`data/metadata/samples.tsv`](../data/metadata/samples.tsv).

- **Biological unit = subject = sample.** The manifest's `subject_id_source`
  column states `internal_one_sample_per_subject` for all eight samples:
  each subject contributes exactly one sample, so there is no
  repeated-measures or within-subject structure to preserve.
- **Design: unpaired, two independent groups, 4 control vs. 4 case.** No
  subject appears in both groups, and no column pairs a control subject to
  a specific case subject (no matched-pair design).
- **Tissue (`eutopic_endometrium`) and menstrual phase (`proliferative`) are
  constant across both groups.** This is a controlled design choice made
  during sample selection, not a blocking variable — it removes a known
  confound rather than requiring stratified or restricted permutation.
- Because the design is unpaired with no blocking factor, **any relabeling of
  the 8 subjects into "4 case / 4 control" is a valid, exchangeable
  permutation under the null hypothesis of no group effect.** What that means
  in practice for exactly this sample size is worked out in S5 below.

**Done when:** We know what can and cannot be shuffled. Met — nothing
here may be paired or blocked; condition labels alone are shuffled.

## S4 — Primary effect

**Source:** [math specification](math_spec.md) and
[ROADMAP.md](../ROADMAP.md) Phase 8.

- **Case summary statistic:** mean `cycle_fraction` across eligible case
  samples, for one gene.
- **Control summary statistic:** mean `cycle_fraction` across eligible
  control samples, for the same gene.
- **Primary effect = case mean − control mean.** This is not a new design
  decision — it is already implemented and exported by the `compare`
  command's per-gene summary (see the demo: `GENE_HIGH_CYCLE` reports a
  case-minus-control difference of `+0.5555555556`). This section exists to
  record that behavior as the frozen statistical definition, not to invent a
  different one.
- **Missing data:** a sample `compare_gene` already excludes (no manifest
  entry, no junction data for this gene, or a failed decomposition check) is
  excluded from both the mean and the group's sample count. It is never
  imputed as zero, consistent with the project-wide missing-vs-zero rule.
- **Variation:** record each group's eligible sample count alongside the
  effect. The specific variation statistic (standard deviation, range, or
  otherwise) is deferred to S12's export design, since the right choice
  depends on real distribution shapes this project does not have yet.

**Done when:** the same effect can be calculated from the TSV. Met — it
already is, by `compare`'s existing output.

## S5 — Permutation design

**Source:** [math specification](math_spec.md) and
[ROADMAP.md](../ROADMAP.md) Phase 8.

- **Observed effect:** the S4 primary effect (case mean − control mean
  `cycle_fraction`) for one gene.
- **Labels that may be shuffled:** the `case`/`control` condition attached to
  each `subject_id`. Each sample's own `cycle_fraction` values are fixed;
  only the group assignment is permuted. This is valid because S1 found no
  pairing or blocking structure that a shuffle would need to respect.
- **The permutation space is small enough to enumerate exactly, not sample.**
  With 4 control and 4 case subjects, the number of distinct ways to choose
  which 4 of the 8 subjects are labeled "case" is

  ```
  C(8, 4) = 70
  ```

  That is the entire space of relabelings under the null — not an estimate.
  **The chosen procedure is exhaustive enumeration of all 70 permutations,
  not Monte Carlo sampling of a subset.** There is no reason to sample when
  the full space is this small to compute.
- **No random seed is needed for this cohort.** S5 ordinarily asks for a
  saved random seed because sampled permutation testing is non-deterministic;
  exhaustive enumeration is not random, so there is nothing to seed. If a
  future, larger cohort ever makes exhaustive enumeration infeasible
  (`C(n, k)` becomes too large to compute directly), this spec must be
  revised to add a sampling procedure and a saved seed at that time — it is
  not an oversight to omit one now.
- **p-value definition:** for a one-sided test, `p` is the fraction of all 70
  permutations whose effect is at least as extreme as the observed effect, in
  the hypothesized direction, counting the observed arrangement itself. Because
  the observed arrangement is always one of the 70 being counted, the
  numerator is always at least 1, so this procedure can never report `p = 0`
  — no separate correction is needed to avoid that failure mode (the kind
  sampled permutation testing needs a "+1" correction for; see S6).
- **One-sided vs. two-sided: open decision, not yet made.** The research
  question ("do endometriosis samples organize splicing differently") has no
  committed direction — cycle fraction could plausibly be higher or lower in
  cases. A two-sided test is the more defensible default given that, but this
  is a scientific judgment call, not a mechanical one, and should be
  confirmed with whoever owns the statistics decision before S6 is
  implemented. For reference: two-sided extremity is counted symmetrically,
  since every permutation's exact complement (swap which 4 subjects are
  "case") produces the exact negation of its effect — extremes always come
  in matched pairs under this design.
- **Consequence for interpretation — the floor on achievable significance:**
  regardless of how strong the real biological effect turns out to be, the
  smallest possible p-value this design can ever produce is:
  - **One-sided:** `1/70 ≈ 0.0143` (only if the real split is the single most
    extreme of all 70 possible relabelings).
  - **Two-sided:** `2/70 ≈ 0.0286` (extremes come in matched pairs; see
    above).

  This is a hard limit imposed by the sample size (S1's 4-vs-4 design), not a
  flaw in the method. It is recorded here so that, once real results exist,
  nobody mistakes "the smallest p-value this test could ever report" for
  "evidence the method doesn't work."

**Done when:** the permutation procedure can be written as numbered
instructions. Met:

1. Fix the observed case/control labeling from `data/metadata/samples.tsv`.
2. Compute the observed effect (S4) for the gene under test.
3. Enumerate all `C(8,4) = 70` ways to relabel the 8 subjects into 4
   case / 4 control.
4. For each relabeling, recompute the same effect using the samples'
   existing `cycle_fraction` values (unchanged) under the new labels.
5. Count how many of the 70 relabelings (including the observed one) are at
   least as extreme as the observed effect, in the direction(s) specified by
   the one-/two-sided decision above.
6. Report `p = count / 70`.
