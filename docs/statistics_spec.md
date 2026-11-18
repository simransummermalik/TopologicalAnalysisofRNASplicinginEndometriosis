# Splice Girl statistics specification

**Spec version:** 0.5 (2026-11-17) — records S1 (biological unit), S4 (primary
effect), S5/S6 (exact permutation design and validation), the configured
statistics CLI, and the [inference rules](#inference-rules-decided-2026-11-17)
(two-sided tests, 4 samples per group, one primary gene-set test) for the study
in [`data/metadata/samples.tsv`](../data/metadata/samples.tsv).

The exact permutation and S8 correction code is connected to the `stats`
command, exports, and offline report. Its biological testing family and final
workflow policy remain open. S2/S3 (eligible real inputs), scientific review of
S7 (whether inference is allowed), and S9–S12 (baselines, sensitivity, and final
study results) are still pending. Synthetic tests establish software behavior
while real-data QC and scientific decisions are completed.

Software update: `stats` now exports same-sample descriptive baselines and
optional, explicitly requested rankings, preserving original evidence columns
in `gene_family.tsv`. `scripts/sensitivity.py` reruns descriptive Hodge effects
across an explicit support/coverage grid, saving every attempted setting and
missing outcome. These implementations do not approve the real-study grid,
family, cutoff, or interpretation. See [running_the_tool.md](running_the_tool.md).

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
- **One-sided vs. two-sided: two-sided (decided 2026-11-17).** The research
  question ("do endometriosis samples organize splicing differently") has no
  committed direction — cycle fraction could plausibly be higher or lower in
  cases. See [inference rules](#inference-rules-decided-2026-11-17). The CLI
  still requires the choice explicitly. For reference: two-sided extremity is counted symmetrically,
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

## S8 — Multiple-testing correction core (in progress)

**Implementation:** [`src/fdr.rs`](../src/fdr.rs), using the
Benjamini–Hochberg method proposed in `ROADMAP.md`. The calculations were
checked against [R's `p.adjust` documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)
and local R `p.adjust(..., method = "BH")` fixture answers.

`fdr::benjamini_hochberg(family_id, tests)` accepts named gene rows with raw
p-values and returns the family label, method, tested/untested counts, and a
row per gene containing both the original and adjusted p-value. It preserves
input order and gene IDs. Untested rows use `None` for both p-values; the
correction denominator counts only supplied valid tests. Missing rows are
never treated as measured zeros. Empty IDs, duplicate genes, and nonfinite or
out-of-range p-values are rejected with an explanation.

The numerical operation sorts the tested p-values, multiplies each by the
number of tests divided by its rank, applies a reverse cumulative minimum,
and caps the answers at one. That final cumulative minimum matters: simply
multiplying each p-value by its rank factor can give incorrectly ordered
adjustments.

| Gene | Raw p-value | BH-adjusted p-value |
|---|---:|---:|
| Toy A | 0.010 | 0.020 |
| Toy B | 0.040 | 0.040 |
| Toy C | 0.030 | 0.040 |
| Toy D | 0.002 | 0.008 |

Eight tests cover these hand-calculated values, reverse-minimum behavior,
ties and row order, untested genes, invalid inputs, empty/singleton families,
boundary values, and an example consuming the exact permutation module's
output. The field is called `adjusted_p_value`: it is a BH adjustment rather
than a probability that an individual gene is a false discovery.

**Decided 2026-11-17:** gene-level BH covers all testable candidate genes and
is exploratory; untested genes stay out of the denominator. The confirmatory
test is a single gene-set test, so no correction across tests is needed. See
[inference rules](#inference-rules-decided-2026-11-17). **Still required:**
review the configured outputs on real data before the study run. The helper does not choose biological genes, a test
direction, or a significance cutoff. Current `compare` and `report` commands
continue to provide descriptive results.

## Inference rules (decided 2026-11-17)

These rules were fixed before any real splicing result existed. Changing them
needs a new dated entry in [`decisions.md`](decisions.md).

### Why the design is shaped by the 70 relabelings

With 4 cases and 4 controls there are 70 relabelings, so a two-sided p-value
can only be 2/70 (0.029), 4/70 (0.057), or larger. Only the single most extreme
split is below 0.05. Two consequences follow:

- **Gene-by-gene correction cannot succeed.** Under BH at 0.05 with `m`
  tested genes, a gene at p = 2/70 passes only if at least
  `0.0286 m / 0.05 ≈ 0.57 m` genes are also at 2/70. With hundreds of genes
  that will not happen, however real the biology.
- **Only one confirmatory test is possible.** Any correction over 2 or more
  tests (Bonferroni or Holm) needs a per-test p of at most 0.025, which is
  below 2/70.

### Rules

1. **Direction:** every test is two-sided (`--alternative two-sided`).
2. **Testable gene:** all 4 cases and all 4 controls must be eligible
   (`--min-per-group 4`). With 3 versus 4 there are 35 relabelings and the
   smallest two-sided p is 2/35 = 0.057, so such a gene could never pass.
3. **Primary (confirmatory) test:** one gene-set test on the full candidate
   set, [`candidate_genes.tsv`](../data/gene_sets/candidate_genes.tsv), at
   α = 0.05. This is the only result that may be called significant.
4. **Set statistic:** the mean over genes of the absolute case-minus-control
   cycle-fraction effect. For each of the 70 relabelings, the same relabeling
   is applied to every gene and the statistic is recomputed; `p` is the
   fraction at least as large as the observed value. Genes changing in
   opposite directions do not cancel. Applying one relabeling to all genes
   keeps the correlation between genes, so the test stays exact. A small
   `p` means the set's genes differ between groups more than under almost
   every relabeling; it does not give a direction. The signed mean effect is
   reported beside it for direction only.
5. **Genes in the set test:** genes that pass rule 2 and have every common
   sample. Others are counted as left out, never imputed.
6. **Secondary (exploratory) results**, reported with raw p-values and no
   significance claim:
   - set tests for the risk set
     ([`endometriosis_risk.tsv`](../data/gene_sets/endometriosis_risk.tsv),
     103 genes), the environmental set
     ([`environment_genes.tsv`](../data/gene_sets/environment_genes.tsv),
     565 genes), and the overlap set
     ([`overlap.tsv`](../data/gene_sets/overlap.tsv), 7 genes);
   - gene-level p-values with BH across all testable candidate genes
     (untested genes excluded from the denominator, as implemented).
7. **Reporting:** state the 2/70 floor, the number of testable genes and genes
   used, and that gene-level results are for ranking and follow-up, not
   discovery claims.

### Commands

```sh
# Primary
splice-girl stats junctions.tsv data/metadata/samples.tsv \
  --genes data/gene_sets/candidate_genes.tsv --family candidates-primary \
  --alternative two-sided --design independent --min-per-group 4 \
  --output runs/stats-candidates
# Secondary: repeat with endometriosis_risk.tsv, environment_genes.tsv and
# overlap.tsv, each with its own --family label and output folder.
```

The primary result is `set_test.tsv` from the candidate run.

### Limits to keep in mind

- Every gene in a set shares the same 8 samples. One unusual sample (depth,
  RNA quality, batch) can shift many genes at once and drive the set test.
  Check `baseline_comparisons.tsv` and the per-sample scores in
  `set_scores.tsv` before interpreting a small set p-value.
- The set statistic weights every gene equally, so many unchanged genes dilute
  a few strongly changed ones.
- Still open: junction support and coverage thresholds (S2/S3), the real-data
  exchangeability review (S7), and the final candidate table (pending spot
  checks and the G1 cold read).

## Configured statistics workflow (R11 and software parts of S7/S12)

```sh
cargo run -- stats examples/demo/all_genes.tsv examples/demo/manifest.tsv \
  --genes examples/demo/stats_genes.tsv --family demo-two-genes \
  --alternative two-sided --design independent --min-per-group 2 \
  --output target/demo-stats
```

The command requires a TSV with a `gene_id` column, a family label, an explicit
alternative, an independent-sample declaration, a minimum eligible sample
count per group (at least two), and a new output directory. It keeps every
requested gene, including absent or untestable ones, and records a reason for
every untested row. Genes with a successful exact test form the BH denominator;
this handling of untested genes must be accepted before the real study run.

Per-gene exclusions remain in `comparisons.tsv`. `statistics.tsv` records raw
and adjusted p-values, the effect, group sizes, extreme/total permutation
counts, tested family size, and reasons. The manifest records settings and
checksums for the junction, sample-manifest, and gene-family inputs. The
offline report displays statistics alongside the existing sample plots and
decomposition figures. If no gene is testable, diagnostic files are saved and
the command exits with a failure code.

The independent-design flag is a caller declaration, not proof that conditions
are exchangeable. Known duplicate subject IDs are rejected when that metadata
column exists. The code does not infer paired or blocked designs, remove batch
confounding, validate tissue matching, apply extra metadata eligibility/QC
fields, or approve biological coverage thresholds. Supply reviewed inputs.

The eight-subject CLI fixture has four cases with cycle fraction 1 and four
controls with fraction 0 for one gene, and a second no-difference gene.
Two-sided raw p-values are `2/70` and `1`; their BH-adjusted values are `4/70`
and `1`. Additional absent/undersampled genes remain untested and do not change
that denominator. Other tests verify both one-sided alternatives, required
configuration, repeated subjects, enumeration limits, and interactive execution.
