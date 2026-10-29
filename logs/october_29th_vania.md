# October 29, 2026 — Vania

**Status:** Pushed

**Commit:** `c69393a` — `Add docs/statistics_spec.md + sketch S1/S4/S5`

## Contribution

Vania sketched the first three parts of **`docs/statistics_spec.md`**, a file
required by Phase 8 task **S1** that did not exist before this contribution:
**S1** (biological unit and permutation restrictions), **S4** (primary
effect), and **S5** (permutation design).

### Why only three of S1–S12

S3 and S6–S12 need real per-sample `cycle_fraction` values to validate
against, baseline, or export — every sample's `junction_path` in
`data/metadata/samples.tsv` is still `pending_D6`, so there is no real
per-sample output yet. S1, S4, and S5 were chosen specifically because they
depend only on the sample *design*, which is already fixed and real, not on
data values that don't exist yet. S2 (freezing eligibility/support rules) was
also deferred, since the coverage/support thresholds it needs are a data-team
decision, not yet made.

### S1 — grounded in the real manifest, not assumed

Read `data/metadata/samples.tsv` directly rather than guessing the study
design. Its `subject_id_source` column states
`internal_one_sample_per_subject` for all eight samples, settling the
question outright: one sample per subject, unpaired, 4 control vs. 4 case,
with tissue and menstrual phase held constant across both groups by design
(a controlled factor, not a blocking variable). Conclusion: any relabeling of
the 8 subjects into 4 case / 4 control is a valid, exchangeable permutation
— with the group-size consequences worked out under S5, per the explicit
agreement to avoid stating the same number in two places that could drift
out of sync.

### S4 — documented existing behavior rather than designing new behavior

The primary effect (case mean minus control mean `cycle_fraction`) is not a
new decision — it is already what the `compare` command's per-gene summary
computes (confirmed against the demo's `GENE_HIGH_CYCLE` example, which
reports exactly this difference). This section records that as the frozen
statistical definition rather than inventing an alternative.

### S5 — the permutation space, worked out exactly

- With 4 control and 4 case subjects, the full permutation space is
  `C(8,4) = 70` distinct relabelings — small enough to enumerate exhaustively
  rather than sample with Monte Carlo methods, and the spec states this as
  the chosen procedure.
- Because exhaustive enumeration always includes the observed arrangement,
  the reported p-value can never be exactly `0` without any separate
  correction — directly satisfying the failure mode S6 exists to catch
  ("the code does not produce a zero p-value merely because too few
  permutations were requested").
- No random seed is needed for this cohort (nothing is sampled), with an
  explicit note that this must be revisited if a future, larger cohort ever
  makes exhaustive enumeration infeasible.
- Recorded the hard floor on achievable significance this sample size
  imposes regardless of effect strength: `p >= 1/70 ≈ 0.0143` one-sided,
  `p >= 2/70 ≈ 0.0286` two-sided (extremes pair with their exact complement
  under label-swapping) — written down now so nobody later mistakes "the
  smallest p-value this design could ever report" for a problem with the
  method.
- Left **one-sided vs. two-sided as an explicitly open decision**, not
  something decided unilaterally: the research question has no committed
  direction, so this is a scientific judgment call for whoever owns the
  statistics decision, not a mechanical one.
- Wrote the full procedure as six numbered steps, satisfying S5's own
  "Done when" requirement directly.

## Current scope

This is a specification document only — no Rust code changed, and no
permutation testing is implemented yet (R12 remains blocked on this and on
real per-sample data). S2, S3, and S6–S12 remain unwritten. The one-sided
versus two-sided decision is explicitly unresolved pending team input.
