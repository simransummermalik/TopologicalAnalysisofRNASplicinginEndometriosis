# November 22–23, 2026 — Vania

**Status:** Pushed

**Commits:**

- `9c2cc19` (Nov 22) — Approved junction thresholds as proposed
- `1dd228d` — Correct threshold timing statement + add post-alignment
  strandedness check
- `6c7dabb` — Update the strandedness check to exclude failed samples
- `5f06284` — Add assemble_runs.py
- `e431a4b` — Primary cohort analysis configuration
- `c8a42d5` — Added interpret_results.py
- `47612bd` — Added junction-usage baseline
- `a06390d` — Proposal and interpretation sign-off update

## Contribution

I approved the junction thresholds, then built everything needed between
Summer's cluster runs and the final interpretation, so the real eight-sample
analysis can be run with a few checked commands. All of it is tested on
synthetic or mock data; no real case/control result was produced.

### Thresholds approved, and an honest correction

- Approved the threshold proposal without changes (Nov 22).
- The approval entry said the values were "chosen before any case/control
  result existed". Summer's November 20 log shows a one-case-versus-one-control
  comparison was viewed on November 20. I added a dated correction: the values
  were written on November 18, before that comparison, and approved
  unchanged; they are now frozen.
- Sign-off for the threshold proposal and the interpretation plan is by
  Summer and me (dated entry in `docs/decisions.md`).

### Strandedness check after alignment

- New `preprocessing/strand_orientation.py`: for primary, uniquely mapped,
  spliced reads, it compares each read's strand with the transcript strand
  HISAT2 records and calls the sample pass or hold for review (95% rule).
- The junction counter now records these counts for every new run, and the
  pipeline writes the call into each run's manifest.
- A held sample is left out of the combined cohort table, and the Rust
  analysis refuses one as a second safeguard. Older runs can be checked from
  their saved BAM.

### Combining separately run samples

- New `preprocessing/assemble_runs.py`. Summer runs each sample as its own
  cluster job, but the analysis needs one checked cohort record. The script
  combines finished runs and checks before writing anything:
  - every stage completed;
  - every sample used the same settings and the same reference files;
  - no count table changed after its run;
  - each sample has a strand check made for that sample.
- The Rust analysis accepts the new record. Tested with fake runs that pass
  the real `splice-girl check --config`.

### Primary analysis configuration

- New `analysis/cohort_primary/analysis.json` with the approved thresholds,
  the 661 candidate genes, the statistics rules and the provenance settings,
  plus a README tracing each value to its decision and giving the cluster
  commands.
- Tested end to end on a mock cohort with the real sample IDs, sample sheet
  and candidate list: it passes `check` and completes `run`, including the
  gene-set test.

### Interpretation (Priority 8)

- New `scripts/interpret_results.py`, run after the analysis:
  - one row per candidate gene with its category, evidence sources and IDs,
    eligibility and result;
  - per-sample edge counts for every tested gene;
  - a summary by category and the exploratory risk, environment and overlap
    set tests beside the confirmatory candidate test.
  - It stops if any gene's category disagrees with the evidence tables.
- Added a junction-usage entropy baseline, run through the identical
  permutation test, to answer whether the Hodge summary adds anything over a
  simpler summary. The script recomputes the primary test and stops if it
  does not match the Rust result.
- New `docs/interpretation_plan.md`, written before any real result:
  - wording for each pre-declared outcome, including a null result;
  - the rule for which five genes are discussed;
  - a table pairing each claim with its required limit.
  The script records which outcome applies.

## Current scope

The code path from cluster runs to interpreted results is complete and
tested. Waiting on: Summer's remaining five samples, Summer's threshold
sign-off, and both sign-offs on the interpretation plan.
