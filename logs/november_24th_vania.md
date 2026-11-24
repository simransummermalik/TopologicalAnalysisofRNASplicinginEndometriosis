# November 24, 2026 — Vania

**Status:** Pushed

**Commits:**

- `7acd044` — Finalize sign-offs + interpretation plan fix
- `371d398` — Spot check E-1
- `2455b53` — Update methods draft
- `82f4435` — Results draft skeleton + corrections
- `9c6e22a` — Updates (priorities and task board)
- `dfcbb88` — Close dataset review + finish off priority 3
- `83629df` — More updates (stale approval wording)
- `2cb7eb0` — README.md update (team, roles, pseudocode)

## Contribution

I finalized the remaining study decisions, checked them against the project
proposal, prepared the methods and results write-ups, and got the GitHub page
ready for submission. No eight-sample case/control result was produced; the
last samples were still being aligned.

### Decisions finalized

- Sign-off for the threshold proposal and the interpretation plan is by me
  alone (dated entry in `docs/decisions.md`). I confirmed the thresholds as
  final and approved the interpretation plan, before any eight-sample result.
- Closed the dataset review: GSE315857 is accepted. The kit and read-length
  disagreement between the GEO record and the paper does not affect any
  setting, because orientation and read length were measured from the data.
  Checksums and QC are checked per run; one subject per run is a stated
  assumption.
- Reviewed and approved the CDH11 technical pilot (Priority 3). Its review
  points (mapping ambiguity, strandedness, thresholds) were already settled.

### Interpretation plan amendment

Checking the plans against the project proposal showed three of its
commitments had no pre-declared analysis. I added an amendment (plan
section 6), implemented it in `scripts/interpret_results.py` with tests, and
approved it before any result existed:

- whether overlap genes have larger effects than other candidates (the
  hypothesis's overlap claim), by a one-sided gene-label permutation;
- which junctions carry each discussed gene's change; the per-junction
  shares add up exactly to the gene's effect, and the script checks this;
- whether cycle-fraction effects simply follow gene expression.

The confirmatory test and the four pre-declared outcomes did not change.

### Evidence check

- Resolved spot check E-1 (MMP2, PCB 126). CTD's "decreased" row is
  zebrafish; the human row says "increased" with three human papers, as our
  table records. Marked pass, with the reviewer's original note kept.

### Write-ups and corrections

- Rewrote the methods draft to the frozen cohort settings; only values from
  the run remain to fill in.
- Added the eight-sample skeleton to the results draft: every pre-declared
  outcome sentence, with blanks naming their source files.
- Corrected the read-depth range to 49.5–89.7 million read pairs (it had been
  given as 49.5–80.0); no setting depends on it.
- Brought the priorities and task board up to date and fixed wording that
  still called settled items provisional.

### Submission

- Added the team (Splice Girls), members and roles, and algorithm pseudocode
  to the README for the November 24 deadline.

## Current scope

All code is finished and tested. Waiting on: Summer's remaining five sample
runs, then filling in the results; the spot checks and G1 cold read by a
teammate; a beginner test; and rehearsal.
