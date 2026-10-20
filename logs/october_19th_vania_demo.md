# October 19, 2026 — Vania (demo package)

**Status:** Pushed

**Commit:** `a1b13cc` — `Add examples/demo: a runnable demo package for the midterm presentation`

## Contribution

Vania built a small runnable demo package for the team's check-up code demo,
separate from and following the same day's **R10** work
([logs/october_19th_vania.md](october_19th_vania.md)).

- Ran the existing reviewed synthetic fixtures (`tests/fixtures/linear.tsv`,
  `cycle.tsv`, `exon_skip.tsv`) and the three deliberately invalid fixtures
  (`zero_support.tsv`, `duplicate_edge.tsv`, `self_loop.tsv`) to confirm their
  outputs match `tests/fixtures/expected.md` live, not just under `cargo test`.
- Ran the new `splice-girl gene` command (R10) against a hand-built
  case/control fixture to confirm the CLI wiring works end to end outside of
  its unit tests.
- Added `examples/demo/`:
  - `junctions.tsv` and `manifest.tsv` — four included samples (two case, two
    control) plus one sample of each exclusion type (`control_03` has no
    junction data for the gene; `orphan_01` has no manifest entry), so a
    single run demonstrates both exclusion reasons `compare_gene` reports.
  - `README.md` — the exact commands, expected output, and what each sample
    is there to prove, organized as three parts: core math, validation
    failure behavior, and the case/control comparison.
- The case/control counts in `junctions.tsv` deliberately reuse the same
  edge weights as `cycle.tsv` and `exon_skip.tsv` (not new values) so the
  resulting `cycle_fraction`s are recognizable during narration. This was a
  conscious tradeoff, discussed explicitly: the demo file is not exercised
  by `cargo test` and duplicates fixture data rather than referencing it, so
  the two copies could drift if the underlying fixtures are ever changed.
- Confirmed both the synthetic `analyze` path and the new `gene` path run in
  well under a second (`time` on the compiled binary), supporting the
  roadmap's "finishes in seconds" live-demo requirement.

## Current scope

This is an illustrative demo package, not the real-data demo described in
`ROADMAP.md` Phase 11 (which uses an actual candidate gene and the real
GSE179640/E-MTAB-15117 samples) — that remains a separate, later milestone.
`examples/demo/` is not referenced by any automated test.
