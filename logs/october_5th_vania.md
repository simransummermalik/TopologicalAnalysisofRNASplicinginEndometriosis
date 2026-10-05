# October 5, 2026 — Vania

**Status:** Pushed.

## Contribution

Vania completed Phase 4 math task **M11** from `delegations_and_tasks/tasks.md`:
the exon-skipping synthetic fixture, its independently derived expected
result, and the documentation explaining why the result holds.

### The fixture

- Added `tests/fixtures/exon_skip.tsv`: the exon-skipping graph `E1 -> E2`,
  `E2 -> E3`, `E1 -> E3` with equal read counts. As a directed graph it has no
  return path to `E1`, but as an undirected graph it is a triangle with one
  independent cycle, so a nonzero cycle-space component is expected despite
  the absence of a directed cycle.

### Independent validation

- Solved the normal equations by hand (B1, the restricted graph Laplacian,
  and the fitted potential with `E1` pinned to zero) for general equal edge
  weight `w`. The cycle component works out to `(w/3, -w/3, w/3)` regardless
  of `w`, giving an expected cycle fraction of exactly `1/9`.
- Independently re-checked that derivation with a second, separate
  calculation (plain linear-algebra solve, not reusing the by-hand algebra)
  before trusting the number.
- Added `exon_skipping_graph_has_nonzero_cycle_component_despite_no_directed_cycle`
  to `tests/synthetic.rs`, asserting the engine's output matches the
  independently derived `1/9` to within `1e-12`.

### Documentation

- Added the exon-skipping example to `docs/math_spec.md` alongside the
  existing linear and cycle examples, explaining why an acyclic directed
  graph can still carry nonzero cycle-space signal.

Validation completed for this contribution:

- `cargo test` passes all six tests, including the new exon-skipping case.
- `cargo fmt --check` passes.
- `cargo clippy --all-targets -- -D warnings` passes.

## Current scope

This contribution adds one of the required Phase 4 synthetic fixtures. It does not implement the remaining math tasks (M10 branching tree, M13 noise, M14 weight sweep) needed to close R9's synthetic end-to-end gate, and it does not touch real RNA-seq data, gene research, statistics, or ranking.
