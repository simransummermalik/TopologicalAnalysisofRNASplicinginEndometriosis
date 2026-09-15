# September 12, 2026 — Vania

**Status:** Pushed
**Commit:** `966b39e` — `Added comments to the files in src.`

## Contribution

Vania added documentation and inline comments across the `splice-girl` Rust
source, with no changes to logic or behavior:

- `src/graph.rs` — module doc explaining the file's role in the parse ->
  graph -> incidence-matrix -> Hodge pipeline, plus doc comments on `Edge`,
  `SpliceGraph`, and `from_junctions` (deterministic node/edge ordering,
  within-group normalization, the open duplicate-edge decision).
- `src/hodge.rs` — module doc summarizing the Hodge decomposition
  (gradient part vs. cycle-space part) and what `cycle_fraction` measures,
  plus doc comments on `Checks`, `Decomposition`, `decompose`,
  `fit_potential` (why the Laplacian is pinned to a root node per
  component), and `connected_components`.
- `src/matrix.rs` — module doc noting the dense, dependency-free linear
  algebra is an MVP choice, plus doc comments on `incidence_matrix`,
  `matrix_vector`, `transpose_vector`, and `solve` (Gaussian elimination
  with partial pivoting).
- `src/parser.rs` — module doc on the TSV input contract, plus doc comments
  on `Junction`, `parse_tsv`, `parse_tsv_text`, and the self-loop rejection
  rationale.
- `src/lib.rs` — module doc on the library's role as the pipeline entry
  point, plus doc comments on `GroupResult` and `analyze_junctions`
  (grouping is always within a single sample/gene, never pooled).
- `src/main.rs` — module doc describing the current single-purpose CLI and
  noting the `analyze`/`gene`/`compare`/`rank` subcommands from
  `ROADMAP.md` are not built yet.

## Current scope

This contribution is documentation only; it does not change the math,
parsing rules, or CLI behavior established in the September 10 prototype.
