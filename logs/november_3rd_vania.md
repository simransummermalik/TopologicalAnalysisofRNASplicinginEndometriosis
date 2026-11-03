# November 3, 2026 — Vania

**Status:** Pushed

**Commit:** `e8e889d` — `Implemented permutation testing in src/permutation.rs`

## Contribution

Vania implemented the exact permutation test specified by **S5** in
`docs/statistics_spec.md`, as a new `src/permutation.rs` module wired into
`src/lib.rs`. This is the first Rust code in the project that performs any
inferential statistics; everything before it (R10's `compare_gene`, R11's
`compare_all_genes`) was explicitly descriptive-only.

### Scope decided before writing code

The module is library-only, not CLI-wired. Wiring it into the `compare`
command would require defaulting to one-sided or two-sided testing, and S5
explicitly left that as an open scientific decision rather than a mechanical
one — so this contribution builds and tests the function itself and leaves
that decision, and the CLI wiring that depends on it, for later.

### Design decisions, each made explicit rather than assumed

- **Reused `compare::IncludedSample`** as the input type instead of
  inventing a new one — the eligible per-sample `(sample_id, condition,
  cycle_fraction)` data already exists there.
- **`C(n, k)` is computed per call**, not hardcoded as the full-cohort `70`.
  A specific gene may have fewer eligible samples than all eight if some
  failed decomposition checks, and the result records
  `permutations_considered` so a reader can tell when fewer than 70 were
  used.
- **Added an explicit combinatorial-explosion guard** (`MAX_PERMUTATIONS =
  1_000_000`): `binomial_coefficient` fails loudly, naming
  `docs/statistics_spec.md`, rather than attempting to enumerate an
  infeasible space if the cohort ever grows. This operationalizes a caveat
  that was previously just prose in the spec.
- **`Sidedness` is a required enum parameter**, not a default, so the
  implementation can't silently resolve the one-sided/two-sided question
  S5 left open.
- **Reused `hodge::TOLERANCE`** for the "at least as extreme" floating-point
  comparison instead of inventing a new tolerance constant.
- **Hand-rolled the `C(n,k)` combination generator** (no external crate),
  consistent with the project's existing dependency-free convention in
  `src/matrix.rs`.

### Verification

- A fully hand-enumerable fixture (2 case, 2 control, `C(4,2) = 6`) where all
  six relabelings' effects were computed by hand in a doc comment and
  asserted exactly, same approach as M11's hand-derived `1/9`.
- A full eight-sample fixture with case values strictly greater than every
  control value, which **independently rederives** the spec's documented
  p-value floor (`1/70` one-sided, `2/70` two-sided) rather than asserting a
  number copied from the spec.
- A no-difference fixture (`p = 1.0`), an order-invariance test (shuffling
  sample order doesn't change the result), and an explicit test that an
  empty case or control group is rejected with a clear error.
- `cargo test` passes all 65 tests (40 lib including the 9 new, plus 9 + 4 +
  12 across the three integration test files).
- `cargo fmt --check` and `cargo clippy --all-targets -- -D warnings` both
  pass.

## Current scope

This implements S5's procedure only. It is not wired into the `compare` CLI
command, does not implement S8's FDR correction, and the one-sided versus
two-sided decision remains explicitly open pending team input. R11 cannot
move from "In progress" to fully "Done" until that decision is made and the
CLI is wired to use this module.
