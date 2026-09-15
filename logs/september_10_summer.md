# September 10, 2026 — Summer

**Status:** Prepared for push; not yet committed when this log was created

## Contribution

Summer prepared the first basic version of **Splice Girl**, the project's Rust
command-line software. This is the synthetic mathematical foundation, not the
complete research tool.

### Rust project setup

- Added `Cargo.toml` and `Cargo.lock` for the `splice-girl` Rust package.
- Added `.gitignore` so Rust's generated `target/` directory is not committed.
- Added the command entry point in `src/main.rs`.
- Added the detailed `assets/splice_girl_banner.txt` startup figure and wired
  it into the CLI with terminal-safe color output.
- Added the `splice-girl --start` interactive menu with the currently available
  junction-TSV analysis action and a quit option.

### Basic input and graph processing

- Added TSV parsing and validation in `src/parser.rs`.
- Added grouping by sample and gene in `src/lib.rs`.
- Added deterministic splice-graph construction in `src/graph.rs`.
- Added incidence-matrix construction and a small dense linear solver in
  `src/matrix.rs`.

The current input columns are:

```text
sample_id	gene_id	from	to	count
```

### Graph-level Hodge calculation

- Added within-sample, within-gene junction normalization.
- Added the gradient projection and cycle-space residual in `src/hodge.rs`.
- Added cycle-energy and cycle-fraction calculations.
- Added checks for reconstruction, cycle-space membership, and orthogonality.
- Added per-edge signal, gradient-component, and cycle-component output.

### Mathematical specification and synthetic validation

- Added `docs/math_spec.md` with the orientation, normalization, equations,
  tolerance, expected examples, and limits of the basic dense solver.
- Added `tests/fixtures/linear.tsv`.
- Added `tests/fixtures/cycle.tsv`.
- Added parser, incidence-orientation, linear-graph, and cycle-graph tests.

Validation completed for this version:

- The linear graph reports a cycle fraction of `0.0`.
- The equal-weight directed cycle reports a cycle fraction of `1.0`.
- All four Rust tests pass.
- Rust formatting and strict lint checks pass.

### Documentation updates

- Updated `README.md` with the Splice Girl name, 2026 project year, basic usage,
  current scope, and test commands.
- Updated `ROADMAP.md` to mark the math and Rust phases in progress and to use
  the settled `splice-girl` command name.
- Updated `roadmap.tex` and rebuilt `roadmap-preview.pdf` with the Splice Girl
  name. The LaTeX preview compiled without formatting warnings.
- Added the dated contribution logs in `logs/`.

## Current scope

This contribution does not include real RNA-seq integration, gene research,
case/control statistics, gene ranking, or the final subcommand interface. Those
remain later project stages.
