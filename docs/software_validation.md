# Splice Girl software validation

**Validation date:** October 28, 2026  
**Package version:** `0.1.0`

This note records installation, compatibility, workflow, test, and performance
checks for the current Rust tool. It is development evidence, not a biological
benchmark.

## Environment

- macOS 26.6.2 on ARM64
- Current compiler: `rustc 1.97.1`
- Declared and tested minimum: `rustc 1.85.0`
- Rust edition: 2024
- External Rust dependencies: none

Rust 1.85 is the minimum because it is the first stable release supporting the
Rust 2024 edition used by this package. The full test suite was compiled and
run with both Rust 1.85.0 and Rust 1.97.1.

## Isolated installation check

The Cargo files, `src/`, and `assets/` were copied into a new temporary source
directory. The following installation completed successfully without editing
source code:

```sh
cargo install --path .
```

The optimized installation took approximately 2.7 seconds on the machine
listed above. The installed executable successfully completed:

- `splice-girl --help`
- `splice-girl --version`
- `splice-girl --start`, followed by `q`
- `splice-girl validate examples/demo/junctions.tsv`
- `splice-girl compare examples/demo/all_genes.tsv examples/demo/manifest.tsv`
- all-gene export to a new result directory
- SVG generation for `case_02` and `GENE_DEMO`

The all-gene export contained all five expected files:

1. `sample_metrics.tsv`
2. `junction_components.tsv`
3. `comparisons.tsv`
4. `gene_summaries.tsv`
5. `run_manifest.json`

## Automated verification

The current suite contains **56 passing tests**:

- 31 library tests for parsing, manifests, graph math, comparisons, outputs,
  validation summaries, and SVG generation
- 9 command-level error tests
- 4 command-level workflow tests
- 12 independent synthetic-math tests

The following checks pass:

```sh
cargo fmt --all -- --check
cargo test
cargo clippy --all-targets -- -D warnings
```

## Small-demo performance

The release executable ran the all-gene comparison on the original one-gene
`examples/demo/junctions.tsv` and `examples/demo/manifest.tsv` 100 times.
Standard output was discarded during timing so terminal rendering did not
dominate the result.

| Measurement | Result |
|---|---:|
| Median runtime | 3.439 ms |
| Mean runtime | 5.675 ms |
| 95th percentile | 4.872 ms |
| Minimum runtime | 2.925 ms |
| Maximum runtime | 186.434 ms |
| Observed child-process peak memory | 1,835,008 bytes (about 1.75 MiB) |

The maximum reflects a single startup outlier. On this tiny input, performance
is not a blocker and there is no evidence that parallel or sparse processing
is currently needed.

## Limits of this validation

- The performance measurement uses synthetic demonstration data. A second
  benchmark is required after a candidate-sized real junction table exists.
- The independent beginner clean-checkout review in task R17 is still pending.
- Successful software validation does not validate RNA-seq preprocessing,
  biological comparability, statistical inference, or gene interpretation.
