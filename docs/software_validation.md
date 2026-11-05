# Splice Girl software validation

**Validation date:** November 5, 2026
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
Rust 2024 edition used by this package. The October 28 suite of 56 tests was
compiled and run with both Rust 1.85.0 and Rust 1.97.1. The current 71-test
suite, including reports and permutation testing, also passes on both versions.
The CI matrix retains Rust 1.85.0 as well as stable.

## Isolated installation check

The installation record below is from October 28, before the report command
was added. The November 3 report tests execute the current Cargo-built binary.

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

The current suite contains **71 passing tests**:

- 41 library tests for parsing, manifests, graph math, comparisons, outputs,
  validation summaries, SVG generation, report escaping/missing values, and
  permutation testing
- 9 command-level error tests
- 9 command-level workflow tests, including report/table agreement, links,
  overwrite protection, missing groups, unsafe identifiers, and menu execution
- 12 independent synthetic-math tests

The following checks pass:

```sh
cargo fmt --all -- --check
cargo test
cargo clippy --all-targets -- -D warnings
```

## November 5 — Rust 1.85 CI correction

The November 3 [GitHub run](https://github.com/simransummermalik/TopologicalAnalysisofRNASplicinginEndometriosis/actions/runs/41695189578)
passed formatting and all 71 tests on Rust 1.85.0, then failed strict Clippy on
`clippy::format_collect` in the SHA-256 hex formatter in `src/output.rs`.
The failure was reproduced locally with that exact toolchain. The formatter
now writes the eight zero-padded words into one preallocated string. The known
SHA-256 test vectors still match, and no lint checks were suppressed.

After the change, formatting, all tests, strict Clippy, release builds, and the
three-run candidate benchmark pass locally on Rust 1.85.0 and Rust 1.97.1.
The GitHub check can verify the fix after the updated code is committed and
pushed; the old run remains a record of the pre-fix commit.

## Offline report workflow

The `report` command uses the same analysis and descriptive comparisons as
`compare`. The command-level acceptance test checks that all four TSV exports
match `compare` after removing their intentionally different run IDs. The demo
report has two genes, eight sample/gene figures, explicit exclusions, and
case-minus-control effects of `+0.5555555556` and `-0.5555555556`.

The report and standalone graphs were checked in headless Chrome. All local
links, manifest output paths, input checksums, inline plot SVGs, and graph SVGs
were also inspected. These are synthetic software acceptance checks; real
biological comparability remains a separate review.

A single release-build report run with 200 synthetic genes, eight samples,
and 4,800 junction rows completed in **3.428 seconds**, producing 1,600 SVGs
and 1,606 total output files. This checks the complete report-writing path,
including figures and filesystem work. It is a single measurement on tiny
three-junction graphs, rather than a general runtime guarantee. Temporary
benchmark inputs and results were removed by the benchmark's temporary-folder
context manager.

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

## Candidate-sized synthetic performance

The repository also includes `scripts/benchmark_candidate.py`. It generates a
temporary, deterministic table with **200 genes, 8 samples, and 4,800 junction
rows**, then runs the same all-gene `compare` command ten times. The table is
synthetic and is used only to check software scale; it is not a substitute for
a real junction table.

The command is:

```sh
cargo build --release
python3 scripts/benchmark_candidate.py
```

On the validation machine listed above, the ten-run output was:

| Measurement | Result |
|---|---:|
| Input | 200 genes × 8 samples × 3 junctions |
| Junction rows | 4,800 |
| Median runtime | 14.229 ms |
| Mean runtime | 62.101 ms |
| Minimum runtime | 10.773 ms |
| Maximum runtime | 504.079 ms |

The maximum is a startup or system-scheduling outlier; the median is the more
useful summary for this small benchmark. The benchmark is intentionally
rerunnable because runtime depends on the computer and compiler. A real-data
benchmark should be added after the first checked candidate junction table
exists.

## Limits of this validation

- The performance measurements use synthetic data. A real-data benchmark is
  required after a candidate-sized checked junction table exists.
- The independent beginner clean-checkout review in task R17 is still pending.
- Successful software validation does not validate RNA-seq preprocessing,
  biological comparability, statistical inference, or gene interpretation.
