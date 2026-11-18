# Splice Girl software validation

**Validation date:** November 17, 2026
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
compiled and run with both Rust 1.85.0 and Rust 1.97.1. The current 100-test
suite, including reports and the full configured statistics workflow, passes
on both versions.
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

On November 17 the current source was copied again into a fresh temporary
source folder and installed with `cargo install --path ... --root ... --locked
--offline`. The installed executable passed version/help/menu checks and a
complete `stats --rank-by adjusted-p` demo, producing its ranking table and
offline report. Installation and output folders were isolated and removed
after the check; the user's installed executable was not replaced.

## Automated verification

The current Rust suite contains **100 passing tests**:

- 59 library tests for parsing, manifests, graph math, comparisons, outputs,
  validation summaries, SVG generation, report escaping/missing values, and
  permutation testing, and BH correction (the four toy adjustment vectors
  were independently checked with R's `p.adjust`)
- 9 command-level error tests
- 9 command-level workflow tests, including report/table agreement, links,
  overwrite protection, missing groups, unsafe identifiers, and menu execution
- 11 statistics CLI acceptance tests covering independently known p-values,
  BH correction, one-sided alternatives, required options, untested rows,
  repeated subjects, enumeration limits, provenance, reports, menu option 7,
  baseline/sample agreement, topology flags, optional rankings, tie handling,
  and exact-ID evidence joins for tested and untested genes
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

The November 16 `stats` workflow adds `statistics.tsv` and a statistical
section to the same offline report. The eight-subject acceptance fixture
reproduces raw p-values `2/70` and `1`, adjusted values `4/70` and `1`, and
keeps absent/undersampled genes untested. The two-gene shipped demo produces
two raw and adjusted values of `1/3` from six exact assignments per gene.
Settings and all three input checksums are retained in the manifest.

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

- The original performance measurements use synthetic data. The real-prefix
  benchmark below covers observed candidates from one control; full-cohort
  performance remains unmeasured.
- The independent beginner clean-checkout review in task R17 is still pending.
- Successful software validation does not validate RNA-seq preprocessing,
  biological comparability, statistical inference, or gene interpretation.

## November 17 — connected software acceptance check

`python3 scripts/check_project.py --with-tools` runs the complete local check.
The default check needs only Cargo/Python; the extra smoke test requires
Cleaver and HISAT2 on `PATH`. The source test helpers use unique temporary
folders even when tests start simultaneously.

- Rust: 95 tests, formatting, strict Clippy and release build pass on stable
  Rust 1.97.1 and the declared minimum 1.85.0.
- Python: 22 preprocessing tests, 56 gene-evidence helper tests and 5 script
  tests pass. Total automated Rust/Python tests: **178**.
  This includes the repository's expanded evidence-helper suite as well as
  the paired-prefix and independent-recount regression tests.
- The actual Cleaver 1.0.1 / HISAT2 2.2.3 smoke test builds the tiny reference,
  cleans and aligns two synthetic paired samples, counts annotated junctions,
  sorts/indexes BAMs through Cleaver, merges samples, and generates the Rust
  report. An initial single-end smoke test passed too. This checks executable
  interoperability rather than only printed command plans.
- The support-grid test reruns the real Rust executable and verifies that
  disappearing genes retain blank effects rather than zero effects.
- A separate NumPy 2.4.6 SVD-based least-squares check agrees with Rust on
  11 sample/gene groups across seven valid mathematical fixtures, including
  disconnected, noisy, cyclic and exon-skipping shapes. The independent BAM
  recount also reproduces the external-tool fixture's junction count.
- The complete ranked demo report was inspected in headless Chrome. Input
  hashes, table run IDs, graph SVGs, preserved gene-family contents, and all
  report file links and section anchors were checked.
- The 200-gene, eight-sample, 4,800-row synthetic comparison benchmark completed
  three runs with a median of 21.817 ms on this machine. Runtime varies.

No full human-reference alignment or final cohort result is implied by these
checks. The tiny reference's intron was corrected to canonical GT/AG for the
splice-aware alignment fixture; its exonic sequence and expected mapping
remain synthetic.

## Real-data technical pilot — November 17

The [CDH11 pilot bundle](../data/pilots/SRR36710387_prefix1m_nov17/README.md)
records an actual alignment against the complete verified human reference,
using the first 1,000,000 cleaned read pairs from control `SRR36710387`.
HISAT2 reports 97.23% overall alignment. The standard converter exports
32,575 annotated junction rows across the prefix.

The predeclared coverage/size rule selected CDH11 (`ENSG00000140937.16`), the
only candidate meeting that pilot rule. Independent indexed-BAM recounting
matches all six base-table edges and 106 supporting junction observations.
At support 3, the retained graph has six nodes, four edges and 103 supporting
observations. Rust and NumPy's independent SVD projection agree edge by edge
to within approximately 4.6e-11. The graph has cycle rank zero; its zero cycle
fraction is expected. Supports 1, 3, 5 and 10 give the same qualitative result,
with a separate Rust run checked against NumPy at each threshold.

The report was inspected in headless Chrome. All recorded input/output hashes,
local report links, SVG XML and stage-completion manifests were checked.
Long genomic node names now use distinct plot IDs with a coordinate key;
the stored node IDs and mathematical calculations are unchanged.

This establishes one bounded real-input technical check, not a cohort
comparison, a nonzero real-data cycle demonstration, or approval of all
preprocessing choices. The pilot records a 73.7% ambiguous-observation exclusion
rate at the annotation-assignment stage, which requires review before expanding.
The recorded single selected-gene Rust analysis/export took about 0.57 seconds
including startup and file writing. The observed-candidate benchmark below
extends this check. Optional macOS kernel-level memory profiling was unavailable;
the actual reference and alignment runs completed and passed verification.

## Final local checks and observed-candidate benchmark — November 17

The current full check passes **100 Rust tests and 88 Python tests (188 total)**,
formatting, strict Clippy and the release build. Rust tests and strict Clippy
also pass with the minimum Rust 1.85.0. The Python count is 27 preprocessing,
56 evidence-helper and five script tests. New regressions distinguish
within-gene exon ambiguity from cross-gene ambiguity, verify mate-count and
condition checks in the storage preflight, and establish that BH significance
is not mathematically impossible at the balanced permutation p-value floor.
The unequal three-versus-four counterexample attains a two-sided p-value of
1/35; the study's four-per-group rule remains its predeclared completeness policy.

`scripts/benchmark_real.py` filtered the actual prefix's junction table to all
observed members of the fixed 661-gene family: **243 genes, 1,388 junction rows**.
The remaining 418 candidates were recorded as unobserved, not assigned zero.
The largest graph had 66 nodes; the maximum edge count was 52. Over five runs,
analysis plus TSV/manifest export took a median **11.194 ms** (range 10.355–16.516
ms, startup included). Maximum measured child-process peak RSS was **4,243,456
bytes (4.05 MiB)**, read through a fresh worker's `getrusage` rather than the
unavailable kernel profiler. Every one of the 243 groups passed the Rust math
checks on every run. These are internal checks across the candidate prefix;
only the selected CDH11 graph has the separate BAM/SVD verification above.

The [benchmark record](../data/pilots/SRR36710387_prefix1m_nov17/candidate_benchmark/benchmark.json)
stores source/binary/script hashes, scope, sizes and each run's measurements;
the exact filtered junction input is saved beside it. This benchmark excludes
alignment, HTML/SVG generation and disease inference. It does not predict
full-cohort runtime or establish biological validity of the other genes.

## Isolated source-package acceptance — November 17

The unpublished allowlisted source archive was extracted into a new temporary
folder. The archive and all 132 included file hashes matched its manifest;
no raw reads, references, alignment files or run folders entered the payload.
The complete default check passed from that extracted source with Cargo
offline: all 188 tests, formatting, strict Clippy, release build and the
synthetic benchmark. An isolated `cargo install --path . --locked --offline`
then passed help, version, menu/quit, synthetic validation/report generation
and live recomputation of the four-edge CDH11 input. The normal installed
executable in the user's account was not replaced.

This run caught a speaking-script mistake: the demo manifest intentionally
contains an absent `control_03`, so validating it with the junction table is
supposed to fail. The opening script now validates the junction table alone,
then lets the report show the intentional exclusion. The corrected commands
were run successfully. This is developer acceptance, not the independent
beginner or biological review requested in the task guide.
