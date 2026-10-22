# October 22, 2026 — Summer

**Status:** Pushed

**Commit:** `a9ed39e` — `adding output.rs`

## Contribution

Today I completed Rust task **R13**, which makes Splice Girl save its analysis
instead of leaving the results only in terminal text. The existing commands
still work the same way, but `analyze` and `gene` now accept an optional
`--output <folder>` argument. A saved analysis includes sample-level Hodge
measurements, every junction's normalized, gradient, and cycle components, and
a JSON run manifest. A `gene` run also saves the included and excluded samples
in a comparison table.

The run manifest records the exact command, Splice Girl version, input paths,
SHA-256 input checksums, normalization and solver settings, numerical
tolerance, output filenames, exclusions, and failed mathematical checks. Each
TSV row carries the same run ID as the manifest. The program refuses to use an
existing output folder so an old result is not silently overwritten.

I also updated the README, demo guide, and data contract with the new command
and file descriptions. Finally, I audited `delegations_and_tasks/tasks.md`
against the files currently in the repository. The status board now marks the
completed setup, math, Rust, usability, pipeline-figure, and demo tasks and
names the actual blocker for unfinished real-data tasks.

## Validation

- `cargo test` passes: **35 tests passed, 0 failed**.
- `cargo clippy --all-targets -- -D warnings` passes.
- The `analyze --output` command produced valid sample and junction TSVs plus
  a JSON manifest from the directed-cycle fixture.
- The `gene --output` command produced all four result files from the demo and
  correctly recorded four included and two excluded samples.
- Python's JSON and TSV readers successfully parsed every generated file.
- Known SHA-256 test vectors for empty input and `abc` pass.

The completed R13 work and task-board audit are recorded in commit `a9ed39e`.
