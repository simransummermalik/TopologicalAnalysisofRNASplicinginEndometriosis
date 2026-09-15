# September 15, 2026 — Summer

**Status:** Splice Girl changes pushed; this log is ready to add  
**Commit:** `bbae10e` — `research log push`

## Contribution

Summer pushed the next Splice Girl user-facing update. Starting the program
now displays the detailed ASCII artwork and opens a small interactive menu
around the one analysis feature currently implemented.

### Startup experience

- Added `assets/splice_girl_banner.txt`, the detailed character-based Splice
  Girl figure used at startup.
- Updated `src/main.rs` to print the banner whenever the executable starts.
- Added terminal-safe soft color for interactive output while keeping redirected
  output as plain text.
- Added `splice-girl --start` support.
- Added a menu with only the current feature:

  ```text
  [1] Analyze a junction TSV
  [q] Quit
  ```

- Added prompts for the TSV path, return-to-menu behavior after analysis, and
  readable errors that allow the user to try again.

The existing direct mode remains available:

```text
splice-girl tests/fixtures/cycle.tsv
```

The interactive mode can be started from the project with:

```text
cargo run -- --start
```

or after a local Cargo install with:

```text
splice-girl --start
```

### Documentation and maintenance

- Updated `README.md` with the startup-menu instructions and the banner asset
  link.
- Updated `logs/september_10_summer.md` to include the startup artwork and
  terminal integration in the September 10 prototype history.
- Fixed a Rust documentation-comment formatting issue in `src/hodge.rs` so the
  strict lint check passes.

### Validation

- `cargo test` passes all four tests.
- `cargo fmt --check` passes.
- `cargo clippy --all-targets -- -D warnings` passes.
- The scripted interactive flow successfully starts with `--start`, analyzes
  `tests/fixtures/cycle.tsv`, reports cycle fraction `1.0000000000`, returns to
  the menu, and exits with `q`.

## Current scope

This push improves startup and access to the existing synthetic junction
analysis. It does not add real RNA-seq data, gene research, case/control
statistics, ranking, or the future `compare` and `rank` subcommands.
