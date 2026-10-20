# October 20, 2026 — Vania

**Status:** Pushed

**Commit:** `49a9542` — `Finished R1: froze the Rust CLI interface in the data contract`

## Contribution

Vania completed Phase 5 Rust task **R1** from `delegations_and_tasks/tasks.md`:
freezing the Rust command-line interface as documented, reviewed interface in
`docs/data_contract.md`.

### Why this was more than originally scoped

R1 was expected to be a short transcription task, but the CLI had grown
substantially since R10 (not through this work — a teammate had added
`analyze` as an explicit alias, `--help`/`-h`/`help`, `--version`/`-V`, and an
entirely new `visualize` command with its own SVG renderer in
`src/visualize.rs`). Freezing the interface meant documenting five real
commands, not two, and none of that growth had a written record anywhere.

### `docs/data_contract.md` (contract version bumped to 0.3)

- Broadened the document's opening sentence, since it previously claimed to
  define only the TSV format and now also covers the CLI.
- Added a `## Rust command-line interface` section with:
  - An **implemented commands** table (`<junctions.tsv>`, `analyze`, `gene`,
    `visualize`, `--start`, `--help`/`-h`/`help`, `--version`/`-V`), noting
    explicitly that ROADMAP.md's directory-based `gene` syntax
    (`--case dir/ --control dir/`) is aspirational, not what is implemented.
  - A **not yet implemented** list (`compare`, `rank`) so the two are never
    confused with real behavior.
  - **Output fields per command**: `analyze`'s per-sample/gene block fields,
    `gene`'s included/excluded tables, and `visualize`'s SVG contents (stated
    explicitly as not a tabular export).
  - **Errors and exit codes**: the `error:`/exit-code convention, and two
    behaviors worth documenting precisely because they could surprise a
    caller — `analyze` and `gene` can both print complete, correct output to
    standard output and still exit with failure (success is not
    all-or-nothing on stdout), and `visualize` refuses to overwrite an
    existing output file.

### Verification

Every claim in the new section was checked against the running binary before
being written down, not assumed from reading the source:

- Ran `cargo run -- --help`, confirmed the printed usage matches exactly.
- Ran a deliberately malformed `gene` invocation, confirmed the full usage
  text is printed as the error body and the exit code is `1`.
- Ran `cargo test` to confirm the documentation-only change didn't break
  anything (it didn't — no source files changed).
- Separately audited `README.md` for consistency with the new documentation.
  Everything in the README's "Try Splice Girl" and "What the output means"
  sections matched exactly (`linear`/`branching` → `0`, `cycle` → `1`,
  `exon_skip` → `0.1111...`, the `visualize` command run verbatim from the
  README produces the exact output described). Two discrepancies were found,
  tracked separately rather than fixed in this commit: the README's claim of
  "29 automated tests" is stale (the real count, confirmed by listing every
  test individually, is 31), and the README's "Technical documentation" table
  still describes `docs/data_contract.md` as covering only "accepted TSV
  formats," which is now incomplete.

## Current scope

This freezes documentation of CLI behavior that already existed; it does not
change any command's behavior. It does not implement the `compare` or `rank`
commands. The two README discrepancies found during verification are fixed
in a separate, following contribution.
