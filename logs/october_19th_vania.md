# October 19, 2026 — Vania

**Status:** Prepared for push; not yet committed when this log was created

## Contribution

Vania completed Phase 5 Rust task **R10** from `delegations_and_tasks/tasks.md`:
one-gene case/control reading, implemented as a new sample-manifest format,
two new library modules, and a CLI command on top of the existing validated
decomposition engine.

### Sample manifest contract

- Added a `## Sample manifest` section to `docs/data_contract.md` (contract
  version bumped to 0.2): a TSV assigning each `sample_id` to a `case` or
  `control` group. Required columns (`sample_id`, `condition`) are found by
  name rather than a fixed header, since real metadata files (GEO,
  ArrayExpress) carry many extra columns this tool does not need. `condition`
  values are matched case-insensitively; column names are not.
- States explicitly that this manifest covers group membership only, not
  sample eligibility (tissue/QC filtering remains a separate, later task).

### `src/manifest.rs` (new)

- `Condition` enum (`Case`/`Control`) with a `Display` impl, and
  `SampleRecord { sample_id, condition }`.
- `parse_manifest` / `parse_manifest_text`, mirroring `parser.rs`'s
  file/text split. Validates: required columns present and not duplicated in
  the header, nonempty unique `sample_id`, `condition` normalizes to `case`
  or `control` or the row is rejected with a line-numbered error naming the
  bad value.
- Seven unit tests covering the valid case (with extra ignored columns),
  case-insensitive condition matching, an invalid condition value, a
  duplicate `sample_id`, a missing required header column, a duplicated
  required header column, and a wrong field count.

### `src/compare.rs` (new)

- `compare_gene(manifest, results, gene_id) -> GeneComparison`: joins the
  manifest against `analyze_junctions`'s output for one gene by `sample_id`,
  preserving sample identity throughout (no pooling).
- Every sample ends up in exactly one of `included` or `excluded`, with one
  of three explicit exclusion reasons: no manifest entry, no junction data
  for this gene (never conflated with a measured zero), or a failed
  decomposition check (`checks.passed == false`) — this last case was not in
  the original task sketch and was added deliberately so a numerically
  unreliable sample is never silently included in a comparison.
- Five unit tests, including one that runs real `parse_manifest_text` +
  `parse_tsv_text` + `analyze_junctions` together across two genes and
  confirms the irrelevant gene never leaks into the comparison.

### CLI wiring (`src/main.rs`)

- New command: `splice-girl gene <GENE_ID> <junctions.tsv> <manifest.tsv>`.
  ROADMAP.md's sketched syntax (`--case dir/ --control dir/`) is
  directory-based; this instead matches the formats actually implemented —
  one junction TSV (already carrying `sample_id`) plus the new manifest TSV.
- New `--start` menu option `[2] Compare one gene across case/control
  samples`, prompting for the gene ID and both paths.
- `print_gene_comparison` reports included/excluded counts, then an
  `included` table (`sample_id`, `condition`, `cycle_fraction`) and an
  `excluded` table (`sample_id`, `reason`). Exit code is `1` only if every
  sample was excluded; `0` otherwise.
- Added `Display` for `Condition` in `manifest.rs` so the CLI prints it
  directly instead of matching on a type it doesn't own.

Validation completed for this contribution:

- `cargo test` passes all 24 tests (12 existing synthetic/math tests plus 7
  new manifest tests plus 5 new compare tests).
- `cargo fmt --check` passes.
- `cargo clippy --all-targets -- -D warnings` passes.
- Manually ran `splice-girl gene` against a hand-built fixture exercising all
  three outcomes at once (an included case sample, an included control
  sample shaped like the M11 exon-skipping triangle — which correctly
  reported `cycle_fraction = 0.1111111111`, matching M11's independently
  derived `1/9` — a manifest-only sample excluded for having no junction
  data, and an unrelated gene's sample correctly appearing nowhere) and
  confirmed exit code `0` when samples are eligible and `1` when none are.

## Current scope

This closes R10's library and CLI work. It does not implement R11
(statistics output interface), R12 (ranking), R13 (result exports/run
manifests), or sample eligibility/QC filtering (Phase 7's C1) — group
membership and numerical-check status are the only exclusion criteria
applied here. All manifests and junction tables used for testing are
synthetic/hand-built; this does not touch the real GSE179640 or E-MTAB-15117
data the team has collected.
