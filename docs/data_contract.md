# Splice Girl input data contract

**Contract version:** 0.4 (2026-10-22) — added machine-readable result bundles.

This version uses text node identifiers for the synthetic prototype. When real
STAR junctions are added, the genome build, annotation release, and genomic
coordinate convention must be recorded in the separate sample/run metadata;
those reference choices are not silently inferred from this five-column file.

This document defines the small TSV formats accepted by the current Splice
Girl prototype and the command-line interface built on top of them. It does
not replace the final biological metadata schema that will be needed for real
RNA-seq analysis.

## Required header

Every input file must begin with this exact tab-separated header, in this
order:

```text
sample_id	gene_id	from	to	count
```

Each later row describes one observed directed splice junction for one sample
and one gene.

| Column | Meaning | Rule |
|---|---|---|
| `sample_id` | Identifier for one biological sample or library | Required and nonempty; keep samples separate |
| `gene_id` | Identifier for the gene containing the junction | Required and nonempty; use a stable ID when available |
| `from` | Start node of the directed junction | Required and nonempty |
| `to` | End node of the directed junction | Required and nonempty; must differ from `from` |
| `count` | Read support for this junction | A finite number greater than or equal to zero |

The prototype accepts numeric decimal values in `count`. Real count data may be
integers, but the parser also accepts decimals so normalized or simulated
fixtures can use the same interface.

## Grouping and validation rules

- Rows are grouped by the pair `(sample_id, gene_id)`.
- Rows from different samples are never pooled together.
- Each sample/gene group must have a positive total count. A group whose total
  support is zero is rejected because its normalized signal is undefined.
- A duplicate directed edge is rejected within one sample/gene group. For
  example, two `E1 -> E2` rows for the same sample and gene are invalid until
  an aggregation rule is approved.
- Self-loops such as `E1 -> E1` are rejected by the current graph model.
- Blank lines are ignored, but a file with no junction rows is rejected.
- Node and edge order is sorted by the text identifiers. Reordering valid input
  rows must not change the result.

The parser reports the line number for malformed rows. The command returns a
failure status instead of silently dropping invalid data.

## Valid example

This is the complete content of `tests/fixtures/cycle.tsv`:

```text
sample_id	gene_id	from	to	count
sample_cycle	GENE_CYCLE	E1	E2	100
sample_cycle	GENE_CYCLE	E2	E3	100
sample_cycle	GENE_CYCLE	E3	E1	100
```

It contains one sample/gene group with total support 300. The normalized edge
signals are each `100 / 300`, and the equal-weight directed cycle should have a
cycle fraction near one.

The chain fixture is also valid:

```text
sample_id	gene_id	from	to	count
sample_linear	GENE_LINEAR	E1	E2	100
sample_linear	GENE_LINEAR	E2	E3	80
```

This graph has no cycle space, so its cycle fraction should be near zero.

## Invalid examples

Negative counts are rejected:

```text
sample_id	gene_id	from	to	count
sample_1	GENE_1	E1	E2	-1
```

Expected error: `count must be finite and nonnegative`.

A self-loop is rejected:

```text
sample_id	gene_id	from	to	count
sample_1	GENE_1	E1	E1	10
```

Expected error: `self-loop E1 -> E1 is not supported`.

Two identical directed edges in one sample/gene group are rejected:

```text
sample_id	gene_id	from	to	count
sample_1	GENE_1	E1	E2	10
sample_1	GENE_1	E1	E2	5
```

Expected error: `duplicate edge E1 -> E2 in one sample/gene group`.

## How to check a file

From the repository root, run:

```sh
cargo test
cargo run -- tests/fixtures/linear.tsv
cargo run -- tests/fixtures/cycle.tsv
```

The automated tests check reconstruction, cycle-space membership, and
orthogonality. The chain and cycle commands provide quick human-readable
checks while the input and output format are still small.

## Sample manifest

A separate TSV assigns each `sample_id` used in the junction table to a case
or control group. This file is independent of the junction table's own format
and is not required for the basic synthetic prototype — only for
sample-preserving comparisons (R10).

Required columns, found by name (extra columns are permitted and ignored;
column order does not matter):

| Column | Meaning | Rule |
|---|---|---|
| `sample_id` | Identifier matching the `sample_id` values used in the junction table | Required, nonempty, unique within the manifest |
| `condition` | Case/control group | Required; normalizes case-insensitively to `case` or `control`; any other value is rejected |

Unlike the junction table's header, the manifest header is matched by column
name rather than an exact fixed header, since real metadata files (for
example GEO or ArrayExpress sample sheets) carry many additional columns this
tool does not need.

The Rust parser uses this manifest for group membership only. The current
real-data manifest at [`data/metadata/samples.tsv`](../data/metadata/samples.tsv)
also carries tissue, QC, processing, and eligibility fields for project
tracking, but the parser ignores those extra fields. Until eligibility
filtering is implemented, only approved samples should be present in the
junction table passed to a comparison.

## Rust command-line interface

The package and command name are settled: `splice-girl`. This section records
the commands that exist today as frozen interface; the multi-gene `compare`
and `rank` commands sketched in ROADMAP.md remain future work and are listed
separately below so the two are never confused.

### Implemented commands

| Command | Behavior |
|---|---|
| `splice-girl <junctions.tsv> [--output <folder>]` | direct-mode shorthand for `analyze` |
| `splice-girl analyze <junctions.tsv> [--output <folder>]` | analyzes every sample/gene group and optionally saves a result bundle |
| `splice-girl gene <GENE_ID> <junctions.tsv> <manifest.tsv> [--output <folder>]` | compares one gene across case/control samples and optionally saves a result bundle |
| `splice-girl visualize <SAMPLE_ID> <GENE_ID> <junctions.tsv> <output.svg>` | renders one sample/gene's splice graph as a standalone SVG |
| `splice-girl --start` | interactive menu wrapping the three commands above |
| `splice-girl --help`, `-h`, `help` | prints usage and the command list |
| `splice-girl --version`, `-V` | prints the Cargo package version |

ROADMAP.md's sketched `gene` syntax is directory-based
(`--case dir/ --control dir/`); the implemented syntax instead takes one
junction TSV (already carrying `sample_id`) plus the manifest TSV described
above, matching the formats this tool actually parses.

### Not yet implemented

- A multi-gene `compare` across case/control sample directories.
- A `rank` command applying candidate-gene membership and a ranking rule.

Neither exists yet; nothing in this document should be read as describing
their behavior.

### Output fields, per command

- **`analyze`** (and its `<junctions.tsv>` alias): one block per sample/gene
  group, with `sample_id`, `gene_id`, `nodes`, `edges`, `signal_energy`,
  `gradient_energy`, `cycle_energy`, `cycle_fraction`, a per-edge component
  table (`from`, `to`, `signal`, `gradient`, `cycle`), and a
  `checks: PASS|FAIL` line reporting the reconstruction, kernel, and
  orthogonality error magnitudes.
- **`gene`**: `gene_id`, the included/excluded sample counts, an `included`
  table (`sample_id`, `condition`, `cycle_fraction`) and an `excluded` table
  (`sample_id`, `reason`).
- **`visualize`**: not a tabular export — writes an SVG containing the same
  node/edge counts, cycle fraction, and checks status, the drawn graph, a
  per-junction table with signed cycle components, and a fixed caveat line
  stating that cycle-space structure alone does not establish circular RNA.
  Standard output only confirms the saved path.

### Saved result bundle

The `analyze` and `gene` commands keep their terminal output by default. Add
`--output <folder>` (or `-o <folder>`) to also save machine-readable files.
The folder must not already exist; Splice Girl refuses to overwrite a previous
run.

| File | Contents |
|---|---|
| `sample_metrics.tsv` | One row per analyzed sample/gene, including graph size, raw total support, energies, cycle fraction, and numerical-check errors |
| `junction_components.tsv` | One row per junction, including raw count, normalized signal, gradient component, and signed cycle component |
| `comparisons.tsv` | `gene` only: one row per included or excluded sample, with condition, cycle fraction, status, and exclusion reason |
| `run_manifest.json` | Run ID, exact command arguments, tool version, input paths and SHA-256 checksums, solver settings, tolerance, output filenames, counts, exclusions, and failures |

The manifest records `random_seed` as `null` because the current commands do
not use randomness. A later statistics command must record its actual seed.
The TSVs include the same `run_id` so rows can be joined back to the manifest.

### Errors and exit codes

- Every error is printed to standard error as `error: <message>`. Exit code
  is `0` on success and nonzero on any error.
- A missing-argument or too-many-argument error prints the full multi-line
  usage/help text as the error body, not a short one-line usage string. This
  is deliberate, not a bug: the full command list stays visible whenever
  invocation fails.
- `analyze` and `gene` can both print complete, correct output to standard
  output and then still exit with a failure status: `analyze` prints each
  sample/gene's full block before checking that group's `checks.passed`, and
  `gene` prints the full `included`/`excluded` tables before failing if
  `included` turned out to be empty. Success on these two commands is not
  all-or-nothing on standard output; a caller scripting against them must
  check the exit code, not assume failure means no output was produced.
- `visualize` refuses to overwrite an existing output file and fails with a
  specific error naming the path, rather than silently replacing a figure.
- `--output` also refuses to use an existing folder. Choose a new run folder
  instead of deleting or silently replacing prior results.

## Future extensions

The current real-data [sample manifest](../data/metadata/samples.tsv) records
subject, condition, tissue, accession, batch, file location, QC status, genome
build, annotation version, and eligibility separately from the five-column
junction table. It preserves the same run-based `sample_id` so each future
junction row can be traced back to its sample.
