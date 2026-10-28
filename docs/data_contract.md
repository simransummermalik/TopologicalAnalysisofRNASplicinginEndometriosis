# Splice Girl input data contract

**Contract version:** 0.5 (2026-10-22) — froze real-data exon-node, coordinate,
and strand conventions.

This version keeps text node identifiers for both synthetic fixtures and real
RNA-seq data. Real junction files use the exact reference, annotation,
coordinate, exon-node, and strand conventions recorded in
[`data_methods.md`](data_methods.md); those choices are not silently inferred
from this five-column file.

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
| `gene_id` | Identifier for the gene containing the junction | Required and nonempty; real data preserves the versioned GENCODE `gene_id` |
| `from` | Upstream exon node in biological transcript order | Required and nonempty |
| `to` | Downstream exon node in biological transcript order | Required and nonempty; must differ from `from` |
| `count` | Read support for this junction | A finite number greater than or equal to zero |

The prototype accepts numeric decimal values in `count`. Real count data may be
integers, but the parser also accepts decimals so normalized or simulated
fixtures can use the same interface.

For real data, each node names the full annotated exon interval:

```text
sequence:exon_start-exon_end:strand
```

The exon coordinates are 1-based and closed. A plus-strand edge normally runs
from lower to higher coordinates; a minus-strand edge runs from higher to lower
coordinates so both follow 5-prime to 3-prime transcript order. Synthetic
fixtures may continue to use short labels such as `E1` because they do not
claim genomic coordinates.

## Grouping and validation rules

- Rows are grouped by the pair `(sample_id, gene_id)`.
- Rows from different samples are never pooled together.
- Each sample/gene group must have a positive total count. A group whose total
  support is zero is rejected because its normalized signal is undefined.
- A duplicate directed edge is rejected within one sample/gene group. For
  example, two `E1 -> E2` rows for the same sample and gene are invalid until
  they are aggregated. The preprocessing converter performs this aggregation
  before writing the final table.
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

A real plus-strand row follows the same contract with coordinate-based exon
nodes:

```text
sample_id	gene_id	from	to	count
SRR36710387	ENSG00000123456.7	chr1:100-199:+	chr1:300-399:+	12
```

On the minus strand, the higher-coordinate exon appears in `from` and the
lower-coordinate exon appears in `to`.

## Invalid examples

Negative counts are rejected:

```text
sample_id	gene_id	from	to	count
sample_1	GENE_1	E1	E2	-1
```

Expected error: `column 'count': expected a finite nonnegative number`.

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
the commands that exist today as frozen interface. Multi-gene descriptive
comparison is implemented; inferential statistics and ranking remain separate
future work so they cannot be confused with the current summaries.

### Implemented commands

| Command | Behavior |
|---|---|
| `splice-girl <junctions.tsv> [--output <folder>]` | direct-mode shorthand for `analyze` |
| `splice-girl analyze <junctions.tsv> [--output <folder>]` | analyzes every sample/gene group and optionally saves a result bundle |
| `splice-girl validate <junctions.tsv> [manifest.tsv]` | audits row, sample, gene, support, manifest-matching, graph, and mathematical-check status without printing every edge |
| `splice-girl gene <GENE_ID> <junctions.tsv> <manifest.tsv> [--output <folder>]` | compares one gene across case/control samples and optionally saves a result bundle |
| `splice-girl compare <junctions.tsv> <manifest.tsv> [--output <folder>]` | compares every observed gene using sample-level descriptive summaries and optionally saves a result bundle |
| `splice-girl visualize <SAMPLE_ID> <GENE_ID> <junctions.tsv> <output.svg>` | renders one sample/gene's splice graph as a standalone SVG |
| `splice-girl --start` | interactive menu wrapping analysis, validation, one-gene comparison, all-gene comparison, and visualization |
| `splice-girl --help`, `-h`, `help` | prints usage and the command list |
| `splice-girl --version`, `-V` | prints the Cargo package version |

ROADMAP.md's sketched `gene` syntax is directory-based
(`--case dir/ --control dir/`); the implemented syntax instead takes one
junction TSV (already carrying `sample_id`) plus the manifest TSV described
above, matching the formats this tool actually parses.

### Not yet implemented

- Inferential p-values or multiple-testing correction.
- A `rank` command applying approved candidate-gene and ranking rules.

The implemented `compare` command is descriptive. It does not claim statistical
significance or rank biological importance.

### Output fields, per command

- **`analyze`** (and its `<junctions.tsv>` alias): one block per sample/gene
  group, with `sample_id`, `gene_id`, `nodes`, `edges`, `signal_energy`,
  `gradient_energy`, `cycle_energy`, `cycle_fraction`, a per-edge component
  table (`from`, `to`, `signal`, `gradient`, `cycle`), and a
  `checks: PASS|FAIL` line reporting the reconstruction, kernel, and
  orthogonality error magnitudes.
- **`validate`**: compact counts for rows, samples, genes, sample/gene groups,
  total support, group-support range, and mathematical checks. If a manifest
  is supplied it also reports case/control counts and sample IDs missing from
  either input. A mismatch produces a failure exit code.
- **`gene`**: `gene_id`, the included/excluded sample counts, an `included`
  table (`sample_id`, `condition`, `cycle_fraction`) and an `excluded` table
  (`sample_id`, `reason`).
- **`compare`**: one deterministic row per observed gene containing status,
  case/control sample counts, group means, group medians, case-minus-control
  mean difference, and exclusion count. Missing groups are `NA`; they are
  never converted to zero. No p-value, q-value, or ranking is calculated.
- **`visualize`**: not a tabular export — writes an SVG containing the same
  node/edge counts, cycle fraction, and checks status; aligned observed,
  gradient, and cycle graph panels; a table with all three values per junction;
  and a fixed caveat stating that cycle-space structure alone does not
  establish circular RNA. The panels preserve node positions, edge order,
  stored direction, and one absolute-value width scale. Solid lines are
  positive and dashed lines are negative. Standard output only confirms the
  saved path.

### Saved result bundle

The `analyze`, `gene`, and `compare` commands keep their terminal output by
default. Add `--output <folder>` (or `-o <folder>`) to also save
machine-readable files.
The folder must not already exist; Splice Girl refuses to overwrite a previous
run.

| File | Contents |
|---|---|
| `sample_metrics.tsv` | One row per analyzed sample/gene, including graph size, raw total support, energies, cycle fraction, and numerical-check errors |
| `junction_components.tsv` | One row per junction, including raw count, normalized signal, gradient component, and signed cycle component |
| `comparisons.tsv` | `gene` and `compare`: one row per included or excluded sample/gene decision, with condition, cycle fraction, status, and exclusion reason |
| `gene_summaries.tsv` | `compare` only: one descriptive row per gene with group sizes, means, medians, ranges, case-minus-control effect, status, and interpretation note |
| `run_manifest.json` | Run ID, exact command arguments, tool version, input paths and SHA-256 checksums, solver settings, tolerance, output filenames, counts, exclusions, and failures |

The manifest records `random_seed` as `null` because the current commands do
not use randomness. A later statistics command must record its actual seed.
The TSVs include the same `run_id` so rows can be joined back to the manifest.

### Errors and exit codes

- Every error is printed to standard error as `error: <message>`. Exit code
  is `0` on success and nonzero on any error.
- Input errors name the file. Row errors also name the line and column when
  that information is available, show the invalid value when safe, and suggest
  the next correction.
- A missing-argument or too-many-argument error prints the full multi-line
  usage/help text as the error body, not a short one-line usage string. This
  is deliberate, not a bug: the full command list stays visible whenever
  invocation fails.
- `analyze`, `gene`, `validate`, and `compare` can print useful output to standard
  output and then still exit with a failure status: `analyze` prints each
  sample/gene's full block before checking that group's `checks.passed`, and
  `gene` prints the full `included`/`excluded` tables before failing if every
  sample was excluded or if no eligible case or control sample remains.
  `validate` prints its audit before failing on a mismatch, and `compare`
  prints every gene before failing when no gene has both groups. Success on
  these commands is not all-or-nothing on standard output; a
  caller scripting against them must check the exit code, not assume failure
  means no output was produced.
- `visualize` refuses to overwrite an existing output file and fails with a
  specific error naming the path, rather than silently replacing a figure.
- `--output` also refuses to use an existing folder. Choose a new run folder
  instead of deleting or silently replacing prior results.
- `tests/cli_errors.rs` checks missing and empty files, malformed headers,
  invalid counts and manifests, missing sample metadata, zero-support groups,
  insufficient comparison groups, protected output folders, and nonzero exit
  statuses. The matrix unit test separately checks safe solver failure.

## Future extensions

The current real-data [sample manifest](../data/metadata/samples.tsv) records
subject, condition, tissue, accession, batch, file location, QC status, genome
build, annotation version, and eligibility separately from the five-column
junction table. It preserves the same run-based `sample_id` so each future
junction row can be traced back to its sample.
