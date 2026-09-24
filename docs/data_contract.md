# Splice Girl input data contract

This document defines the small TSV format accepted by the current Splice Girl
prototype.It does not replace the final biological metadata schema that will be
needed for real RNA-seq analysis.

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

## Future extensions

Real-data analysis will add separate metadata files for subject, condition,
tissue, accession, batch, genome build, and annotation version. Those fields
should not be squeezed into this five-column junction table. The future
metadata schema must preserve the same `sample_id` so each junction row can be
traced back to its sample.
