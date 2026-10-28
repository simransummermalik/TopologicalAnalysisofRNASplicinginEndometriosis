# Splice Girl demo package

Small, hand-built files for demonstrating the tool in seconds, not minutes.
Nothing here is real biological data — see `docs/data_contract.md` and
`tests/fixtures/expected.md` for the reviewed synthetic fixtures this demo
reuses the shapes of.

## Part 1 — core math on known shapes

These already live in `tests/fixtures/` and are covered by `cargo test`; running
them directly shows the decomposition output a reader can sanity-check by eye.

```sh
cargo run -- tests/fixtures/linear.tsv      # a chain: cycle_fraction 0
cargo run -- tests/fixtures/cycle.tsv       # a pure directed cycle: cycle_fraction 1
cargo run -- tests/fixtures/exon_skip.tsv   # acyclic directed graph, nonzero cycle_fraction (1/9)
```

The `exon_skip.tsv` case is the one worth narrating: there is no directed path
back to `E1`, but the *undirected* shape is a triangle, so the decomposition
still finds cycle-space signal. See `docs/math_spec.md` for the derivation.

## Part 2 — validation fails loudly

```sh
cargo run -- tests/fixtures/zero_support.tsv    # error: positive total support required
cargo run -- tests/fixtures/duplicate_edge.tsv  # error: duplicate edge
cargo run -- tests/fixtures/self_loop.tsv       # error: self-loop not supported
```

## Part 3 — case/control comparison (R10)

`junctions.tsv` and `manifest.tsv` in this folder are built specifically for
this demo, reusing the same known shapes from Part 1 so the numbers are
recognizable:

- `case_01` is shaped like `cycle.tsv` (pure circulation) → `cycle_fraction = 1.0`
- `case_02` is shaped like `exon_skip.tsv` → `cycle_fraction = 1/9`
- `control_01` and `control_02` are tree-shaped (a chain and a branch) → `cycle_fraction = 0.0`
- `control_03` is listed in the manifest but has no junction rows for `GENE_DEMO`
- `orphan_01` has junction rows for `GENE_DEMO` but no manifest entry

```sh
cargo run -- gene GENE_DEMO examples/demo/junctions.tsv examples/demo/manifest.tsv
```

Expected output:

```text
gene: GENE_DEMO
included: 4 sample(s), excluded: 2 sample(s)

included:
  sample_id	condition	cycle_fraction
  case_01	case	1.0000000000
  case_02	case	0.1111111111
  control_01	control	0.0000000000
  control_02	control	0.0000000000

excluded:
  sample_id	reason
  control_03	no junction data for gene 'GENE_DEMO'
  orphan_01	no manifest entry for this sample_id
```

`control_03` and `orphan_01` demonstrate the two distinct exclusion reasons
`compare_gene` (`src/compare.rs`) reports — nothing is silently dropped in
either direction.

Save the same run as reusable files:

```sh
cargo run -- gene GENE_DEMO examples/demo/junctions.tsv examples/demo/manifest.tsv \
  --output demo_results
```

This creates `sample_metrics.tsv`, `junction_components.tsv`,
`comparisons.tsv`, and `run_manifest.json`. Use a new output-folder name for
each run because Splice Girl protects existing results from being overwritten.

## Part 4 — graph visualization

Create a browser-ready picture of `case_02`'s exon-skipping graph:

```sh
cargo run -- visualize case_02 GENE_DEMO examples/demo/junctions.tsv splice-girl-demo.svg
```

Open `splice-girl-demo.svg` in any browser. The plain academic figure shows the
observed, gradient, and cycle signals in three aligned versions of the same
graph. Edge IDs and node positions stay fixed, line width uses one shared
absolute-value scale, and dashed lines mark negative components. The table
reports the read support and all three mathematical values for each junction.
Running the command with another valid sample or gene generates a figure from
that selected result. Splice Girl refuses to overwrite a figure that already
exists, so choose a new output name or remove the old demo SVG before running
the command again.

## Part 5 — compact validation audit

Validate the junction file through the same parser, graph construction, Hodge
solver, and mathematical checks used by the analysis:

```sh
cargo run -- validate examples/demo/junctions.tsv
```

That command passes. Adding the manifest intentionally produces a failed audit
because this demo contains two visible mismatch examples:

```sh
cargo run -- validate examples/demo/junctions.tsv examples/demo/manifest.tsv
```

The report names `orphan_01` as junction data without metadata and
`control_03` as metadata without junction data. A failed validation command
returns a nonzero exit code after printing the audit.

## Part 6 — all-gene descriptive comparison

`all_genes.tsv` contains two synthetic genes with opposite group patterns.
Compare every gene in one command:

```sh
cargo run -- compare examples/demo/all_genes.tsv examples/demo/manifest.tsv
```

`GENE_HIGH_CYCLE` has a case-minus-control mean difference of `+0.5555555556`.
`GENE_LOW_CYCLE` reverses the same graph shapes and has a difference of
`-0.5555555556`. This demonstrates deterministic multi-gene output and effect
direction. Neither value is a p-value or a claim of biological significance.

Save the full five-file result bundle with:

```sh
cargo run -- compare examples/demo/all_genes.tsv examples/demo/manifest.tsv \
  --output demo_all_genes
```

The bundle adds `gene_summaries.tsv` to the sample metrics, junction
components, inclusion/exclusion table, and reproducibility manifest.

## Timing

All demo commands run in well under a second on a development machine (the
startup ASCII banner printing accounts for most of the visible delay).

## Scope

This is an illustrative demo package for showing the tool working end to end.
It is not the real-data demo described in `ROADMAP.md` Phase 11, which will use
an approved candidate gene and the selected GSE315857 samples. That remains a
separate, later milestone.
