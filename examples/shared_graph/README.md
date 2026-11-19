# Eight-sample shared-graph example

All counts, genes and samples here are invented. S0–S3 are labeled case and
S4–S7 control solely to exercise the four-versus-four software workflow.
These are not GSE315857 results.

## Using the installed Rust command

```sh
splice-girl example --output demo-inputs
```

Read `demo-inputs/README.md` and run its commands from that new folder. The
installed `prepare-shared` and `sensitivity` commands include this entire
analysis workflow without Python or a repository checkout. See the
[installed guide](INSTALLED.md) for all commands and expected results.

## Reproduce the earlier Python workflow

From the repository root:

```sh
cargo build --release
python3 scripts/shared_graph.py \
  --junctions examples/shared_graph/junctions.tsv \
  --manifest examples/shared_graph/samples.tsv \
  --genes examples/shared_graph/genes.tsv \
  --complete-counts \
  --minimum-edge-count 5 --minimum-samples-per-edge 4 \
  --minimum-gene-support 50 --output-dir target/shared-demo
```

`--complete-counts` declares that these example inputs contain the complete
counts for all eight samples. On real data, verify the upstream completed
exports and settings before making that declaration.

Open `target/shared-demo/gene_eligibility.tsv` in a text editor or spreadsheet:

| Gene | Expected status | Why |
|---|---|---|
| TRIANGLE | eligible | All eight samples have at least 50 reads on the same three-edge graph. |
| TREE | no_cycle_in_shared_graph | Two connected edges have no cycle space. |
| LOW_COVERAGE | incomplete_coverage | S7 has only one retained read. The 10,000-read edge found in S0 alone does not enter the shared graph. |
| ABSENT | no_shared_edges | The requested gene has no input junctions. |

In `junctions.tsv`, each TRIANGLE sample has the same three edges. S4–S7
have a zero on A→C because that edge had no reads in those otherwise covered
samples. LOW_COVERAGE/S7 has no rows: its result is missing. TREE and ABSENT
also have no analysis rows, and their exclusion reasons remain in the audits.

Run the existing Rust statistics command on the prepared table and the
**eligible** family:

```sh
target/release/splice-girl stats target/shared-demo/junctions.tsv \
  target/shared-demo/samples.tsv --genes target/shared-demo/eligible_genes.tsv \
  --family synthetic-shared --alternative two-sided --design independent \
  --min-per-group 4 --output target/shared-demo-stats
```

Only TRIANGLE is tested. Its control cycle fraction is approximately 2/3 on
the shared graph. This differs from analyzing each control's two observed
edges alone, which would produce zero. The zero shortcut is part of the
common mathematical graph; it does not assert observed biological use of
that connection in the control.

For all seven sensitivity scenarios and optional statistics:

```sh
python3 scripts/sensitivity.py --mode shared \
  --junctions examples/shared_graph/junctions.tsv \
  --manifest examples/shared_graph/samples.tsv \
  --genes examples/shared_graph/genes.tsv --complete-counts \
  --primary-edge-count 5 --primary-sample-count 4 --primary-gene-support 50 \
  --minimum-junction-counts 3,5,10 --minimum-samples-per-edge 2,4,6 \
  --minimum-gene-support 20,50,100 \
  --run-stats --family synthetic-shared --alternative two-sided \
  --design independent --min-per-group 4 \
  --output-dir target/shared-demo-sensitivity
```

Read `sensitivity.tsv` and each scenario's eligibility tables. Increasing the
edge cutoff to 10 or requiring six samples removes TRIANGLE's shortcut and
therefore its cycle. A coverage cutoff of 100 leaves no eligible genes. Those
scenarios are recorded and statistics are skipped rather than invented.

Use a new output folder for a repeat run. See the [full guide](../../docs/shared_graphs.md)
for real-input requirements, file definitions and interpretation limits.
