# Shared graphs and threshold sensitivity

Implemented November 18, 2026. The software implements the
[threshold proposal](threshold_proposal.md) with explicit caller settings.
Implementation and synthetic checks do not approve the proposed values for
the biological study. Choose and record the primary settings before seeing
case/control effects; report the declared sensitivity grid without choosing
the most favorable result afterward.

## What the extra step does

The existing Rust commands build a graph from every row supplied for one
sample/gene, including rows with count zero. `scripts/shared_graph.py` prepares
those rows so covered samples share the same nodes and edges for each gene:

1. Read complete counts from every cohort sample and the requested gene list.
2. For each edge of a gene, count how many samples meet the edge read cutoff.
   Select the edge if enough samples support it. Condition labels are not used.
3. Calculate cycle-space dimension as edges minus nodes plus connected
   components. A value of at least one is required; a directed return path
   is not required. A simple exon-skipping graph qualifies.
4. Sum each sample's counts **only on the selected edges**. Counts on rejected
   edges cannot make a sample pass the coverage threshold.
5. For covered cyclic sample/gene graphs, write all selected edges. A selected
   edge absent from that sample's complete export gets count zero. Below-
   coverage sample/gene graphs get no analysis rows and an exclusion reason.
6. Write an eligible gene family containing only cyclic genes for which every
   manifest sample passes coverage. The original family annotations are kept.
7. Pass these tables to the existing Rust math, comparison, reporting and
   statistics commands. The math and original input files are not rewritten.

The eligibility rule is deliberately stricter than a minimum group size:
every manifest sample must be covered. With the study's four cases and four
controls, that implements the required common eight-sample cohort.

## Before using real counts

- Use the combined five-column `junctions.tsv` from successful exports of
  **all** intended samples, with read-counter minimum junction support 1.
- Confirm matching reference, node rule, RF setting and anchor/filter settings
  across exports. Rebuild the map for the adopted shared-splice-site node rule
  before counting; do not mix the old pilot's node IDs with new counts.
- Use only approved samples in the manifest; extra metadata is preserved.
  The tool rejects junction samples absent from the manifest and manifest
  samples with no positive input counts anywhere, which may indicate a
  missing or failed export. A truly empty sample requires upstream review.
- Supply `--complete-counts` only after checking completed upstream manifests
  and counts. This is a caller declaration, not automatic proof of complete
  sequencing, unbiased detection, or absence of the RNA molecule.
- Counts must be nonnegative integers; duplicate sample/gene edges, self-loops,
  malformed inputs and invalid thresholds are rejected before exporting.
- Pass the complete candidate family with `--genes` to retain absent genes in
  the audit. Without it, the audit contains only originally observed genes.

The missing seven local sample pairs and storage/compute limitations remain
separate prerequisites for the full cohort. The original one-control pilot
cannot establish an eight-sample shared graph by itself.

## Start with the small example

The [eight-sample synthetic walkthrough](../examples/shared_graph/README.md)
includes copyable commands and expected answers. It downloads nothing and
uses no human reads. The scripts use Python's standard library and the
existing Rust executable, built with `cargo build --release`.

For a real approved cohort, replace the example paths. The parameter structure
is:

```sh
python3 scripts/shared_graph.py \
  --junctions COMPLETE_COHORT_COUNTS.tsv --manifest COHORT_SAMPLES.tsv \
  --genes data/gene_sets/candidate_genes.tsv --complete-counts \
  --minimum-edge-count 5 --minimum-samples-per-edge 4 \
  --minimum-gene-support 50 --output-dir target/cohort_shared_review
```

These numbers demonstrate the proposal; they are not silently installed as
production defaults. Every threshold is a required argument. The script
protects existing output folders and records checksums, settings and the
completion/failure status in `run_manifest.json`.

## Files to read

| File | Meaning |
|---|---|
| `edges.tsv` | Every observed candidate edge, how many samples pass its support cutoff, and whether it was selected. |
| `sample_eligibility.tsv` | Every requested gene × manifest sample; retained support, coverage pass, cycle rank, zero-edge count and status. |
| `gene_eligibility.tsv` | Shared edge count, cycle rank, covered/required sample counts and final eligibility, including absent genes. |
| `junctions.tsv` | Covered cyclic sample/gene graphs, retaining selected zero-count edges. This may contain partially covered genes for descriptive review. |
| `samples.tsv` | A copy of the full original sample manifest, including subject metadata. |
| `eligible_genes.tsv` | Fully covered cyclic genes only, retaining input family columns. Use this family for the study's statistics handoff. |
| `run_manifest.json` | Input/output SHA-256 checksums, implementation checksum, explicit settings and status. |

An all-excluded result is a valid preparation outcome. The junction and family
files then contain only headers: inspect the eligibility tables and skip Rust
analysis instead of inventing values. Rust rejects an empty input table.

When partially covered genes remain, ordinary `compare` can show descriptive
results for their covered samples. They are excluded from `eligible_genes.tsv`.
Do not substitute the unfiltered original candidate list or reduce the sample
minimum to make them eligible for the predeclared eight-sample test.

## Seven sensitivity scenarios

`scripts/sensitivity.py --mode shared` runs the explicit primary setting,
then changes only one setting at a time. With the proposed grids
`3,5,10` edge counts, `2,4,6` supporting samples and `20,50,100` coverage,
there are seven unique scenarios. The primary is not repeated when it occurs
inside the grids. The old per-sample Cartesian grid remains available in
the default `--mode sample` for reproducing earlier exploratory work.

Each shared scenario saves its prepared tables, Rust comparison bundle and
logs. The top-level `sensitivity.tsv` retains every requested gene in every
scenario, including those with no shared edges or cycles. An effect is blank
unless the gene has full cohort coverage and a cyclic shared graph. Different
scenarios may test different subsets; eligible gene counts and identities
must be reported with any changes in statistics.

Statistics are opt-in through `--run-stats` plus explicit `--genes`, `--family`,
`--alternative`, `--design independent` and `--min-per-group`. The runner
passes each scenario's eligible family to Rust, preserving annotations, exact
gene-set tests, gene-level BH and reports. It skips statistics if no genes
qualify. These are separate sensitivity outputs; no automatic correction
across scenarios or selection of a winning scenario is performed. The primary
study direction remains two-sided and the sample minimum four per group.

Anchors cannot be varied using a count table that already discarded short
anchors. An anchor sensitivity check requires separate complete count exports
with each predeclared anchor setting from the same alignments (and settings
recorded). This runner records that anchors were not evaluated.

## Interpretation changes to explain

Cycle fraction is now measured relative to the cohort's shared topology.
A sample with zero reads on a shortcut can have a nonzero cycle-space
component because the shortcut exists in the common graph. This does not
establish that the sample biologically used the shortcut. The toy walkthrough
contains an explicit example whose cycle fraction is 2/3 under this setup.
Do not compare these values directly with the earlier per-sample-topology
pilot as though they were produced by identical analysis definitions.

Selecting edges without condition labels makes the selection unchanged under
permutation of those labels; it does not remove sequencing-depth differences,
cell composition, batch effects or unequal detection. A fixed graph improves
structural comparability while those QC and interpretation checks remain.

The implementation checks include label/row-order invariance, inclusive
threshold boundaries, disconnected cycle rank, measured-zero preservation,
coverage exclusions, missing exports, all-excluded scenarios and actual Rust
statistics handoff. See `scripts/test_shared_graph.py` and the dated
[software validation record](software_validation.md).
