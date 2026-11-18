# Run Splice Girl from start to finish

Run commands from the repository root. No command here pushes to GitHub.
Result folders must be new; change their names for a repeat run.

## First: a complete demonstration with synthetic data

```sh
cargo build --release
cargo run --release -- stats examples/demo/all_genes.tsv \
  examples/demo/manifest.tsv --genes examples/demo/stats_genes.tsv \
  --family synthetic-demo --alternative two-sided --design independent \
  --min-per-group 2 --rank-by absolute-effect --output target/demo_complete
```

Open `target/demo_complete/report.html` in a browser. This needs no web server
or internet. It uses made-up junction counts for two genes and four measured
samples; the fifth manifest sample is deliberately missing and remains an
explicit exclusion. This demonstrates the software, not a disease finding.
The two raw and adjusted p-values are both 1/3; the large-looking toy effects
are not evidence of statistical significance.

## What the output files tell you

| File | Question it answers |
|---|---|
| `report.html` | What happened, including samples, graphs, tests and exclusions? |
| `sample_metrics.tsv` | What was each sample/gene's Hodge measurement? |
| `junction_components.tsv` | Which exact junctions contributed to each component? |
| `comparisons.tsv` | Which samples were included or excluded, and why? |
| `gene_summaries.tsv` | What are the case/control means and their difference? |
| `statistics.tsv` | What are the exact permutation p-values and BH-adjusted values? |
| `set_test.tsv` | Does the requested family as a whole differ? One exact gene-set test (mean absolute effect when two-sided) |
| `set_scores.tsv` | What is each sample's mean cycle fraction over the genes in the set test? |
| `baseline_metrics.tsv` | What are each eligible sample's simpler graph/count measurements? |
| `baseline_comparisons.tsv` | Do simpler measurements also differ, and do observed topologies match? |
| `gene_family.tsv` | Which original gene list and evidence annotations were supplied? |
| `rankings.tsv` | In what explicitly requested order should tested genes be inspected? |
| `run_manifest.json` | Which files, hashes, settings and command produced this run? |

Baselines are raw junction support, observed nodes, junctions, and cycle rank
(edges minus nodes plus connected components). They use the same included
samples as the Hodge comparison. Missing measurements remain missing. Raw
support is sensitive to sequencing depth and is not a gene-expression assay.
`different_review_required` flags differing observed graph topologies; this is
a prompt for coverage/annotation review, not proof of biological change.

Ranking is optional. Omit `--rank-by` to retain alphabetical tables and create
no ranking file. `absolute-effect` sorts tested genes by effect magnitude;
`adjusted-p` sorts by adjusted p-value, then effect magnitude. Exact ties use
gene ID. Untested genes are retained without ranks. Ordering never selects a
significance cutoff or approves a gene family.

When the supplied gene-family table contains symbols, risk/environment flags
or categories, these now appear directly in the statistics and ranking tables,
joined by exact versioned gene ID. The report shows the supplied symbol and
category too. Missing annotations stay blank; the full original family table
is still preserved for evidence tracing.

## Check robustness to support thresholds

```sh
python3 scripts/sensitivity.py \
  --junctions examples/demo/all_genes.tsv \
  --manifest examples/demo/manifest.tsv \
  --minimum-junction-counts 1,30 --minimum-gene-support 0,150 \
  --output-dir target/demo_sensitivity
```

These four combinations are demonstration settings. The tool filters junctions,
then checks total retained support per sample/gene and reruns the Rust Hodge
analysis. `sensitivity.tsv` retains every originally observed gene in every
scenario, including genes with no retained data; missing effects stay blank.
Each scenario saves its filtered inputs, logs and available result bundle;
the unfiltered reference bundle is saved too. The run manifest records the
grid, binary/input hashes and output checksum. This is descriptive: it does
not select the most favorable threshold, perform extra significance tests,
or replace rerunning alignments to test anchor/aligner sensitivity.

## Then: real reads

Follow the [Cleaver/HISAT2 guide](../preprocessing/HISAT2_GUIDE.md) to prepare
the verified reference and process a local sample sheet. That produces the
same junction/manifest contract as the demo. Inspect read and alignment QC,
sample labels, exclusions and graph comparability before inference.

For smaller intermediate files, add `--storage-mode stream-bam` to the
preprocessing plan/run command. It pipes alignments into sorted BAM without
saving SAM. The adapter checks estimated disk space before starting and
retains stage logs/partial files after failures or interruptions. The
[storage guide](../preprocessing/HISAT2_GUIDE.md#save-disk-space-with-streaming-bam)
explains the reserve, memory controls and retry behavior. No existing data
are automatically deleted.

The gene evidence helpers are documented in
[gene_set_methods.md](gene_set_methods.md). The repository now contains
source-backed risk/environment extractions and a 661-gene candidate table;
see [the evidence status](../data/gene_sets/README.md) for completed reviews
and remaining independent spot checks. The software preserves the supplied
gene annotations; it does not invent genetic-risk or environmental-response
evidence.

## One command to check the software

```sh
python3 scripts/check_project.py
```

This runs formatting, Rust tests, strict Clippy, release build, preprocessing
tests, gene-set tests, sensitivity tests and the synthetic benchmark.
Add `--with-tools` to also execute the tiny actual Cleaver/HISAT2 batch test.
External tools must be installed first. The default check needs Cargo and
Python, with no Galaxy account or full human dataset.

## What still requires the group

The software can process inputs, but the real eight-sample analysis is not
finished: the full local FASTQ collection, reviewed preprocessing/coverage
choices and comparable cohort junction tables are still needed. The human
index, candidate evidence tables and bounded CDH11 pilot now exist; independent
evidence sign-offs remain. The inference policy is recorded in the statistics
specification, while its real-data assumptions still need review. The
[assistant-run beginner walkthrough](../beginner_test/README.md) records the
installation/demo checks; an independent human usability review, final biological
interpretation, license choice and final presentation are separate remaining
tasks. No final disease or exposure claim follows from the synthetic demo.
