# Try Splice Girl after installation

These eight samples and four genes are entirely invented. The example checks
software behavior; it is not an endometriosis result. S0–S3 are cases and
S4–S7 controls. No Python or repository checkout is needed for these commands.

Open a terminal **inside this example folder**. All paths below are relative
to your current folder. Keep output folder names new for each run.

## 1. Prepare comparable graphs

```sh
splice-girl prepare-shared \
  --junctions junctions.tsv --manifest samples.tsv --genes genes.tsv \
  --complete-counts --minimum-edge-count 5 --minimum-samples-per-edge 4 \
  --minimum-gene-support 50 --output prepared
```

This keeps an edge if at least four samples have five reads on it. Coverage
requires at least 50 reads on the selected shared graph. `--complete-counts`
declares that every sample's count export is complete, with read-counter
support cutoff 1. On real data, verify that declaration against upstream logs.

Read `prepared/gene_eligibility.tsv`. TRIANGLE is eligible in all eight
samples. TREE has no cycle; LOW_COVERAGE lacks coverage in S7; ABSENT has no
junctions. Their exclusions remain visible instead of becoming fake results.
The numerical thresholds here exercise a study proposal and do not constitute
approval for another dataset.

## 2. Analyze the eligible genes

```sh
splice-girl stats prepared/junctions.tsv prepared/samples.tsv \
  --genes prepared/eligible_genes.tsv --family synthetic-shared \
  --alternative two-sided --design independent --min-per-group 4 \
  --output statistics
```

Open `statistics/report.html` directly in your browser. Only TRIANGLE is
tested. Its controls have cycle fraction approximately 2/3 on the common
triangle, even though their shortcut count is zero. A zero-count shared edge
is part of the mathematical topology, not evidence that the sample used it.

## 3. Check the seven sensitivity settings

```sh
splice-girl sensitivity \
  --junctions junctions.tsv --manifest samples.tsv --genes genes.tsv \
  --complete-counts \
  --primary-edge-count 5 --primary-sample-count 4 --primary-gene-support 50 \
  --minimum-junction-counts 3,5,10 --minimum-samples-per-edge 2,4,6 \
  --minimum-gene-support 20,50,100 \
  --run-stats --family synthetic-shared --alternative two-sided \
  --design independent --min-per-group 4 --output sensitivity
```

Read `sensitivity/sensitivity.tsv`. Each scenario records every requested gene,
including exclusions. Three scenarios have no eligible genes and skip
statistics; missing effects remain blank. Each usable scenario has its own
comparison tables and optional statistics/report. Anchor sensitivity needs
separate count exports and is not evaluated by these grids.

For descriptive sensitivity only, omit `--run-stats`, `--family`,
`--alternative`, `--design`, and `--min-per-group` together.

## Help and requirements

```sh
splice-girl prepare-shared --help
splice-girl sensitivity --help
```

These commands start from standardized junction counts. Processing raw FASTQ
still requires the separately documented Cleaver/HISAT2 preprocessing tools.
Select primary settings before inspecting biological case/control effects;
report sensitivity results without selecting the most favorable scenario.
