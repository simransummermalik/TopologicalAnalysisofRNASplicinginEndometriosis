# Real biological data: current packaged tool check

Completed November 18, 2026. **Open [the report](checked/report.html).**

This run used real human RNA-seq data, not the invented demo. It reused the
saved alignment from the first **1,000,000 cleaned read pairs of control
SRR36710387**, about 2.02% of that sample. No reads were downloaded or realigned.
The executable came from the tested, installed `splice-girl` crate archive.

## What ran and what passed

1. Check the BAM, annotation map, original junction export and candidate list
   against their previously recorded SHA-256 hashes.
2. Select all observed candidate genes, without choosing genes for a high Hodge
   score: **243 of the 661 candidates**, containing **1,388 junction edges**.
3. Run the current Rust tool's input validation and Hodge analysis.
4. Read the BAM again with a separate junction parser. All **19,983
   read-junction observations** match the exported counts exactly; there are
   **zero edge-count differences**. These observations are not independent
   subjects or read pairs.
5. For every observed gene, build its incidence matrix separately in NumPy and
   solve the rectangular least-squares problem using SVD. All **243 genes pass**
   the reconstruction, cycle-kernel, orthogonality and Rust-agreement checks.
   The largest exported-component difference is **4.96 × 10⁻¹¹**, below the
   declared **10⁻⁸** tolerance.
6. Save plain graphs and an offline report. The graph examples were selected by
   highest read support within the cyclic and acyclic categories, with gene ID
   breaking ties. They are examples, not a ranking of disease associations.

Six graphs have cycle rank 1. Five have a cycle fraction above 10⁻⁸; MTSS2 has
a graph loop but a cycle fraction consistent with numerical zero. This is why
having an alternative path and measuring a circulating signal are different.

| Gene | Cycle fraction | Interpretation in this run |
|---|---:|---|
| TNC | 0.001429 | Nonzero cycle component |
| PCSK5 | 0.030893 | Nonzero cycle component |
| MTSS2 | approximately 0 | A graph loop, with essentially no cycle component |
| MYOF | 0.029077 | Nonzero cycle component |
| MME | 0.173225 | Nonzero cycle component |
| FLNA | 0.085446 | Nonzero cycle component |

These fractions describe the mathematical decomposition on each observed
graph. They do not establish circular RNA or identify endometriosis genes.

## What this establishes

The current installed tool handles this real junction table, including real
graphs with nonzero cycle signals. Its exported numbers agree with a separate
numerical implementation, and a separate BAM recount agrees with every
candidate edge count. This extends the earlier independent check of CDH11,
which had no graph cycles.

Both checks reuse the same reference and underlying sequencing data. They
check implementation consistency, not an external biological ground truth.

## What it does not establish

- This is **one control, with no case samples**. No disease effect, p-value,
  adjusted p-value or case/control ranking was calculated.
- The 418 unobserved candidates remain labeled **not observed in this prefix**.
  They were not given zero cycle fractions or treated as absent from the donor.
- This deliberately reproduces the historical **exon-interval node rule and
  unstranded alignment**, with anchor 8 and minimum junction support 1.
  Historical ambiguity exclusions remain. It does not validate the newer
  RF/shared-splice-site, common-cohort graph workflow.
- A deterministic prefix is not the full sample or a random subsample.

The full eight-sample study still needs the other local read inputs and enough
storage/compute. Approximately 400 MB was free when this check began, so a new
human alignment or full-cohort download was not started.

## Files

- [Offline report](checked/report.html): summary and all 243 observed genes.
- [Independent checks](checked/independent_checks.tsv): per-gene numerical errors,
  graph size, support and cycle fraction.
- [Candidate coverage](checked/candidate_coverage.tsv): all 661 candidates,
  including those not observed.
- [Recount differences](checked/recount_differences.tsv): header only, because
  every count agrees.
- [FLNA graph](checked/cyclic.svg) and [COL1A1 graph](checked/acyclic.svg).
- [Rust metrics](checked/rust/sample_metrics.tsv) and
  [signed edge components](checked/rust/junction_components.tsv).
- [Validation manifest](checked/validation.json): commands, hashes, limits,
  actual executable identity and outputs.

## Repeat the exact check

The report opens without installing anything. Recomputing it needs the local
historical BAM/map/export, a Python environment with NumPy, samtools and the
reviewed installed binary. Large raw files, alignments and indexes stay outside
Git. The script refuses an existing output folder.

From the repository root, using the same tools as this run:

```sh
python data/pilots/SRR36710387_real_validation_nov18/reproduce.py \
  --binary target/crate-busy-fix/installed/install/bin/splice-girl \
  --samtools /tmp/splice-girl-hisat2-tools/bin/samtools \
  --output-dir data/pilots/SRR36710387_real_validation_nov18/repeat-01
```

Use a new output name for each repeat. Another installation can supply its own
binary/samtools paths; their actual hashes and commands are recorded each time.
Nothing was committed, pushed or published by the assistant.
