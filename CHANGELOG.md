# Changelog

## Unreleased

- Exclude paired samples with missing or inconclusive strand checks from
  sequential batches, preserving their files and recording why they were
  excluded. Configured Rust analysis also refuses recorded inconclusive
  checks. The assembly workflow already required a passing check.
- Reject invalid strand agreement thresholds rather than interpreting
  values outside 0--1 or nonfinite values.
- Add focused `--help` and `-h` for analyze, validate, gene, compare, report,
  visualize, stats and cite, including report/SVG viewing instructions.
- Update the crate installation guide to show the published registry command.

## 0.2.1 — 2026-11-19

- Fix the configured-analysis report footer to use the package's declared
  license instead of stale MIT wording.

## 0.2.0 — 2026-11-19

- Set the source release license to Creative Commons Attribution 4.0
  International (CC BY 4.0), permitting commercial and noncommercial reuse
  with attribution under its terms.
- Preserve the published `0.1.0` crate under its original MIT license.
- Update crate metadata, citation metadata, and license notices.
- Record authorization from all four named contributors for the license update.

## 0.1.0 — 2026-11-19

First experimental crates.io release of the Splice Girl Rust analysis CLI.
This release does not claim full-cohort validation or provide the raw FASTQ
preprocessing tools in the installed crate; those remain repository workflows.

- Added strict saved JSON configuration, read-only `check` and configured `run`
  with frozen file hashes, completed upstream provenance and merged/per-sample
  count checks. Reference/map identities and recorded command settings must agree.
- Selected MIT with a separate citation request; shipped CITATION.cff, `cite`,
  a focused crate README/payload, exact-archive acceptance and a pending human
  reviewer exercise. Kept the ASCII banner in interactive startup only.
- Added installed native Rust `prepare-shared`, `sensitivity` and `example`
  commands, replacing repository/Python dependencies for the shared-graph
  workflow. Added binary-only acceptance tests, seven-scenario Python parity,
  executable/input/output checksums and an isolated installed-guide check.
- Added explicit shared-topology preparation with edge/support/coverage audits,
  zero-count selected edges, cycle eligibility and an annotated full-cohort
  statistics family. Added shared one-factor sensitivity with optional Rust
  statistics, failure provenance and an eight-sample synthetic walkthrough.
  The proposed biological thresholds still await team sign-off.
- Added optional `--storage-mode stream-bam` for single-sample and sequential
  batch preprocessing: supervised HISAT2/samtools pipes avoid SAM files while
  preserving junction-counting rules and Cleaver cleaning/indexing.
- Added disk estimates/reserve checks, bounded sort-memory options, partial
  output handling, per-process logs/exit codes and interruption propagation.
- Added pipeline-failure/signal regressions and actual-tool comparisons of
  retained-SAM/streaming-BAM counts, QC and disk-spilling sort behavior.

- Added a diagnostic junction-assignment/strand audit, a paired-input/storage
  preflight, an observed-candidate real-prefix benchmark, scientific pilot
  figures and an unpublished source-package preparation command.
- Added permutation/BH regressions and corrected the explanation of small
  sample p-value resolution without changing the study's frozen testing policy.
- Added methods/results drafts and a shared speaking plan; compiled the simple
  roadmap section with working CLI examples and no assigned team roles.

- Joined supplied gene symbols, risk/environment/overlap flags and categories
  directly into statistics, ranked results and the HTML report, retaining
  annotations for untested genes and leaving missing metadata blank.
- Added a checksummed paired-FASTQ prefix extractor and independent real-pilot
  BAM/SVD verifier. Reference preparation now supports lower-memory index
  builds, cleanup of its expanded FASTA copy after success, and rejection of
  partial or still-running indexes at alignment time.

- Added the runnable Cleaver/HISAT2 reference and sample adapter, SAM junction
  counting, sequential batch processing, QC and checksummed provenance.
- Added same-sample baselines, topology flags, optional deterministic ranking,
  preserved gene-family annotations, and descriptive support-grid sensitivity.
- Added `scripts/check_project.py`, paired-read external-tool smoke testing,
  sensitivity CI checks, and beginner workflow guides. Fixed a parallel-test
  temporary-directory collision.

- Added `stats` and menu option 7: explicit-family exact permutation tests,
  BH correction, per-gene statistics TSV, and a statistical section in the
  offline report. Test direction, independent design, sample minimum, family
  input/label, and output destination are required.
- Recorded untested genes and their reasons, repeated-subject rejection,
  assignment counts, family sizes, and gene-list checksums/settings.

- Added `report <junctions.tsv> <manifest.tsv> --output <folder>` and menu
  option 6: an offline academic report, sample-level comparison plots,
  exclusions, and links to each sample/gene's Hodge decomposition SVG.
- Kept report and comparison exports numerically identical and recorded all
  report files in the run manifest. Existing output folders remain protected.
- Added workflow checks for report/table agreement, missing groups, escaped
  input identifiers, figure links, and interactive report generation.

### Prototype baseline included in 0.1.0

The experimental crate contains the reusable Splice Girl Rust CLI.

- Added strict five-column junction TSV parsing with line-specific errors.
- Added deterministic directed splice graphs and Hodge gradient/cycle
  decomposition with reconstruction, kernel, orthogonality, and energy checks.
- Added `validate`, one-gene `gene`, all-gene `compare`, and `visualize` CLI
  workflows plus the interactive `--start` menu.
- Added case/control manifest matching, explicit exclusions, descriptive
  group summaries, and five-file reproducible result bundles.
- Added synthetic fixtures for linear, branching, exon-skipping, cyclic, noisy,
  disconnected, duplicate, self-loop, and zero-support inputs.
- Added beginner documentation, a source map, decisions log, release checklist,
  and automated Rust checks.

This version is an experimental software release, not a claim that the GSE315857
biological analysis is complete. Real preprocessing, final junction tables,
approved statistical inference, and biological interpretation remain separate
study milestones.
