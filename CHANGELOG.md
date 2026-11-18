# Changelog

## Unreleased

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

## 0.1.0 — release candidate

The first release candidate contains the reusable Splice Girl Rust prototype.

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

This version is a software release candidate, not a claim that the GSE315857
biological analysis is complete. Real preprocessing, final junction tables,
approved statistical inference, and biological interpretation remain separate
study milestones.
