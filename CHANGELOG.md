# Changelog

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
