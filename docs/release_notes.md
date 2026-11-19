# Splice Girl 0.1.0 — experimental crate release candidate

Prepared November 17, updated November 19, 2026. This is the experimental
0.1.0 crate candidate. MIT with a separate citation request was selected
November 18. Package acceptance and Cargo's publish dry-run pass; no tag or
release publication has been performed.

## Included software

- Rust CLI for validated junction tables, graph Hodge
  decomposition, descriptive comparisons, explicit-family exact tests,
  candidate-set testing, BH correction and optional rankings.
- Saved signed components, same-sample count/topology baselines, provenance,
  three-panel graph figures and plain offline reports.
- Installed native `prepare-shared`, `sensitivity` and `example` commands:
  shared topology, eligibility audits, one-factor grids, optional statistics
  and embedded synthetic inputs work without Python or the checkout.
- Saved JSON configuration, read-only `check`, configured `run`, frozen input
  identities and completed pipeline/reference/map provenance checks.
- MIT license, CITATION.cff and `cite`, a focused crate README/payload, actual
  64-file archive installation acceptance and a pending independent-reader
  exercise.

The Cargo package is the Rust junction-analysis CLI. The repository also has
Python/Cleaver/HISAT2 preprocessing, validation and plotting workflows; these
are not installed by `cargo install splice-girl` and are outside the crate
archive.

## Requirements

Rust 1.85 or newer for building the CLI. Locked `serde` and `serde_json`
dependencies handle configuration/provenance JSON. Python is needed only for
the repository preprocessing/helper scripts (Python 3.11+ for the release
packer); Cleaver, HISAT2 and compatible BAM tools are needed only for read
processing. NumPy is optional for independent SVD verification and Matplotlib
for scientific figures. The installed core CLI, reports and synthetic examples
can run offline after building.

## Current limits

The real-input validation is one deterministic control prefix. Full-cohort
preprocessing, comparability, eligibility, case/control results and biological
interpretation remain pending. The high exon-pair ambiguity rate and library
methods discrepancy require review. Independent evidence and beginner reviews,
the main Overleaf source and final scientific/release review are still outstanding.

The exact crate archive passed offline package tests and isolated installed-
guide acceptance; see the [November 19 archive record](software_validation.md#010-crate-archive-acceptance--november-19-2026).

## Prepare a local source archive

```sh
python3 scripts/prepare_release.py --output-dir target/release-review
```

This creates an unpublished source archive and SHA-256 manifest, recording
the base commit and uncommitted files. It does not change Git history or
upload anything. Large inputs, references, alignments and generated run
folders are excluded. Build and run the demo from the extracted snapshot
before treating it as a release candidate. Follow the
[release checklist](release_checklist.md) before any publication.
