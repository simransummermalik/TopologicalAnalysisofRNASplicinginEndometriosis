# Splice Girl crate releases

## 0.2.2 — patch release, published November 24, 2026

Adds focused command help, accepts checked cohort bundles assembled from
separate sample jobs, and refuses recorded inconclusive strand checks during
configured analysis. Installation instructions now lead with crates.io.
The license remains CC BY 4.0. Repository preprocessing also has stricter
strand-check handling; Python preprocessing scripts are outside the crate.

The exact archive is reviewed with:

```sh
python3 scripts/check_crate.py --output-dir target/crate-review-0.2.2
```

The 64-file archive passed its Rust tests and the isolated installed-tool
walkthrough. Cargo confirmed `Published splice-girl v0.2.2` on crates.io.
Install or update with `cargo install splice-girl --version 0.2.2 --locked`.
Start the menu with `splice-girl --start`.

## 0.2.1 — patch release, published November 19, 2026

This CC BY 4.0 patch release fixes the configured-analysis report footer so it
uses the package's declared license instead of stale MIT wording. The 64-file
crate archive passed Cargo's dry-run verification and the
[GitHub checks](https://github.com/simransummermalik/TopologicalAnalysisofRNASplicinginEndometriosis/actions/runs/44615757528),
including the extracted crate acceptance workflow. Install it with
`cargo install splice-girl --version 0.2.1`.

## 0.2.0 — license update, published November 19, 2026

The source release changed the license from MIT to Creative Commons
Attribution 4.0 International (CC BY 4.0), with approval from all four named
contributors. Version 0.1.0 remains under its original MIT license.

## Historical: 0.1.0 experimental release

Prepared November 17, updated November 19, 2026. This is the experimental
0.1.0 crate release, published to crates.io on November 19. MIT with a
separate citation request was selected November 18. Package acceptance and
Cargo's publish dry-run passed before publication. The release commit is on
GitHub `main`; no release tag was created.

This section records the original 0.1.0 release and its contents.

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
interpretation remain pending. Residual junction-assignment and coverage
questions require review. Independent evidence and beginner reviews, the main
Overleaf source and final scientific review are still outstanding.

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
