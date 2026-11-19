# Splice Girl 0.1.0 — release review draft

Prepared November 17, updated November 18, 2026. This is an unpublished review
candidate. MIT with a separate citation request was selected November 18.
No tag or release publication has been performed.

## Included software

- Rust CLI for validated junction tables, graph Hodge
  decomposition, descriptive comparisons, explicit-family exact tests,
  candidate-set testing, BH correction and optional rankings.
- Saved signed components, same-sample count/topology baselines, provenance,
  three-panel graph figures and plain offline reports.
- Sequential Cleaver/HISAT2 adapter, checked exon-pair assignment, independent
  prefix/recount validation and a diagnostic ambiguity/strand audit.
- Explicit support/coverage sensitivity runner, synthetic acceptance examples,
  real-prefix benchmark, storage preflight and scientific plotting helper.
- Installed native `prepare-shared`, `sensitivity` and `example` commands:
  shared topology, eligibility audits, one-factor grids, optional statistics
  and embedded synthetic inputs work without Python or the checkout.
- Saved JSON configuration, read-only `check`, configured `run`, frozen input
  identities and completed pipeline/reference/map provenance checks.
- MIT license, CITATION.cff and `cite`, a focused crate README/payload, actual
  archive installation acceptance and a pending independent-reader exercise.

## Requirements

Rust 1.85 or newer for building the CLI. Serde/serde_json handle strict JSON
configuration/provenance parsing; numerical Hodge code remains unchanged.
Dependencies are locked and downloaded for the initial build. Python is
needed for preprocessing/helper scripts (Python 3.11+ for the release packer).
Cleaver, HISAT2 and compatible BAM
tools are needed for processing reads. NumPy is optional for independent SVD
verification; Matplotlib is optional for scientific figures. The core CLI,
saved reports and shipped synthetic examples can run offline after building.

## Current limits

The real-input validation is one deterministic control prefix. Full-cohort
preprocessing, comparability, eligibility, case/control results and biological
interpretation remain pending. The high exon-pair ambiguity rate and library
methods discrepancy require review. Independent evidence and beginner reviews,
the main Overleaf source and final scientific/release review are still outstanding.

The local source snapshot passed the full offline check and isolated install/
demo commands. See the [acceptance record](software_validation.md#isolated-source-package-acceptance--november-17).

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
