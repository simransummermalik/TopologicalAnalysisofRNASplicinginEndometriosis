# Splice Girl release checklist

This checklist describes the smallest honest software release. It does not
claim that the biological study is finished.

The [release notes](release_notes.md) and `scripts/prepare_release.py` prepare
an unpublished, allowlisted source archive with a checksum manifest. A license
choice and final review remain necessary before publication; no release tag
has been created.

## Included now

- Rust CLI source and dependency-free `Cargo.lock`.
- Synthetic fixtures and runnable demo inputs in `examples/demo/`.
- Parser, graph, Hodge, comparison, validation, export, and SVG tests.
- Beginner installation and command instructions in [`README.md`](../README.md).
- Input and output contracts in [`docs/data_contract.md`](data_contract.md).
- Mathematical limits and validation identities in [`docs/math_spec.md`](math_spec.md).
- A run manifest with input checksums for saved result bundles.

## Before tagging a version

1. Run `cargo fmt --all -- --check`.
2. Run `cargo test` and `cargo clippy --all-targets -- -D warnings`.
3. Run the Python preprocessing unit tests with
   `python3 -m unittest discover -s preprocessing/tests`.
   The combined command `python3 scripts/check_project.py` runs these plus
   gene-set and sensitivity tests, the Rust checks and benchmark. Use
   `--with-tools` for the actual tiny Cleaver/HISAT2 batch integration test.
4. Run the demo commands from the README in a fresh checkout.
5. Run `cargo install --path .` from a clean source copy and check
   `splice-girl --version`.
6. Confirm that no raw reads, references, BAM files, or generated run folders
   are staged.
7. Update the version, changelog, and software-validation date together.
8. Record the license decision before publishing a release artifact.

## Current scope boundary

The current release candidate supports validated synthetic analysis,
descriptive comparisons, and explicitly configured exact unpaired permutation
tests with BH correction and reporting. The real study still needs completed
splice-aware alignment, checked junction tables, approved eligibility/design
and testing-family rules, real baseline/sensitivity results, and biological
interpretation. Reference preparation, sequential Cleaver/HISAT2 processing,
baseline exports, optional deterministic rankings and support-grid software
are implemented; see [the run guide](running_the_tool.md).
