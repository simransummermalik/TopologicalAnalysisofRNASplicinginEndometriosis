# Splice Girl release checklist

This checklist describes the smallest honest software release. It does not
claim that the biological study is finished.

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
4. Run the demo commands from the README in a fresh checkout.
5. Run `cargo install --path .` from a clean source copy and check
   `splice-girl --version`.
6. Confirm that no raw reads, references, BAM files, or generated run folders
   are staged.
7. Update the version, changelog, and software-validation date together.
8. Record the license decision before publishing a release artifact.

## Current scope boundary

The current release candidate supports validated synthetic analysis and
descriptive sample-preserving comparisons. It does not claim completed
splice-aware alignment, final real-data junction tables, approved inferential
statistics, multiple-testing correction, or biological interpretation.
