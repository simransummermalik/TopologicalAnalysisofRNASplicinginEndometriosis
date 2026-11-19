# Splice Girl release checklist

This checklist describes the smallest honest software release. It does not
claim that the biological study is finished.

The [release notes](release_notes.md) describe the 0.1.0 experimental crate.
MIT and a separate citation request were selected November 18. On November 19,
the exact 64-file crate archive, extracted tests/install walkthrough, and
`cargo publish --dry-run --locked` passed. No upload or release tag has been
performed. Independent reader feedback and full biological validation remain
open and are clearly identified as limits of this experimental release.

## Included now

- Rust CLI source and locked JSON-parsing dependencies in `Cargo.lock`.
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

MIT is recorded in Cargo.toml and LICENSE; CITATION.cff provides the requested
software citation. Citation in papers is requested rather than imposed as an
additional MIT condition.

## Check the exact crate payload

```sh
python3 scripts/check_crate.py --output-dir target/crate-review
```

This packages without uploading, checks the archive against the software-only
allowlist, extracts it, runs its Rust tests and installs that extracted source.
It then runs all seven installed-guide commands outside the checkout with
PATH empty and retains logs/artifacts. Research data, PDFs, class logs, website
files, preprocessing scripts and generated runs are outside the crate payload.
The stable CI job repeats this check; no publishing credentials are used.

The installed shared-graph workflow has its own reproducible developer check:

```sh
python3 scripts/check_installed_cli.py --output-dir target/installed-review
```

It installs into a private review folder, executes every command in the
embedded guide from a temporary folder outside the checkout with PATH empty,
and retains reports, checksums and logs. This verifies that installed
preparation/sensitivity do not depend on Python or repository scripts at
runtime. The exact-crate check above covers distribution payloads; independent
human feedback and scientific validation remain separate.

## Current scope boundary

The current release candidate supports validated synthetic analysis,
descriptive comparisons, and explicitly configured exact unpaired permutation
tests with BH correction and reporting. The real study still needs completed
splice-aware alignment, checked junction tables, approved eligibility/design
and testing-family rules, real baseline/sensitivity results, and biological
interpretation. Reference preparation, sequential Cleaver/HISAT2 processing,
baseline exports, optional deterministic rankings and support-grid software
are implemented; see [the run guide](running_the_tool.md).
