# Splice Girl

An experimental research CLI for graph Hodge analysis of RNA splice-junction
counts. It prepares a common graph across samples, separates normalized edge
signal into gradient and cycle-space components, and exports comparisons,
explicitly configured permutation tests, SVG figures and offline reports.

## Install

Install from crates.io with Rust 1.85 or newer:

```sh
cargo install splice-girl --locked
```

Or build the current source checkout:

```sh
cargo install --path . --locked
```

No Python or source checkout is needed to run the installed junction-analysis
commands. Building resolves the JSON/configuration dependencies recorded in
Cargo.lock; analysis itself does not require network access.

## Try the example

```sh
splice-girl example --output demo-inputs
```

Open `demo-inputs/README.md` and follow its commands. The inputs are invented:
eight samples exercise measured zeros, missing coverage and cycle eligibility.
Existing output folders are protected.

Or use its saved configuration directly:

```sh
splice-girl check --config demo-inputs/analysis.json
splice-girl run --config demo-inputs/analysis.json
```

Preflight writes nothing. The run saves the configuration, checked file hashes,
eligibility decisions and `demo-inputs/configured-results/report.html`.
For biological pipeline inputs, the configuration guide explains the required
completed batch/reference/map provenance and consistency checks.

The installed tool includes `analyze`, `validate`, `gene`, `compare`, `stats`,
`report`, `visualize`, `prepare-shared`, `sensitivity`, `check`, `run`, `cite` and the interactive
`--start` menu. Use `splice-girl --help` for the command list.
For help with one command, use `splice-girl analyze --help` or
`splice-girl report --help` in the current source version.

## Inputs and interpretation

Junction TSV columns are exactly `sample_id`, `gene_id`, `from`, `to`, `count`.
Comparison workflows use a separate manifest with `sample_id` and `condition`
(`case` or `control`). Counts remain separate for each biological sample.
Raw FASTQ cleaning/alignment is a separate Cleaver/HISAT2 preprocessing workflow.

A cycle-space component is a mathematical projection relative to the supplied
graph. It does not establish circular RNA, observed use of a zero-count edge,
disease causation or a diagnostic result. Biological usefulness and full-cohort
validation are still under investigation. Independent researcher usability
review remains outstanding; automated acceptance checks are developer checks.

See the [project documentation](https://github.com/simransummermalik/TopologicalAnalysisofRNASplicinginEndometriosis/tree/main/docs)
for input contracts, equations, processing provenance, scope and limitations.
Report reproducible software problems through the
[issue tracker](https://github.com/simransummermalik/TopologicalAnalysisofRNASplicinginEndometriosis/issues).

## License and citation

Splice Girl 0.2.0 and later source releases are licensed under CC BY 4.0.
Commercial and noncommercial reuse are allowed with attribution and notice of
changes under the license terms. The published 0.1.0 crate remains MIT-licensed.

If this software contributes to a paper, please cite it and record the version
or commit used. Run `splice-girl cite` for the supplied citation metadata,
also shipped as `CITATION.cff`. The paper-citation request is separate from
the license's attribution requirement.
