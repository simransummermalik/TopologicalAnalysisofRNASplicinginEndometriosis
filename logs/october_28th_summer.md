# October 28, 2026 — Summer

**Status:** Prepared for push

Today I made the Rust tool much more complete instead of only adding visual
polish. I added a `validate` command that runs the real parser, graph builder,
Hodge decomposition, and mathematical checks, then prints a short audit of the
rows, samples, genes, read support, and any mismatch between the junction file
and sample manifest. This lets us catch a broken input before starting a larger
analysis.

I also added an all-gene `compare` command. It keeps every biological sample
separate and reports case/control sample counts, means, medians, ranges, the
case-minus-control mean difference, exclusions, and a clear status for every
observed gene. It does not make up p-values or rankings before our statistics
rules are approved. The command can save five traceable files, including a new
`gene_summaries.tsv`, detailed sample decisions, junction components, sample
metrics, and the run manifest with checksums and settings.

To prove the all-gene workflow actually works, I added a two-gene synthetic
demo with opposite effect directions and command-level tests. I expanded the
interactive menu, help text, README, demo guide, data contract, priorities,
task board, and SQL-style project map. The project now has 56 passing tests.
The full suite passes on both Rust 1.85 and the current Rust compiler, strict
Clippy passes, an isolated `cargo install --path .` works, and the installed
program successfully validates, compares, exports, and makes an SVG. I also
recorded the install check and small-demo timing and memory results in
`docs/software_validation.md`.

