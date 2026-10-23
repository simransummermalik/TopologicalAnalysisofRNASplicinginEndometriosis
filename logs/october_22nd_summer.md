# October 22, 2026 — Summer

**Status:** R13 and D4 pushed; reference work prepared for push

**Commits:** `a9ed39e` — `adding output.rs`; `ef9a204` — `created the full GSE315857`

## Contribution

Today I completed Rust task **R13**, which makes Splice Girl save its analysis
instead of leaving the results only in terminal text. The existing commands
still work the same way, but `analyze` and `gene` now accept an optional
`--output <folder>` argument. A saved analysis includes sample-level Hodge
measurements, every junction's normalized, gradient, and cycle components, and
a JSON run manifest. A `gene` run also saves the included and excluded samples
in a comparison table.

The run manifest records the exact command, Splice Girl version, input paths,
SHA-256 input checksums, normalization and solver settings, numerical
tolerance, output filenames, exclusions, and failed mathematical checks. Each
TSV row carries the same run ID as the manifest. The program refuses to use an
existing output folder so an old result is not silently overwritten.

I also updated the README, demo guide, and data contract with the new command
and file descriptions. Finally, I audited `delegations_and_tasks/tasks.md`
against the files currently in the repository. The status board now marks the
completed setup, math, Rust, usability, pipeline-figure, and demo tasks and
names the actual blocker for unfinished real-data tasks.

## Validation

- `cargo test` passes: **35 tests passed, 0 failed**.
- `cargo clippy --all-targets -- -D warnings` passes.
- The `analyze --output` command produced valid sample and junction TSVs plus
  a JSON manifest from the directed-cycle fixture.
- The `gene --output` command produced all four result files from the demo and
  correctly recorded four included and two excluded samples.
- Python's JSON and TSV readers successfully parsed every generated file.
- Known SHA-256 test vectors for empty input and `abc` pass.

The completed R13 work and task-board audit are recorded in commit `a9ed39e`.

## D4 sample manifest

I also made the final working sample table for GSE315857 at
`data/metadata/samples.tsv`. It puts all eight samples in one place and records
which four are controls and which four are endometriosis cases, their GEO and
SRA accessions, their matched proliferative eutopic-endometrium tissue, where
their reads or reports are stored, and how far each sample has moved through
quality control. All eight have raw FastQC/MultiQC results. SRR36710387 is
marked as the completed Cleaver pilot, while the other seven are still waiting
for Cleaver. The genome, annotation, alignment, and junction fields are clearly
marked as pending instead of being guessed. I checked that the sample IDs and
accessions are unique and that the table has four cases and four controls.

## Exact human genome FASTA

I selected and downloaded the current GENCODE v50 GRCh38.p14 primary-assembly
human genome FASTA. The compressed file is stored locally at
`references/gencode_v50/GRCh38.primary_assembly.genome.fa.gz` and remains
ignored by Git because it is 845,635,028 bytes. Its MD5 exactly matches the
official GENCODE value, its gzip integrity check passes, and I recorded both
MD5 and SHA-256 checksums in `docs/data_methods.md`. I also selected the exact
matching GENCODE v50 primary-assembly GTF, but correctly left it marked as
pending because it has not yet been downloaded or checked. The sample manifest
now records the exact genome build for all eight samples. D5 is in progress
rather than blocked.
