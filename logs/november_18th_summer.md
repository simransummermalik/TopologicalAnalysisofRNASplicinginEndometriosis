# November 18, 2026 — Summer

**Prepared locally; no commit or push performed by the assistant.**

I added an optional storage-saving mode to preprocessing so HISAT2 can send
alignments straight into sorted BAM instead of writing a huge SAM file first.
Cleaver still handles cleaning and indexing. The runner now checks estimated
disk space before starting, records each stage's logs and exit codes, and
stops its child processes if a tool fails or the run is interrupted. Partial
files stay labeled until the streaming run succeeds, and the batch will not
publish a complete cohort result after a failed sample. I tested the new mode
against the original one using synthetic reads: the counts and QC matched
exactly. A larger synthetic test also forced the sorter to use temporary disk
files. The complete check passes 200 Rust/Python tests, and I saved the test
logs, source hashes and usage instructions. These tests did
not process our human read files or change personal files or laptop settings;
the full eight-sample analysis still needs the remaining inputs, reviewed
settings and sufficient storage.

I also built the software for Vania's shared-graph proposal. It takes complete
counts from the cohort, chooses the same edges for each sample of a gene,
keeps real zero counts separate from samples with too little coverage, and
records why each sample or gene was left out. It only passes fully covered
genes with a shared graph cycle into the statistics family. I added a small
eight-sample made-up example and a runner for the seven proposed sensitivity
settings, including optional Rust statistics and reports. The example and new
tests pass, and the full project check now passes 216 tests. I also fixed the
mocked batch tests so low laptop disk space does not stop them before they
reach the behavior they test; the real disk protection stays in place. The
group still needs to approve the actual study cutoffs before applying them to
the human cohort. No commit or push was performed by the assistant.

I moved shared-graph preparation and sensitivity into the Rust command itself,
so someone who installs Splice Girl can use them without finding our Python
scripts or staying in the project folder. There are now `prepare-shared`,
`sensitivity` and `example` commands. The example writes the small made-up
inputs and a step-by-step guide, and sensitivity can run the existing Rust
statistics and reports directly. I installed it into a separate review folder
and ran the guide outside the repository with PATH empty; everything passed.
I checked that all seven scenarios match the earlier Python results, including
exclusions and statistics. The complete check now passes 226 tests, and the
Rust tests and Clippy also pass on Rust 1.85. The usual installed command was
not replaced, and I did not commit, push or publish anything.
