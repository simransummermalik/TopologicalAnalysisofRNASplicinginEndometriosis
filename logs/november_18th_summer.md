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
