# Streaming-BAM validation — November 18, 2026

**PASS.** The full check passes **200 tests**: 100 Rust, 39 preprocessing,
56 gene-evidence and five script tests. Formatting, strict Clippy and the
release build also pass. See [full_checks.log](full_checks.log) and the
source/log hashes in [verification.json](verification.json).

## What was actually run

```sh
python3 scripts/check_project.py --with-tools
```

This ran the full project checks and actual Cleaver 1.0.1, HISAT2 2.2.3 and
samtools 1.25 on tiny synthetic paired reads against the tiny reference.
Both samples were processed in both storage modes and passed through the
Rust report command. Stream mode produced final indexed BAMs, junction
tables and complete manifests, with no saved SAM or leftover `.partial`
BAM. Its junction tables and every QC counter matched retained-SAM mode.

A second external-tool test used **12,000 synthetic pairs per sample** and
a 1 MiB sort-memory setting to force actual disk spill:

```sh
python3 -c 'import sys; sys.path.insert(0, "scripts"); import check_project; check_project.smoke(12000)'
```

The [spill-test log](sort_spill_smoke.log) confirms temporary sort files were
used for both CASE and CONTROL. Each sample's single synthetic junction had
24,000 supporting alignments: both mates count under the unchanged read-level
policy. Counts and all QC counters agreed across modes. The 1 MiB setting is
deliberately for this test, not a recommended setting for the human cohort.
All reference/read/alignment outputs for these tests were temporary synthetic
files, not the project's saved human inputs.

## Failure and interruption checks

The added runtime tests exercise:

- Producer failure even when the pipe's consumer succeeds.
- Consumer failure while a producer is still writing, without hanging.
- Keyboard interruption and an actual SIGTERM sent to an isolated test parent;
  the parent records interruption and stops both streaming child processes.
- Low space before any tools/run folder are started, and a reserve reached
  between stages.
- A batch whose total estimate exceeds the available space even though one
  sample alone would fit.
- Failures during alignment, BAM decoding/counting and indexing: partial files
  remain, with no successful final bundle published.
- A batch interrupted after one sample: completed work is retained and the
  cohort table is not published.
- Exact equivalence of file and stdin junction conversion, including QC.

Orchestration tests use artificial format bytes where only failure propagation
is under test; the actual-tool checks above verify real BAM handling. A separate
real CLI refusal deliberately requested an impossible free-space reserve;
it exited 1 before tools started and created no output folder. See
[disk_refusal.log](disk_refusal.log).

## Limits and use

Use `--storage-mode stream-bam` on an existing single-sample or batch command.
The [usage guide](../../HISAT2_GUIDE.md#save-disk-space-with-streaming-bam)
explains both modes, sort memory, disk estimates, logs, partial names and retries.
The original retained-SAM mode remains the default. A plan is read-only.

This validates implementation on synthetic reads. It does not change the
frozen annotation/counting rules, approve real sample eligibility, establish
a biological result or guarantee that the complete cohort fits this laptop.
No personal files or system settings were changed. No commit or push performed.
