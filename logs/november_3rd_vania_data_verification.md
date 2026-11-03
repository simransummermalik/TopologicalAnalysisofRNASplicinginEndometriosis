# November 3, 2026 — Vania (data verification)

**Status:** Prepared for push; not yet committed when this log was created

## Contribution

I ran a first verification pass over the GSE315857 pilot sample's
existing data-collection outputs, per `dataneeded.md`, ahead of helping with
data collection and verification more broadly. This is a read-only check —
no files were changed, no new data was collected.

### What was checked, and how, against the pilot sample (`SRR36710387`)

- **Sample identity.** Confirmed the post-cleaning MultiQC data folder
  (`data/galaxy_data/october_20/cleaver_pilot/post_cleaning_multiqc/`)
  genuinely references `SRR36710387` internally, not just in its folder
  name, ruling out an accidental copy from another run.
- **Cleaning direction.** Read `SRR36710387.fastp.json` directly: Q30 rate
  improved from `0.884059` before filtering to `0.890305` after, and Q20
  from `0.953679` to `0.958907`. This independently supports the manifest's
  `cleaver_status: complete_postclean_qc_pass` claim rather than trusting
  the label alone.
- **Internal arithmetic.** `fastp.json`'s `filtering_result` categories
  (`passed_filter_reads`, `low_quality_reads`, `too_short_reads`, etc.) sum
  to exactly `98,954,118`, matching `before_filtering.total_reads` exactly.
- **Cross-document consistency.** `docs/data_methods.md`'s GTF row counts
  ("11,254,876 total; 5,087,789 exon rows") are repeated identically in
  `logs/october_22nd_summer.md` — no drift between the log that recorded the
  original work and the doc meant to be the frozen reference.
- **Report scope.** Confirmed `controls_multiqc_html.html` (cited as the
  pilot's `raw_qc_report`) legitimately covers `SRR36710387` alongside
  `SRR36710388` and `SRR36710389` — it is a combined 4-control report, not a
  mismatched per-sample path.

### A real limit found, not a gap in the check

Could not re-verify the recorded FASTA/GTF SHA-256 checksums or the pilot's
raw/cleaned FASTQ files, because `data/raw/`, `preprocessing/runs/`, and
`references/gencode_v50/` do not exist on this checkout at all — not hidden
by `.gitignore`, genuinely absent. Those large files only exist on whichever
machine actually downloaded and processed them. File-level checksum
re-verification needs to happen there, not from this environment.

### One open question, raised rather than resolved

`fastp.json` reports `after_filtering.total_reads = 98,953,584` but
`filtering_result.passed_filter_reads = 98,953,851` — a 267-read gap between
two numbers that sound like they should match. This may be a known fastp
behavior (e.g. a secondary trimming step not counted in `filtering_result`)
rather than a problem, but it wasn't confirmed either way and is worth
asking whoever ran Cleaver rather than guessing.

## Current scope

This is a verification note only. It does not clean, align, or collect any
new sample data, and does not change the status of any `dataneeded.md`
checklist item. The next unchecked items remain: Cleaver cleaning and
checksums for the remaining seven samples, the HISAT2 index, and all of
alignment/junction derivation.
