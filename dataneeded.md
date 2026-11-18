# Data needed for Splice Girl

This checklist follows the current GSE315857 analysis. Large FASTQ, reference,
index, SAM, and BAM files stay outside Git. Their identifiers, paths, versions,
and checksums belong in the project records.

## Sample information

- [x] Eight stable SRA-run sample IDs.
- [x] Four `control` and four `case` labels.
- [x] GEO sample, SRA experiment, and SRA run accessions.
- [x] Tissue and menstrual phase: proliferative-phase eutopic endometrium for
  every sample.
- [x] Paired-end library layout and the current read locations.
- [x] Separate internal subject labels for all eight biological subjects.
- [x] Public-metadata gaps, including unavailable patient IDs and batch, are
  written explicitly instead of guessed.
- [x] Final project genome build: GENCODE v50 GRCh38.p14 primary assembly.
- [x] Matching GENCODE v50 comprehensive primary-assembly GTF downloaded,
  checksum-verified, and checked against the FASTA sequence names.
- [x] Final junction coordinate, exon-node, and strand conventions recorded in
  [`docs/data_methods.md`](docs/data_methods.md#frozen-junction-and-graph-conventions).

The complete table is [`data/metadata/samples.tsv`](data/metadata/samples.tsv).
Per-subject ages and original clinical IDs are not available in the public
sample records; the study reports an overall age range of 18–45 years.

## Raw reads and quality control

- [x] All eight paired-end SRA runs were retrieved in Galaxy.
- [x] Raw FastQC and MultiQC reports exist for all eight samples.
- [x] One control pilot, `SRR36710387`, has a recorded Cleaver cleaning run,
  local cleaned paired reads, checksums, and post-cleaning FastQC/MultiQC.
  Its original raw files are not currently in `data/raw/`.
- [ ] Run the same reviewed Cleaver cleaning on the remaining seven samples.
- [ ] Record cleaned-file checksums for the remaining seven samples.
- [ ] Review the complete post-cleaning cohort together before alignment.

The raw and cleaned reads are intentionally ignored by Git. The small reports
are under [`data/galaxy_data/october_20`](data/galaxy_data/october_20).

## Human reference files

- [x] Human genome FASTA: checksum-verified GENCODE v50
  `GRCh38.primary_assembly.genome.fa.gz`.
- [x] Matching GENCODE v50 comprehensive primary-assembly GTF annotation.
- [x] Exact genome and annotation versions recorded in the manifest; source
  URLs, sizes, and checksums recorded in
  [`docs/data_methods.md`](docs/data_methods.md).
- [x] HISAT2 index built from that same genome (November 17). HISAT2 is the
  splice-aware Bowtie-family aligner; plain Bowtie2 should not be used for
  intron-spanning RNA-seq reads.
- [x] FASTA and GTF source URLs, sizes, checksums, and integrity checks.
- [x] HISAT2 2.2.3 and the exact executed index-building command, source/output
  hashes and splice-site map recorded in the [pilot bundle](data/pilots/SRR36710387_prefix1m_nov17/README.md).

## Alignment and splice junctions

**Completed bounded pilot:** the first one million cleaned pairs from the
local control have a full-reference alignment, BAM/index, 32,575 junction rows
and an independently verified CDH11 result. See the [review](data/pilots/SRR36710387_prefix1m_nov17/README.md).
The boxes below still refer to the complete sample/cohort, not its prefix.
Review the large annotation-ambiguity exclusion rate before expanding.

- [ ] Align the cleaned paired reads for the pilot sample with HISAT2.
- [ ] Save the alignment summary, SAM/BAM output, sorted BAM, and BAM index.
- [ ] Derive splice-junction counts without losing the sample ID.
- [ ] Map junction coordinates to unique gene/exon pairs using the matching
  GTF and record all unannotated or ambiguous exclusions.
- [ ] Repeat the approved alignment and junction process for the remaining
  eligible samples.
- [ ] Run junction QC and record exclusion reasons.

For every retained sample, preprocessing must produce a tab-separated junction
file with these five columns:

```text
sample_id    gene_id    from    to    count
```

`sample_id` must match the SRA-run value in the sample manifest. `count` is the
number of reads supporting that directed splice junction. Samples must stay
separate throughout preprocessing.

## What the final analysis uses

Splice Girl analyzes the cleaned junction TSVs, not FASTQ files directly. The
Rust tool will build one graph per sample and gene, run the Hodge decomposition,
and compare the resulting measurements only after preprocessing, alignment,
junction QC, and eligibility checks are complete.
