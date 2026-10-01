# RNA-seq preprocessing pipeline notes

These notes describe the work before Splice Girl's graph and Hodge analysis.
The goal is to turn raw sequencing reads into a small, traceable junction table
that follows the [Splice Girl input contract](../docs/data_contract.md).

The initial one-sample runner is documented in
[`preprocessing/README.md`](README.md). It uses the existing RAW Lab
[Cleaver-X](https://github.com/raw-lab/cleaver-x) command-line tool for read
cleaning and BAM indexing. STAR performs splice-aware alignment, optional
Bowtie2 filtering removes PhiX when appropriate, and project-specific helpers
create the final junction table.

## Current flow

```text
raw FASTQ reads
      |
      v
trim, filter, and write JSON QC (Cleaver-X: cleaver fastp)
      |
      v
remove PhiX/control reads when the experiment requires it (Bowtie2)
      |
      v
align cleaned reads to the matching human reference (STAR)
      |
      v
coordinate-sorted BAM and splice-junction table (STAR)
      |
      v
build BAM index (Cleaver-X: cleaver samtools index)
      |
      v
annotate and filter STAR splice junctions (project helper)
      |
      v
quality-controlled junction TSV for Splice Girl
```

## What each step means

1. **Sample record:** Give every sample a stable `sample_id`. Record whether
   reads are paired-end, the tissue, condition, accession, and original files.
2. **Reference preparation:** Record the human genome build and annotation
   release. Build or obtain the matching STAR genome index.
3. **Read cleaning:** Use `cleaver fastp` to trim and filter FASTQ records. Keep
   its JSON report and exact command. Paired-end overlap detection is enabled
   by the current runner; final thresholds still need dataset-specific review.
4. **PhiX filtering:** Decide from the library method and QC evidence whether
   PhiX removal is needed. Cleaver-X read cleaning does not itself remove PhiX;
   the optional workflow step uses Bowtie2 and a matching PhiX index.
5. **Splice-aware alignment:** Align the cleaned reads with STAR. Inspect the
   STAR logs, especially the input-read and mapping counts, before continuing.
6. **BAM indexing:** STAR creates a coordinate-sorted BAM. The workflow runs
   `cleaver samtools index` so tools can access it by genomic position.
7. **Feature extraction:** Convert `SJ.out.tab` coordinates and read support to
   gene-labelled directed junctions using an annotation-derived gene map.
8. **Final QC:** Check unannotated and ambiguous junction counts, support
   thresholds, duplicate rows, sample labels, and missing values.
9. **Splice Girl export:** Write exactly these tab-separated columns:

   ```text
   sample_id    gene_id    from    to    count
   ```

## Files to keep

For each sample, retain or record the location of:

- original FASTQ files and checksums;
- Cleaver-X JSON QC and cleaned FASTQ files;
- reference build, annotation release, and index version;
- STAR logs, sorted BAM, and BAM index;
- STAR splice-junction output;
- standardized Splice Girl junction TSV and conversion QC;
- tool versions, commands, date, and filtering settings.

Do not pool samples during preprocessing. Each row must remain traceable to one
sample, one gene, and one directed junction. Keep missing measurements distinct
from measured zeros.

## Decisions still needed

- Which public dataset and tissue design will be used?
- Which human genome build and annotation release match that dataset?
- Is PhiX filtering appropriate for these libraries?
- Which read-cleaning and minimum junction-support thresholds are justified?
- Does the final question require back-splice detection, which STAR's standard
  `SJ.out.tab` export does not represent as ordinary forward splice junctions?
- Which aligner setup will be used if STAR reports zero input reads on Apple
  Silicon?

Record these choices before processing the full dataset. The current Rust CLI
starts after this workflow, when it receives the standardized junction TSV.
