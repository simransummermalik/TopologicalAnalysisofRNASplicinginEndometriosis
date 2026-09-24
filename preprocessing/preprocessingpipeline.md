# RNA-seq preprocessing pipeline notes

These are planning notes for the step before Splice Girl. The goal is to turn
raw sequencing reads into a small, traceable junction table that follows the
[Splice Girl input contract](../docs/data_contract.md).

This is not implemented yet. The exact tools and settings should be confirmed
with Dr. White and the data team before anyone processes the full dataset.

## Proposed flow

```text
raw FASTQ reads
      |
      v
quality check (FastQC or confirmed equivalent)
      |
      v
trim and clean reads (tool to be confirmed)
      |
      v
remove or filter PhiX/control reads when needed
      |
      v
align cleaned reads to an indexed human reference genome
      |
      v
SAM alignment
      |
      v
sorted and indexed BAM (samtools or confirmed equivalent)
      |
      v
extract gene/exon/splice-junction features
      |
      v
quality-controlled junction TSV for Splice Girl
```

## What each step means

1. **Sample manifest:** Give every sample a stable `sample_id` before touching
   the files. Record whether reads are paired-end, the tissue, condition,
   accession, and original filenames.
2. **Reference preparation:** Record the human genome build and annotation
   release. Build or obtain the matching searchable genome index.
3. **Raw-read quality check:** Inspect read quality, length, adapter content,
   and obvious contamination before trimming.
4. **Trim and clean:** Remove adapters and low-quality bases with the tool that
   the data team approves. Save its version and settings.
5. **PhiX/control filtering:** Decide whether PhiX reads are present and how
   they will be removed. Do not assume that a trimming tool automatically does
   PhiX filtering.
6. **Splice-aware alignment:** Align cleaned reads to the matching genome index
   with the confirmed RNA-seq aligner. Record whether the settings preserve
   the junction information needed by this project.
7. **SAM to BAM:** Convert alignment output to sorted, indexed BAM files. Keep
   the sample ID attached to every file.
8. **Feature extraction:** Derive the gene, exon, or splice-junction features
   needed for the graph input. Record the extraction tool, annotation, and
   filtering settings.
9. **Final QC:** Check missing fields, low support, duplicate junctions,
   sample joins, and whether back-splice detection is supported.
10. **Splice Girl export:** Write rows with exactly these columns:

    ```text
    sample_id    gene_id    from    to    count
    ```

    The separators in the real file must be tab characters.

## Files to keep

For each sample, keep a small manifest that points to:

- original FASTQ files and their checksums;
- raw and cleaned quality reports;
- cleaned FASTQ files or their documented storage location;
- the reference build, annotation, and index version;
- sorted/indexed BAM files or their documented storage location;
- extracted junction features;
- the final standardized junction TSV;
- tool versions, commands, dates, and filtering settings.

Do not pool samples during preprocessing. A row must remain traceable to one
sample, one gene, and one directed junction. A missing measurement must remain
different from a measured zero.

## Decisions still needed

- Which public dataset and tissue design will be used?
- Which human genome build and annotation release match that dataset?
- Which trimming/cleaning tool will be used?
- Which tool will handle PhiX filtering, if it is needed?
- Which splice-aware aligner will be used: STAR, HISAT2, or another approved
  tool?
- Which feature extractor will produce junction counts?
- Can the chosen workflow detect back-splice junctions?
- Which support and coverage filters will be applied?

These choices belong in the data methods record before full-scale processing.
The current Rust CLI starts after this pipeline, when it receives the
standardized junction TSV.
