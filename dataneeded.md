# Data needed for Splice Girl

This is the full data checklist. The raw FASTQ files, genome files, BAM files,
and STAR index should stay outside Git because they are very large. Keep their
file paths and checksums in the project notes.

## Needed for the first small test

- [ ] **Human genome FASTA** — the DNA sequence that the RNA reads will be
  matched against.
- [ ] **Matching GTF annotation** — the file that tells us where genes, exons,
  and transcripts are located.
- [ ] **STAR genome index** — the computer-ready version of the genome FASTA
  and GTF. It must use the same genome build as the FASTA and GTF.
- [ ] **Junction gene map** — made from the matching GTF with:
  `python3 preprocessing/build_junction_gene_map.py`.
- [ ] **B-01 read 1 and read 2** — the paired raw FASTQ files for one control
  sample.
- [ ] **DE-01 read 1 and read 2** — the paired raw FASTQ files for one
  deep-infiltrating-endometriosis sample.
- [ ] **Pilot sample sheet** —
  `data/metadata/E-MTAB-15117_pilot.tsv`, which connects each sample name to
  its condition, tissue, and FASTQ download links.

## Metadata needed for every sample

- [ ] A stable sample ID.
- [ ] The condition or group: control, deep-infiltrating endometriosis,
  ovarian endometriosis, or peritoneal endometriosis.
- [ ] The tissue or lesion site.
- [ ] The study, sample accession, and sequencing run accession.
- [ ] The individual ID, age, sex, and other available variables that could
  affect a comparison.
- [ ] The FASTQ file names and whether the sample is paired-end.
- [ ] The genome build, GTF release, and STAR index version used for that
  sample.

## Needed after the pilot works

- [ ] The remaining E-MTAB-15117 paired FASTQ files, if the group approves
  processing the full study.
- [ ] The final comparison groups chosen by the research team.
- [ ] Any batch, library, or sequencing-center information available for
  checking technical differences.
- [ ] Optional comparison data from GSE179640, if the group decides to use the
  tissue-matched endometrium dataset as a second analysis.

## Files Splice Girl must receive

For each sample, preprocessing must produce a tab-separated junction file with
these five columns:

```text
sample_id    gene_id    from    to    count
```

The `count` value is the number of reads supporting that splice junction. Each
row must remain connected to its sample and gene; samples must not be pooled
during preprocessing.

## Quality and record files

- [ ] Cleaver read-cleaning report.
- [ ] STAR alignment log and `SJ.out.tab` file.
- [ ] Sorted BAM file and BAM index for each sample.
- [ ] Junction-conversion quality report.
- [ ] Checksums or recorded download locations for raw files and references.
- [ ] A note stating how many reads were kept, aligned, and converted into
  junctions.

## What the final analysis uses

The final Hodge analysis uses the cleaned junction TSV files, not the raw FASTQ
files directly. It compares the graph measurements for the chosen sample
groups after the preprocessing and quality checks pass.
