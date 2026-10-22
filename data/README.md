# Data directory

The current Splice Girl study is
[`GSE315857`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315857):
paired-end bulk RNA-seq from proliferative-phase eutopic endometrium.

- Four samples are from women with confirmed endometriosis.
- Four control samples are from women who were laparoscopically negative for
  endometriosis.
- All eight samples use the same tissue and menstrual phase.

The current sample table is
[`metadata/samples.tsv`](metadata/samples.tsv), with a plain-language guide in
[`metadata/README.md`](metadata/README.md). It is the authoritative manifest
for the GSE315857 analysis.

## What is stored here

- `metadata/samples.tsv`: one row per current sample, including accessions,
  group, tissue, file locations, processing status, and eligibility.
- `galaxy_data/october_20/`: raw FastQC and MultiQC reports for all eight
  samples plus the tracked results of the one-sample Cleaver pilot.
- `raw/`: local multi-gigabyte FASTQ files. These are ignored by Git.
- `geo/` and the other files in `metadata/`: records of datasets considered
  earlier. They are retained for provenance and are not current GSE315857
  inputs.

Raw and cleaned FASTQ files are too large for the repository. Their locations
are recorded in the manifest; only small metadata and quality-control reports
belong in Git. Junction files have not been generated yet because the matching
project genome and annotation must be finalized first.
