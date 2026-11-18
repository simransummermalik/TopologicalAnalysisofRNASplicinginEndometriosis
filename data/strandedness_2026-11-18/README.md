# Library strandedness check: all eight GSE315857 runs

November 18, 2026. Result: **every run is RF** (read 1 antisense, read 2 sense
to the transcript), the orientation of dUTP-style stranded libraries.

| Sample | Group | Read 1 antisense | Read 2 sense | Reads matched (read 1 / read 2) |
|---|---|---:|---:|---:|
| SRR36710390 | control | 99.11% | 99.10% | 5,376 / 5,363 |
| SRR36710389 | control | 98.87% | 98.83% | 6,120 / 6,043 |
| SRR36710388 | control | 99.15% | 99.14% | 6,121 / 6,131 |
| SRR36710387 | control | 99.33% | 99.05% | 5,252 / 5,239 |
| SRR36710386 | case | 99.49% | 99.08% | 5,836 / 5,893 |
| SRR36710385 | case | 98.73% | 98.54% | 4,820 / 4,849 |
| SRR36710384 | case | 98.14% | 97.95% | 4,520 / 4,479 |
| SRR36710383 | case | 98.47% | 98.86% | 5,937 / 5,955 |

Full counts are in [`results.tsv`](results.tsv).

## Method

- **Reads:** the first 200,000 reads of each mate of each run, streamed from
  the ENA mirror of the SRA runs
  (`http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR367/0NN/<run>/<run>_<mate>.fastq.gz`).
  These are raw reads before cleaning. Only these few megabytes were read; no
  full file was downloaded or kept.
- **Reference:** nine RefSeq mRNAs of highly expressed genes, fetched from NCBI
  on November 18, 2026: ACTB `NM_001101.5`, GAPDH `NM_002046.7`, EEF1A1
  `NM_001402.6`, RPLP0 `NM_001002.4`, TPT1 `NM_003295.4`, FTL `NM_000146.4`,
  VIM `NM_003380.5`, B2M `NM_004048.4`, EEF2 `NM_001961.4`. Saved as
  [`reference_mrnas.fa`](reference_mrnas.fa), SHA-256
  `d284d63f9c2d33cf9b464e712ee4e12bb7ede269fce9a8142a1c1356df83501e`.
- **Classification:** [`preprocessing/strand_kmer_check.py`](../../preprocessing/strand_kmer_check.py)
  compares 25-mers from each read with both strands of these mRNAs. It needs
  no aligner or genome index.

```sh
python3 preprocessing/strand_kmer_check.py \
  --mrna data/strandedness_2026-11-18/reference_mrnas.fa \
  --sample SRR36710390 SRR36710390_1.fq SRR36710390_2.fq \
  ...   # one --sample per run
```

## How it fits with other evidence

- The [pilot alignment audit](../pilots/SRR36710387_prefix1m_nov17/assignment_audit/README.md)
  found 99.72% RF orientation for SRR36710387 from HISAT2 spliced alignments.
  This check gives 99.33% for the same run by an independent method.
- The GEO records name the NEBNext rRNA Depletion Kit v2 but not the
  library-preparation kit; the associated paper describes a KAPA mRNA kit.
  Neither source states the orientation. The reads do, consistently for all
  eight runs and both groups.

The first reads of a file are not a random sample; the check is about
library chemistry, which is the same for every read in a run.
