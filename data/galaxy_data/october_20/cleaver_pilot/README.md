# Cleaver pilot: SRR36710387

This folder records the October 21, 2026 Cleaver pilot for one paired-end
control sample. The raw and cleaned FASTQ files remain in ignored local storage
because they are several gigabytes each. The small JSON report is kept here so
the run can be reviewed through Git.

## Command

```text
cleaver fastp \
  -i data/raw/SRR36710387_1.fastq \
  -I data/raw/SRR36710387_2.fastq \
  -o preprocessing/runs/SRR36710387/trimmed/SRR36710387.R1.fastq.gz \
  -O preprocessing/runs/SRR36710387/trimmed/SRR36710387.R2.fastq.gz \
  -j preprocessing/runs/SRR36710387/qc/SRR36710387.fastp.json \
  -2 -f 1 -F 1
```

The command removed one base from the 5-prime end of both mates, used
paired-end overlap analysis to detect adapters, and otherwise used Cleaver's
default quality and length filters.

## Results

| Measurement | Result |
|---|---:|
| Input read pairs | 49,477,059 |
| Output read pairs | 49,476,792 |
| Removed read pairs | 267 |
| Read-pair retention | 99.999460% |
| Q20 before / after | 95.37% / 95.89% |
| Q30 before / after | 88.41% / 89.03% |
| Reads with adapter trimming | 8,529 (0.008619%) |
| Low-quality reads reported | 237 |
| Reads with too many N bases | 0 |
| Reads reported as too short | 30 |

Both cleaned gzip files passed an integrity check. These results support using
conservative cleaning: the reads were already high quality and adapter
contamination was minimal.

## Post-cleaning quality check

FastQC 0.13.0 was run locally on both cleaned mates, and MultiQC 1.35 combined
the two reports. The complete report is
[`post_cleaning_multiqc/SRR36710387_post_cleaver_multiqc.html`](post_cleaning_multiqc/SRR36710387_post_cleaver_multiqc.html).

```text
fastqc --threads 2 --outdir post_cleaning_fastqc \
  SRR36710387.R1.fastq.gz SRR36710387.R2.fastq.gz

multiqc post_cleaning_fastqc \
  --outdir post_cleaning_multiqc \
  --filename SRR36710387_post_cleaver_multiqc.html
```

| Check | R1 | R2 |
|---|---:|---:|
| Cleaned reads | 49,476,792 | 49,476,792 |
| Median length | 100 bp | 100 bp |
| Mean quality at the new first base | Q34.27 | Q35.73 |
| GC content | 49% | 50% |
| Per-base sequence quality | PASS | PASS |
| Adapter content | PASS | PASS |
| Per-base N content | PASS | PASS |
| Per-base sequence content | FAIL | WARN |
| Sequence duplication | FAIL (59.16%) | FAIL (59.34%) |

The original first base had much lower mean quality (Q20.84 for R1 and Q19.59
for R2). After the one-base trim, the new first base is above Q34 in both mates,
and the FastQC per-base quality check changes from FAIL to PASS. Adapter and N
content also pass. The sequence-content warning/failure and duplication signal
remain. Those signals can occur in RNA-seq libraries and will be monitored
across samples; they do not by themselves justify discarding this sample. The
higher post-trim duplication estimate should also be interpreted cautiously
because trimming the variable first base can make previously distinct reads
look identical, and the raw and post-clean reports were produced with different
FastQC versions.

This pilot passes the post-cleaning quality check and is ready for alignment.
The same settings should still be checked on the remaining seven samples before
the full cohort proceeds.

## Local output checksums

```text
f94bdb4fb8f283ab03a4861e3bf44404bec0835f91fa8d93a1322e41f8e8acc2  SRR36710387.R1.fastq.gz
b0dacc0f6fcd403db34f4034b577a4ab14110592a1d0e5ca7b6d10641495e88e  SRR36710387.R2.fastq.gz
```
