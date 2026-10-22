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
contamination was minimal. A post-cleaning FastQC report is still needed before
the settings are applied to every sample.

## Local output checksums

```text
f94bdb4fb8f283ab03a4861e3bf44404bec0835f91fa8d93a1322e41f8e8acc2  SRR36710387.R1.fastq.gz
b0dacc0f6fcd403db34f4034b577a4ab14110592a1d0e5ca7b6d10641495e88e  SRR36710387.R2.fastq.gz
```
