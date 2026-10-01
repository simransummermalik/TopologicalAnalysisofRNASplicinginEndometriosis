# Cleaver preprocessing

Cleaver is Splice Girl's small preprocessing wrapper. It turns one raw FASTQ
sample into cleaned reads, a sorted/indexed alignment, a STAR splice-junction
table, and the five-column TSV accepted by Splice Girl.

```text
FASTQ
  -> fastp trimming and before/after QC report
  -> optional Bowtie2 PhiX removal
  -> STAR splice-aware alignment
  -> sorted BAM and STAR SJ.out.tab
  -> annotated Splice Girl junction TSV
```

Cleaver is the project wrapper name. It is not the unrelated Bioconductor
`cleaver` protein-cleavage package.

## What is working now

- `cleaver.py check` reports missing external programs.
- `cleaver.py plan` validates the paths and prints the exact commands without
  processing the reads.
- `cleaver.py run` executes one single-end or paired-end sample.
- `build_junction_gene_map.py` builds the required junction-to-gene mapping
  from adjacent exons in the same transcript in a GTF file.
- `star_junctions.py` converts annotated STAR junctions into Splice Girl's
  input format and writes a small QC JSON file.
- The included smoke tests exercise the map builder, strand-aware conversion,
  filtering, and command planning without requiring large reference files.

The external bioinformatics programs are not installed automatically. The
provided `environment.yml` records the required packages, but creating the
environment should wait until installation is approved.

## 1. Check the computer

From the repository root:

```sh
python3 preprocessing/cleaver.py check
```

If PhiX filtering has been approved for the dataset:

```sh
python3 preprocessing/cleaver.py check --with-phix
```

## 2. Build the annotation map

The GTF must match the genome build used to create the STAR index.

```sh
python3 preprocessing/build_junction_gene_map.py \
  --gtf references/annotation.gtf \
  --output references/junction_gene_map.tsv
```

The map uses STAR's one-based intron coordinates. Junctions assigned to more
than one gene are treated as ambiguous and omitted by the converter.

## 3. Preview one sample

Paired-end example:

```sh
python3 preprocessing/cleaver.py plan \
  --sample-id C01 \
  --read1 data/raw/C01_R1.fastq.gz \
  --read2 data/raw/C01_R2.fastq.gz \
  --star-index references/star_index \
  --gene-map references/junction_gene_map.tsv \
  --output-dir preprocessing/runs/C01
```

Add `--phix-index references/phix/phix` only when PhiX filtering is needed and
the matching Bowtie2 index exists. The default count uses uniquely mapped
reads only. Multimapping support is excluded unless
`--include-multimapping` is explicitly supplied.

## 4. Run one approved sample

After reviewing the plan, replace `plan` with `run`. Cleaver refuses to
overwrite its main outputs.

```sh
python3 preprocessing/cleaver.py run \
  --sample-id C01 \
  --read1 data/raw/C01_R1.fastq.gz \
  --read2 data/raw/C01_R2.fastq.gz \
  --star-index references/star_index \
  --gene-map references/junction_gene_map.tsv \
  --output-dir preprocessing/runs/C01
```

The final input will be:

```text
preprocessing/runs/C01/splice_girl/C01.junctions.tsv
```

Run it through the current Rust program with:

```sh
cargo run -- preprocessing/runs/C01/splice_girl/C01.junctions.tsv
```

## Output directories

| Directory | Contents |
|---|---|
| `trimmed/` | Adapter- and quality-trimmed FASTQ files |
| `cleaned/` | Reads not aligned to PhiX, when that optional step is enabled |
| `aligned/` | STAR logs, sorted BAM, and `SJ.out.tab` |
| `qc/` | fastp HTML/JSON and junction-conversion QC JSON |
| `splice_girl/` | Five-column junction TSV |
| `run_manifest.json` | Input paths, options, commands, tool paths/versions, and status |

Raw reads, reference indexes, BAM files, and run directories are deliberately
excluded from Git because they can be very large. Keep their checksums and
storage locations in the sample manifest.
