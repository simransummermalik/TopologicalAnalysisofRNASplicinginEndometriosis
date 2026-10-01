# RNA-seq preprocessing

This folder contains Splice Girl's small workflow for turning one raw FASTQ
sample into the junction table accepted by the Rust analyzer.

It uses [Cleaver-X](https://github.com/raw-lab/cleaver-x) as an external
bioinformatics tool. Cleaver-X supplies the `cleaver fastp` read-cleaning step
and `cleaver samtools index` BAM-indexing step. The local Python runner is named
`run_preprocessing.py`; it coordinates tools but is not Cleaver-X itself.

```text
raw FASTQ
  -> Cleaver-X read trimming, filtering, and JSON QC
  -> optional Bowtie2 PhiX removal
  -> STAR splice-aware alignment
  -> Cleaver-X BAM indexing
  -> project helper converts STAR junctions to a Splice Girl TSV
```

## What is implemented

- `run_preprocessing.py check` reports missing external programs.
- `run_preprocessing.py plan` validates inputs and prints the exact commands
  without processing reads.
- `run_preprocessing.py run` executes one single-end or paired-end sample and
  records the commands, settings, tool versions, and status in a manifest.
- `build_junction_gene_map.py` builds a junction-to-gene map from adjacent
  exons in the same GTF transcript.
- `star_junctions.py` converts annotated STAR junctions into Splice Girl's
  five-column input and writes conversion QC as JSON.

## Install the command-line tools

Cleaver-X's official installation command is:

```sh
cargo install cleaver
cleaver doctor
```

The remaining dependencies are listed in `environment.yml`:

```sh
conda env create -f preprocessing/environment.yml
conda activate splice-girl-preprocessing
```

The runner expects `cleaver` and `STAR` on `PATH`. Bowtie2 is required only
when `--phix-index` is used.

On Apple Silicon, STAR 2.7.11b may complete while reading zero input sequences.
Check STAR's `Log.final.out` before trusting an alignment. This upstream issue
is tracked at <https://github.com/alexdobin/STAR/issues/2632>. The project has
not yet validated a replacement aligner, so a zero-read STAR run must be
treated as failed preprocessing.

## 1. Check the computer

From the repository root:

```sh
python3 preprocessing/run_preprocessing.py check
```

If PhiX filtering will be used:

```sh
python3 preprocessing/run_preprocessing.py check --with-phix
```

## 2. Build the annotation map

The GTF must match the genome build used to create the STAR index.

```sh
python3 preprocessing/build_junction_gene_map.py \
  --gtf references/annotation.gtf \
  --output references/junction_gene_map.tsv
```

The map uses STAR's one-based intron coordinates. The converter omits
junctions that map ambiguously to more than one gene.

## 3. Preview one sample

```sh
python3 preprocessing/run_preprocessing.py plan \
  --sample-id C01 \
  --read1 data/raw/C01_R1.fastq.gz \
  --read2 data/raw/C01_R2.fastq.gz \
  --star-index references/star_index \
  --gene-map references/junction_gene_map.tsv \
  --output-dir preprocessing/runs/C01
```

For paired-end reads, the runner asks Cleaver-X to detect adapters from read
overlap. Add `--phix-index references/phix/phix` only when PhiX filtering is
needed and the matching Bowtie2 index exists. Junction counts use uniquely
mapped reads unless `--include-multimapping` is supplied.

## 4. Run one sample

After checking the printed plan, replace `plan` with `run`:

```sh
python3 preprocessing/run_preprocessing.py run \
  --sample-id C01 \
  --read1 data/raw/C01_R1.fastq.gz \
  --read2 data/raw/C01_R2.fastq.gz \
  --star-index references/star_index \
  --gene-map references/junction_gene_map.tsv \
  --output-dir preprocessing/runs/C01
```

The runner refuses to overwrite its main outputs. When every command succeeds,
the Splice Girl input is:

```text
preprocessing/runs/C01/splice_girl/C01.junctions.tsv
```

Analyze it with:

```sh
cargo run -- preprocessing/runs/C01/splice_girl/C01.junctions.tsv
```

## Output directories

| Directory | Contents |
|---|---|
| `trimmed/` | Cleaver-X-cleaned FASTQ files |
| `cleaned/` | Reads not aligned to PhiX, if that optional step is enabled |
| `aligned/` | STAR logs, sorted BAM, BAM index, and `SJ.out.tab` |
| `qc/` | Cleaver-X JSON and junction-conversion JSON |
| `splice_girl/` | Five-column junction TSV |
| `run_manifest.json` | Inputs, settings, commands, tool versions, outputs, and status |

Raw reads, references, indexes, BAM files, and run directories are excluded
from Git because they can be large. Keep their checksums and storage locations
in the sample record.
