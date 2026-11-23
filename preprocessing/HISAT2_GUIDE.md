# Cleaver and HISAT2: the runnable preprocessing workflow

This adapter implements the HISAT2 method already named in the project's data
methods. Cleaver is the external RAW Lab tool, not the name of our pipeline.
The old `run_preprocessing.py` remains a separate STAR prototype.

```text
one sample's FASTQ file(s)
  -> cleaver fastp: clean reads and save QC
  -> HISAT2: align reads using the matching genome and annotated splice sites
  -> sam_junctions.py: count eligible annotated junctions
  -> cleaver samtools: sort and index the alignment
  -> junctions.tsv: Splice Girl's five-column input
```

## 1. Install and check tools

From the repository root, with Cargo and Conda installed:

```sh
cargo install --git https://github.com/raw-lab/cleaver-x --bin cleaver
conda env create -f preprocessing/environment.yml
conda activate splice-girl-preprocessing
python3 preprocessing/hisat2_pipeline.py check
```

Install Cleaver-X from the RAW Lab Git repository; the similarly named crate on
crates.io is a library and does not install the `cleaver` command. The
environment supplies HISAT2 2.2.3 and samtools. Cleaver must also be on `PATH`.
This guide does not require Galaxy or download any reads automatically.

## 2. Test the complete workflow without human data

```sh
python3 scripts/check_project.py --with-tools
```

This runs the Rust and Python checks and an actual external-tool smoke test.
The test builds a tiny synthetic reference, processes two separate synthetic
paired-end samples, merges their junction tables without pooling their counts, and
generates a Rust comparison report. Temporary outputs are removed afterward.
Its anchor/support choices are test settings, not approved study thresholds.
The smoke test checks both retained-SAM and streaming-BAM modes and requires
identical junction tables and QC counters.

## 3. Prepare the real reference once

Use the verified matching GENCODE v50 primary-assembly FASTA and GTF recorded
in [data_methods.md](../docs/data_methods.md). Building the full human index
requires substantially more resources than the smoke test.

```sh
python3 preprocessing/hisat2_pipeline.py plan-reference \
  --fasta references/gencode_v50/GRCh38.primary_assembly.genome.fa.gz \
  --gtf references/gencode_v50/gencode.v50.primary_assembly.annotation.gtf.gz \
  --output-dir references/hisat2_v50 --threads 2
```

Inspect the plan, then replace `plan-reference` with `prepare-reference` to
execute it. The new folder contains the eight-part index (`genome` is its
prefix), exon map, splice-site hints, log, and checksummed reference manifest.
An existing output folder is never overwritten. Use a new folder for a retry.

On a machine with limited disk space, `--discard-expanded-reference` removes
only this run's expanded FASTA copy after a successful build. The original
compressed reference remains intact. `--build-bmaxdivn` exposes HISAT2's positive
suffix-sort bucket divisor for a lower-memory build. See the
[bounded real-data pilot](../docs/real_pilot.md) for recorded example choices;
neither option changes the biological reference.

The exon map uses the `shared-splice-site` node rule (decided 2026-11-18; see
[data_methods.md](../docs/data_methods.md#graph-nodes)). `--node-rule
exon-interval` reproduces the earlier rule. A reference folder built before
2026-11-18 has an `exon-interval` map but valid index and splice-site hints,
which do not depend on the rule. Rebuild only the map, in seconds rather than
the hour an index takes, and pass it with `--gene-map`:

```sh
python3 preprocessing/build_junction_gene_map.py \
  --gtf references/gencode_v50/gencode.v50.primary_assembly.annotation.gtf.gz \
  --output references/hisat2_v50/junction_exon_map_shared_site.tsv
```

## 4. Write a local sample sheet

Make a tab-separated file with these four columns. The paths below are
examples: replace them with files you actually have. Paths are resolved
relative to the sheet's folder. Keep one row per biological sample.

```text
sample_id	condition	read1	read2
SRR36710390	control	raw/SRR36710390_R1.fastq.gz	raw/SRR36710390_R2.fastq.gz
SRR36710386	case	raw/SRR36710386_R1.fastq.gz	raw/SRR36710386_R2.fastq.gz
```

Add all eight rows when their files are available. `condition` must be `case`
or `control`; leave `read2` blank only for genuinely single-end data. A file
cannot be reused as another sample or mate. The rich study metadata remains
in `data/metadata/samples.tsv`; a batch sheet does not establish study eligibility.

## 5. Plan and run the batch

For GSE315857 the strandedness is **RF** (decided 2026-11-18 from all eight
runs; see [data_methods.md](../docs/data_methods.md#library-strandedness-2026-11-18)).
For another dataset, check it with `preprocessing/strand_kmer_check.py` and
the library protocol. Choose the count/anchor thresholds through the group's
methods review. Do not infer these from
whether an example happens to run. `FR`/`RF` describe paired library
strandedness; `unstranded` is an explicit choice, not an automatic guess.

Use `python3 preprocessing/batch_preprocessing.py --help` for all arguments.
For example, the command structure is:

```text
python3 preprocessing/batch_preprocessing.py plan
  --samples LOCAL_SAMPLE_SHEET.tsv
  --index references/hisat2_v50/genome
  --gene-map references/hisat2_v50/junction_exon_map.tsv
  --splice-sites references/hisat2_v50/splice_sites.tsv
  --strandedness RF
  --min-anchor 8
  --min-unique 1
  --output-dir preprocessing/runs/cohort_run_01
  --threads 2
```

Join those lines with shell continuation backslashes when running. Replace
the capitalized placeholders first. `plan` checks every input and prints
commands without starting jobs; change it to `run` to execute the batch.
Samples run sequentially to avoid starting eight alignments at once.
For already cleaned reads, use a separate sheet of cleaned paths and pass
`--skip-cleaning`; the manifest records that choice.

Every sample gets its own logs, QC, BAM/index, junction table and manifest.
The batch publishes a combined `junctions.tsv` and `samples.tsv` only after
every sample completes. A failure stops the batch, records the failed and
unstarted samples, and leaves existing completed results for inspection.
The combined table preserves sample IDs and does not manufacture missing zeros.

Anchor 8 and minimum support 1 are the
[thresholds](../docs/threshold_proposal.md) provisionally approved on
2026-11-22; filtering by edge support happens later, in the shared-graph step.

### Check strandedness after alignment

Each sample's `run_manifest.json` has a `strand_check` entry. `pass` means at
least 95% of read-1 spliced unique alignments are antisense, as expected for
RF. `hold_for_review` means the sample must not enter the analysis until the
team has looked at it; the run also prints a warning.

The batch command leaves a held sample out of the combined `junctions.tsv`
and `samples.tsv`. It lists the sample under `excluded_samples` in
`batch_manifest.json` with the reason and the counts, and keeps the sample's
own folder for review. If every sample is held, nothing is published. The
Rust analysis also refuses any sample whose run manifest says
`hold_for_review`. Putting a reviewed sample back requires a new dated
decision and a new run.

For a sample aligned before 2026-11-23, check its saved BAM:

```sh
samtools view preprocessing/runs/SAMPLE/aligned.bam \
  | python3 preprocessing/strand_orientation.py --sam - --expected RF \
      --sample-id SAMPLE > SAMPLE_strand_check.json
```

The command prints the counts and the call as JSON. It exits with 0 on pass
and 3 on hold. A full BAM takes a few minutes to read.

### Save disk space with streaming BAM

Add these options to either the single-sample or batch `plan`/`run` command:

```text
--storage-mode stream-bam
--sort-memory-mb 256
--min-free-gb 2
```

For example, after preparing a real sample sheet and choosing the reviewed
reference/strand/threshold arguments, append `--storage-mode stream-bam` to
the existing command. Keep `--skip-cleaning` only for already cleaned inputs.
Planning prints the real pipe connections and the disk estimate without
starting processing. Use the same options with `run` when ready.

The two storage choices are:

| Mode | What happens |
|---|---|
| `retain-sam` (default) | Existing workflow: save SAM, count junctions, then sort/index BAM using Cleaver. |
| `stream-bam` | Clean with Cleaver, pipe HISAT2 stdout into samtools sort, check the BAM, stream its records through the same junction counter, then index with Cleaver. No SAM file is created. |

Standalone samtools is already included in `environment.yml`. The streaming
sort uses its memory setting and disk-spill prefix; Cleaver remains the actual
read-cleaning and BAM-indexing tool. `--sort-memory-mb` is an approximate
**per-thread** sort-memory setting in MiB, not a limit on HISAT2 or the whole
pipeline. Sort spill files stay in the new run's output folder. See the
[samtools sort manual](https://www.htslib.org/doc/samtools-sort.html).

In stream mode, BAM/count outputs have `.partial` names until alignment,
BAM checking, full record decoding/counting and indexing succeed. Only then
are the final `aligned.bam`, `aligned.bam.bai` and `junctions.tsv` published.
Always check that `run_manifest.json` says `complete`; a file's presence alone
does not establish a completed run. The quick BAM check covers header/EOF,
while the subsequent record-decoding stage reads through the BAM; quickcheck
alone does not detect all internal corruption. See its
[manual](https://www.htslib.org/doc/samtools-quickcheck.html).

### Disk checks, errors and interrupted runs

Both modes now check the output filesystem **before creating the run folder
or starting tools**. The batch checks the summed sample estimate first, then
each sample checks again before its own run. The estimate uses the input
read-file sizes and these additional-output multipliers:

| Mode | Already cleaned (`--skip-cleaning`) | Cleaning needed |
|---|---:|---:|
| `retain-sam` | 12 × input bytes | 13 × input bytes |
| `stream-bam` | 3 × input bytes | 4 × input bytes |

The estimate includes an allowance for outputs/sort spill; original input
files already occupy disk and are not counted as new allocations. The
`--min-free-gb` reserve is added on top, using GiB (1024 cubed bytes), and is
rechecked between stages. Its default is 2 GiB. It can be adjusted explicitly
but does not disable the output estimate. These conservative heuristics are
planning checks, not guarantees against every compression ratio, mapping
pattern or another process filling the disk. The input files must already
exist; downloads and new reference builds need their own storage allowance.

On a tool failure or Control-C/SIGTERM interruption, the runner stops its
started subprocesses, records `failed` or `interrupted`, and retains logs and
partial files in that run's folder. It checks **both sides of each pipe**, so
a failed HISAT2 process cannot be hidden by a successful sorter, or a failed
BAM reader by a successful counter. Stage records include commands, process
exit codes and separate log filenames. A batch does not publish its complete
cohort bundle when a sample fails or is interrupted.

Choose a **new output folder** for a retry. Automatic resume is not implemented.
No old reads, references or prior run folders are deleted. Samtools manages
its own sort scratch files inside the new run. Reserve checks happen between
stages; an abrupt power loss/forced kill cannot execute Python's interruption
handler, and a completely full disk can prevent final log updates.

### Combine samples that were run as separate jobs

On a cluster it is often easier to run each sample as its own job
(`hisat2_pipeline.py run` per sample) than one long batch. When every sample
has finished, combine them with `assemble_runs.py`. It realigns and recounts
nothing; it checks the finished runs and writes the same `batch_manifest.json`,
`junctions.tsv` and `samples.tsv` as the batch command, so the configured
analysis (`splice-girl check --config`) can verify every sample back to its
reference.

Write a run sheet, one row per sample:

```text
sample_id	condition	run_dir	strand_check
SRR36710390	control	runs/SRR36710390	strand/SRR36710390.json
SRR36710386	case	runs/SRR36710386
```

`run_dir` is the sample's finished run folder. `strand_check` is the JSON from
[the strand check](#check-strandedness-after-alignment) for runs made before
2026-11-23; leave it blank for newer runs, which record the check themselves.
Paths are relative to the sheet. Then:

```sh
python3 preprocessing/assemble_runs.py \
  --samples runs.tsv --output-dir cohort \
  --strandedness RF --min-anchor 8 --min-unique 1
```

Before writing anything it stops if any sample:

- did not complete every stage, or used different strandedness, anchor,
  support or duplicate settings;
- used a different gene map, splice-site file or index from the others;
- has a `junctions.tsv` that changed after its run finished;
- has no strand check, or a strand-check file made for another sample.

The strand call is recomputed from the recorded counts. A sample below 95%
is left out of `junctions.tsv` and `samples.tsv`, listed under
`excluded_samples` with its counts, and reported as a warning. The run
folders are not changed or copied, so keep them in place: the analysis
rechecks them through the absolute paths in `batch_manifest.json`.

## 6. Analyze the combined output

```sh
cargo run --release -- report \
  preprocessing/runs/cohort_run_01/junctions.tsv \
  preprocessing/runs/cohort_run_01/samples.tsv \
  --output target/cohort_report_01
```

Open `target/cohort_report_01/report.html`. Use `stats` only with a reviewed
gene family and explicit statistical settings; see the
[software guide](../docs/running_the_tool.md). A successful command is not a
substitute for alignment, read-quality, coverage, or cohort review.

## Counting rules and limits

- CIGAR `N` yields intron boundaries in 1-based closed coordinates. Insertions
  and clipping do not consume reference positions; deletions do.
- A junction needs the requested contiguous aligned anchor on each side.
- Only primary, non-QC-failed alignments with `NH:i:1` are counted. Transcript
  direction comes from `XS:A:+/-`, not the read's reverse-complement flag.
- HISAT2 runs with `--dta-cufflinks` and known splice-site hints. Unknown-strand,
  nonunique, ambiguous and unannotated observations are excluded with QC counts.
  This deliberately restricts the output to the frozen annotation map.
- The count unit is supporting read alignments. Two mates supporting the same
  junction can contribute two observations; these are not deduplicated
  fragment or molecule counts. Duplicate-flagged records are retained unless
  `--exclude-duplicates` is requested; this does not mark duplicates itself.
  Proper-pair status is not an additional filter.
- This workflow does **not** detect back-splice junctions or establish circular
  RNA. It can still supply alternative-path graphs for Hodge analysis. D7's
  detection limitation is documented; a circular-RNA analysis would need a
  separate reviewed detector and input adapter.
  A missing back-splice in these outputs is therefore not evidence that the
  biological sample lacks back-splicing.
- The default mode retains SAM alongside BAM; `--storage-mode stream-bam`
  avoids writing SAM while retaining BAM, counts, QC and provenance. FASTQs,
  references, indexes, BAMs and run folders stay outside Git.
- Reference preparation does not certify that arbitrary FASTA/GTF files match;
  use the verified pair. Index/map/hint hashes record exactly what was used.

Coordinate and tag definitions follow the [HISAT2 manual](https://daehwankimlab.github.io/hisat2/manual/)
and [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).
