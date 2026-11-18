# Run the bounded real-data pilot

This pilot uses one real control sample, `SRR36710387`. It is a technical check
of alignment, junction counting and mathematics. It cannot estimate a disease
effect, validate the entire sample, or replace the eight-sample study.
The [written plan](../data/pilots/SRR36710387_prefix1m_nov17/plan.md) records
the exploratory choices and gene-selection rule before the Hodge result.

**Node rule:** this pilot used the earlier `exon-interval` node rule. Since
2026-11-18 the default is `shared-splice-site`
([decision](decisions.md)). To reproduce the pilot exactly, add
`--node-rule exon-interval` to `prepare-reference` in step 1.

**Completed run:** [CDH11 review and findings](../data/pilots/SRR36710387_prefix1m_nov17/README.md)
and [generated report](../data/pilots/SRR36710387_prefix1m_nov17/verified/report.html).
The technical checks pass; the review records the zero-cycle result and the
substantial ambiguity exclusions that still need a team decision.

Run from the repository root. Outputs must use new folder names for a repeat.
The commands below use the locally available cleaned files; they do not
download reads. Cleaver and HISAT2 must be on `PATH`.

Reuse an existing reference bundle when its manifest is `complete` and its
recorded source hashes match the verified files; skip step 1 in that case.
The human index is built once, not once per sample. A verified subset can also
be reused. New alignment and verification runs need new output folders.

## 1. Prepare the verified full reference

```sh
python3 preprocessing/hisat2_pipeline.py prepare-reference \
  --fasta references/gencode_v50/GRCh38.primary_assembly.genome.fa.gz \
  --gtf references/gencode_v50/gencode.v50.primary_assembly.annotation.gtf.gz \
  --output-dir references/hisat2_v50_pilot_nov17 \
  --threads 2 --build-bmaxdivn 16 --discard-expanded-reference
```

This is the expensive step. `--build-bmaxdivn` changes the suffix-sort memory
tradeoff, not the reference sequence. The expanded FASTA copy is removed only
after successful indexing; the original compressed FASTA is preserved. Failed
builds retain their logs and working copy. Use a complete reference bundle
with a `complete` manifest before alignment; eight filenames alone are not
proof that a still-running build is ready.

## 2. Verify and extract the first million paired reads

```sh
python3 preprocessing/subsample_pairs.py \
  --read1 preprocessing/runs/SRR36710387/trimmed/SRR36710387.R1.fastq.gz \
  --read2 preprocessing/runs/SRR36710387/trimmed/SRR36710387.R2.fastq.gz \
  --pairs 1000000 \
  --expected-sha256-r1 f94bdb4fb8f283ab03a4861e3bf44404bec0835f91fa8d93a1322e41f8e8acc2 \
  --expected-sha256-r2 b0dacc0f6fcd403db34f4034b577a4ab14110592a1d0e5ca7b6d10641495e88e \
  --output-dir preprocessing/runs/SRR36710387_pilot_nov17/subset
```

The extractor checks source hashes, synchronized read names and complete
four-line FASTQ records. It saves exactly the requested number of pairs and a
checksummed manifest. This is a deterministic prefix, not random subsampling.
The original input contains 49,476,792 cleaned pairs; the prefix is about 2.02%.

## 3. Align and count annotated junctions

```sh
python3 preprocessing/hisat2_pipeline.py run \
  --sample-id SRR36710387 \
  --read1 preprocessing/runs/SRR36710387_pilot_nov17/subset/R1.fastq.gz \
  --read2 preprocessing/runs/SRR36710387_pilot_nov17/subset/R2.fastq.gz \
  --index references/hisat2_v50_pilot_nov17/genome \
  --gene-map references/hisat2_v50_pilot_nov17/junction_exon_map.tsv \
  --splice-sites references/hisat2_v50_pilot_nov17/splice_sites.tsv \
  --strandedness unstranded --min-anchor 8 --min-unique 1 \
  --skip-cleaning --threads 2 \
  --output-dir preprocessing/runs/SRR36710387_pilot_nov17/alignment
```

Cleaning is skipped because the original files already passed through Cleaver.
The unstranded alignment and anchor/count values are exploratory pilot settings,
not a decision about the final library protocol or cohort thresholds. Preserve
the alignment summary, exclusion QC and manifests alongside the large files.

## 4. Select and independently verify a real gene

The verifier needs NumPy in addition to the tools above. Use a Python
environment containing NumPy; it remains an optional validation dependency,
not a dependency of the Rust analyzer.

```sh
cargo build --release
python scripts/validate_real_pilot.py \
  --sample-id SRR36710387 \
  --junctions preprocessing/runs/SRR36710387_pilot_nov17/alignment/junctions.tsv \
  --candidates data/gene_sets/candidate_genes.tsv \
  --gene-map references/hisat2_v50_pilot_nov17/junction_exon_map.tsv \
  --bam preprocessing/runs/SRR36710387_pilot_nov17/alignment/aligned.bam \
  --sample-manifest data/metadata/samples.tsv \
  --subset-manifest preprocessing/runs/SRR36710387_pilot_nov17/subset/subset_manifest.json \
  --alignment-manifest preprocessing/runs/SRR36710387_pilot_nov17/alignment/run_manifest.json \
  --plan data/pilots/SRR36710387_prefix1m_nov17/plan.md \
  --min-support 3 --min-anchor 8 --min-edges 3 --max-nodes 12 --min-total 100 \
  --output-dir data/pilots/SRR36710387_prefix1m_nov17/verified
```

Open the generated `report.html`. Verification includes:

1. Source/subset/alignment identity checks and a BAM integrity check.
2. Candidate-gene selection by coverage and graph size, with gene-ID ties.
3. A separate recount from indexed BAM, preserving annotation ambiguities.
4. The actual Rust analysis and three-panel graph.
5. An independent SVD-based least-squares projection, compared edge by edge
   with the Rust output, plus reconstruction, kernel and orthogonality checks.
6. The same gene recalculated in Rust and NumPy at supports 1, 3, 5 and 10.

The report records missing inference explicitly. Passing these checks establishes
agreement on this technical input; it does not approve biological interpretation,
library strandedness, final thresholds, all genes or the remaining samples.

Index-building options follow the [HISAT2 manual](https://daehwankimlab.github.io/hisat2/manual/).
