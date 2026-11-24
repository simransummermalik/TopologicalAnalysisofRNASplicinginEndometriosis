# Methods draft

Updated November 24, 2026 (first assembled November 17). Every setting below
was fixed before the eight-sample results were produced. Completed steps are
written in the past tense; steps that wait on the full cohort are written as
planned. Values that only the cohort run can supply are marked
**[after run: …]**. The controlling records are
[`decisions.md`](decisions.md) and the documents linked in each section; if
this draft disagrees with them, they win.

## Study design and samples

The study used the public bulk RNA-seq dataset GSE315857: proliferative-phase
eutopic endometrium from four women with endometriosis (SRR36710383–86) and
four laparoscopically confirmed controls (SRR36710387–90). Both groups are the
same tissue in the same cycle phase, so the comparison is not a lesion against
a different tissue. Each run comes from a different subject, and the groups
were analyzed as independent samples. Accessions, condition labels and
internal subject IDs are in the [sample manifest](../data/metadata/samples.tsv).
The raw QC summaries list 545,210,220 read pairs across the eight runs
(49.5–80.0 million per run, median read length 100 bases).

The dataset has no patient exposure measurements and no reported sequencing
batch. The written sources also disagree on the library kit and read length;
the [dataset review](dataset_review.md) records this. That disagreement
affects only how the original authors' processing is described: this project's
library orientation was measured directly (see *Library strandedness*).

## Candidate genes

Genes were chosen from two evidence sources, both fixed before any
case/control splicing result existed. Full rules, sources and audits are in
[`gene_set_methods.md`](gene_set_methods.md).

- **Genetic risk (103 genes).** Endometriosis associations from the GWAS
  Catalog and the cross-ancestry loci of Koller et al. 2026 (Supplementary
  Table 4), using the gene each source itself reports. Cross-trait, MTAG and
  non-primary Koller analyses were excluded. Genes were never assigned by
  proximity alone; loci whose authors list several genes keep every listed
  gene, flagged.
- **Environmental response (565 genes).** Curated *Homo sapiens*
  chemical–gene interactions from CTD (release of 2026-10-29) for TCDD and a
  fixed list of dioxin-like compounds, measured responses only. Co-treatment
  rows were excluded because a mixture's effect cannot be credited to the
  dioxin. A gene–chemical pair needed at least two publications, at least one
  of them not a large screen (a publication behind 100 or more evidence rows).
- **Candidate set (661 genes).** The union of both sets, with stable GENCODE
  IDs: 96 risk-only, 558 environment-only and 7 in both (overlap).

Both evidence tables were audited on a random sample. A person other than the
author is checking sampled rows against their sources
([`spot_checks.md`](spot_checks.md)) **[status at submission: …]**.

## Reference and annotation

GENCODE release 50 was used throughout: the GRCh38.p14 primary-assembly FASTA
and the matching comprehensive annotation GTF. File hashes and download dates
are in [`data_methods.md`](data_methods.md). From these, the project built:

- a full HISAT2 index;
- an annotated splice-site file, passed to HISAT2 as known sites;
- a junction-to-graph map: every annotated intron with its gene, strand and
  the two graph nodes it joins (rule under *Junction counting*).

## Read processing and alignment

Each sample was processed as its own job on the UNC Charlotte research
computing cluster (Slurm) with the project's pipeline
(`preprocessing/hisat2_pipeline.py run`; see the
[HISAT2 guide](../preprocessing/HISAT2_GUIDE.md)):

1. **Cleaning:** Cleaver (RAW Lab) running fastp on the paired reads, with
   paired-end adapter detection; the QC report is kept per sample.
2. **Alignment:** HISAT2 2.2.3 with `--dta-cufflinks`, the known splice-site
   file and `--rna-strandness RF`.
3. **Sorting and indexing** with samtools (through Cleaver), either from a
   saved SAM file or streamed directly to BAM; both modes count the same
   records with the same rules.

Each run writes a manifest with tool versions, the exact commands, input and
output checksums and per-stage status. Alignment took 1 h 42 min to 2 h
53 min per sample. **[after run: overall alignment rate per sample, from each
`alignment_summary.txt`.]**

## Library strandedness

The orientation was measured, not taken from the written sources. Before
alignment, the first 200,000 reads of each mate of every run were compared
with both strands of nine highly expressed reference mRNAs (25-mer matching,
`preprocessing/strand_kmer_check.py`). Read 1 was antisense in 98.1–99.5% of
matched reads in every run and both groups: an RF (dUTP-style) library
([check](../data/strandedness_2026-11-18/README.md)).

After alignment, every run is checked again: among primary, uniquely mapped,
spliced reads, the share whose read 1 is antisense to the transcript strand
HISAT2 records (`preprocessing/strand_orientation.py`). A sample below 95%
is held for review and excluded from the cohort table, and the Rust analysis
refuses a held sample as a second safeguard. Runs made before this check was
built are checked from their saved BAM files. **[after run: read-1 agreement
per sample; samples excluded, if any.]**

## Junction counting

Junction reads were counted from the alignments (`preprocessing/sam_junctions.py`):

- **Reads counted:** primary alignments that map uniquely (`NH:i:1`), have a
  known transcript strand (`XS` tag) and span an intron with at least 8
  aligned bases on each side. Duplicate reads were kept. Every junction with
  at least one read was written; read-support thresholds are applied later,
  across samples.
- **Coordinates:** 1-based closed intron intervals, with sequence names exactly
  as in GENCODE.
- **Graph nodes (`shared-splice-site` rule):** within a gene and strand,
  annotated exons that share a start or end coordinate are one node. A
  junction read sees only its two splice sites, so exons differing only at
  their far ends are the same exon to the data. A retained-intron exon may
  not merge two nodes that each hold a normal exon. Node IDs are
  `chromosome:lowest-highest:strand`, from the annotation alone, so they are
  identical in every sample. Under this rule 98.5% of annotated candidate
  junctions map to one node pair, against 69.6% under one node per exon
  interval ([review](data_methods.md#node-rule-review-2026-11-18)).
- **Edges:** directed in transcript order (upstream to downstream exon).
  Junctions joining the same node pair, such as alternative 5′ or 3′ splice
  sites, were summed into one edge.
- **Exclusions, each counted in the run's QC:** junctions with no annotation
  match, with matches in more than one gene, or with both ends in one node.
  No junction was assigned by guessing.

## Combining samples

The per-sample runs were combined into one cohort table by
`preprocessing/assemble_runs.py`, which writes nothing unless every check
passes: every stage completed; every sample used the same settings and the same
reference, splice-site and map files (checked by hash); no count table changed
after its run; and every sample has its own strand-check result. A sample held
by the strand check is listed as excluded with the reason.

## Shared graphs and gene eligibility

To keep sequencing depth from changing a gene's topology, each gene was given
**one shared graph** used by every sample, chosen without the group labels
([threshold proposal](threshold_proposal.md); approved 2026-11-22, confirmed
2026-11-24):

- **Shared edges:** an edge is kept if it has at least 5 reads in at least 4
  of the 8 samples. Because labels are not used, the selection is the same
  under every relabeling in the permutation test.
- **Sample coverage:** a sample is covered for a gene if it has at least 50
  reads on the gene's shared edges. A shared edge with no reads in a covered
  sample is a measured zero; an uncovered sample is missing, not zero.
- **Gene eligibility:** the shared graph must have at least one independent
  cycle (cycle rank = edges − nodes + connected components ≥ 1), and all
  eight samples must be covered. A gene without a cycle has a cycle fraction
  of zero by definition and is not tested.

**[after run: candidate genes with shared edges, with a cycle, and eligible
(tested), by category.]**

## Graph decomposition

Each sample's graph for one gene was analyzed separately
([`math_spec.md`](math_spec.md)). Each edge count was divided by the gene's
total count in that sample, giving the edge signal *F*. With the incidence
matrix *B*₁ (−1 at an edge's source node, +1 at its target), least squares
gave the gradient component *F*grad = *B*₁ᵀφ\*, the part of the signal
explained by node potentials. The remainder *F*cycle = *F* − *F*grad lies in
the cycle space, ker(*B*₁). The **cycle fraction** is
‖*F*cycle‖² / ‖*F*‖². Every result was checked numerically for
reconstruction, zero divergence of the cycle part and orthogonality of the two
parts. The full simplicial decomposition (curl and harmonic parts) was not
used, as the proposal allowed.

Because every covered sample shares the gene's graph, a sample with no reads
on an edge can still have a cycle component on that edge
([shared graphs](shared_graphs.md#interpretation-changes-to-explain)); this
does not show that the sample used that junction.

## Statistics

The rules were fixed on 2026-11-17 ([statistics specification](statistics_spec.md#inference-rules-decided-2026-11-17)):

- **Effect per gene:** mean case cycle fraction minus mean control cycle
  fraction.
- **Primary (confirmatory) test:** one test on the candidate set. The
  statistic is the mean over tested genes of the absolute effect, so genes
  changing in opposite directions do not cancel. All 70 ways of choosing 4 of
  the 8 samples as cases are applied to every gene at once (exact, no random
  seed). Two-sided; significance at p < 0.05. With 4 versus 4, the smallest
  possible p-value is 2/70 ≈ 0.029.
- **Exploratory:** the same set test on the risk, environment and overlap
  sets; gene-level exact tests with Benjamini–Hochberg adjustment across
  tested genes. These are reported with raw p-values and never called
  significant.
- **Sensitivity:** the primary analysis plus six runs that each change one
  threshold (edge reads 3/5/10, samples per edge 2/4/6, coverage
  20/50/100). All seven are reported; none is chosen afterwards.

## Interpretation rules

How results are reported was fixed before the results were produced
([interpretation plan](interpretation_plan.md), approved 2026-11-24 with its
amendment). `scripts/interpret_results.py` applies it:

- **Baseline.** The same set test was run on a summary with no topology:
  junction-usage entropy (normalized Shannon entropy of the read shares over
  the gene's shared edges), on the same genes, samples and counts. The script
  first reproduces the Rust primary result and stops if it does not match.
- **Outcome.** One of four pre-declared outcomes, with fixed wording: no
  difference detected; a difference also seen in the baseline; a difference
  not seen in the baseline; or not tested.
- **Expression check.** Each gene's junction reads per million assigned
  junction reads (log2) went through the same set test and were correlated
  with the cycle-fraction effects, to show whether the effects follow
  expression.
- **Overlap comparison.** A one-sided gene-label permutation tested whether
  overlap genes had larger absolute effects than the other tested candidates;
  it requires at least 3 tested overlap genes.
- **Genes discussed.** Only the five tested genes with the largest absolute
  effect are named, as illustrations. For each, the per-junction share of
  the cycle-fraction change is reported; these shares add up to the gene's
  effect.

## Timing of decisions

The node rule, strandedness and thresholds were written on 2026-11-18. On
2026-11-20 an exploratory comparison of one case and one control (per-sample
graphs, no shared edges, no statistics) was viewed. The thresholds were then
approved on 2026-11-22 without any change and frozen; any later change would
have to be reported as made after case/control data were seen
([decisions, 2026-11-23 correction](decisions.md)). The interpretation plan
was approved on 2026-11-24, before any eight-sample result was produced.

## Software and validation

Splice Girl is a Rust command-line program with no external crates; the
preprocessing and interpretation scripts are Python. The method was first
checked on synthetic graphs with known answers (a linear chain, exon
skipping, a back-splice cycle, added noise and an increasing back-splice
weight), as the proposal required. On real data, one control's first million
read pairs were aligned as a technical pilot, and one gene (CDH11) was
recounted independently from the BAM and recomputed with a NumPy SVD
projection; both matched ([software validation](software_validation.md)).
That pilot used earlier settings (unstranded alignment and one node per exon
interval) and is not part of the cohort result. Every run records input
hashes, parameters, tool versions and exclusions. Reads, alignments and the
reference stay outside Git.

## Limitations

- Four cases and four controls from one dataset: only large, consistent
  differences are detectable, and results describe association, not cause.
- Only uniquely mapped reads at annotated junctions were counted; ambiguous
  and unannotated junctions were excluded. Sequencing depth differs between
  samples (49.5–80.0 million pairs).
- Back-splice junctions were not detected by this pipeline. Cycles in the
  real graphs are alternative paths between annotated exons (for example exon
  skipping), not circular RNA.
- Environmental membership comes from curated responses to dioxin-like
  compounds, mostly in non-reproductive tissue; no exposure was measured in
  these subjects. Risk membership comes from reported GWAS genes; the
  associated variant may act through another gene.
- Cell composition, batch and clinical differences between donors were not
  measured and cannot be modeled with eight samples.

## To fill in after the cohort run

- Alignment rate and read-1 strand agreement per sample; any excluded sample.
- Candidate genes with shared edges, with a cycle and tested, by category.
- Primary p-value and outcome; baseline and expression results.
- Sensitivity: whether the primary outcome changed in any scenario.
- Spot-check and G1 cold-read status at submission.
