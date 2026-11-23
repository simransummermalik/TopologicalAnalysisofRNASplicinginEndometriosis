# Data methods

This file records the reference and preprocessing choices that must stay the
same across every GSE315857 sample. Large reference files are stored locally in
the ignored `references/` directory; this document is the reviewable record
that belongs in Git.

## Human genome reference

Splice Girl uses the **GENCODE release 50 GRCh38.p14 primary assembly**. GENCODE
lists release 50 as the current GRCh38 release and provides a primary-assembly
FASTA whose sequence names match its primary-assembly annotation.

| Field | Recorded value |
|---|---|
| Provider | GENCODE / EMBL-EBI |
| GENCODE release | 50 |
| Assembly | GRCh38.p14 |
| Scope | Primary assembly: chromosomes and scaffolds; alternate loci, haplotypes, and assembly patches excluded |
| Filename | `GRCh38.primary_assembly.genome.fa.gz` |
| Local path | `references/gencode_v50/GRCh38.primary_assembly.genome.fa.gz` |
| Download date | October 22, 2026 |
| Compressed size | 845,635,028 bytes |
| Uncompressed size | 3,151,417,447 bytes |
| Official and observed MD5 | `da1a11258be075cfa7af718162c894e7` |
| Observed SHA-256 | `b760d18dbb651dd14dfc290083371b3ef3bff122d43a9cefb13ca4ecf38f05ca` |
| Gzip integrity | PASS |

Exact download URL:

```text
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh38.primary_assembly.genome.fa.gz
```

The original FASTA remains compressed and ignored by Git. On November 17 a
temporary expanded copy was used to build the full HISAT2 index, then removed
after successful indexing. The original compressed file was preserved.

## Matching annotation

Splice Girl uses the **GENCODE release 50 comprehensive primary-assembly GTF**.
It was downloaded from the same GENCODE release as the genome FASTA and checked
before being recorded in the sample manifest.

| Field | Recorded value |
|---|---|
| Annotation | GENCODE v50 comprehensive primary-assembly GTF |
| Filename | `gencode.v50.primary_assembly.annotation.gtf.gz` |
| Local path | `references/gencode_v50/gencode.v50.primary_assembly.annotation.gtf.gz` |
| Download date | October 22, 2026 |
| Compressed size | 124,650,284 bytes |
| Official and observed MD5 | `289b91e5e95e8b0450d223246f10a12e` |
| Observed SHA-256 | `89bbad69a8c89fee5fadec0a2f14a098752d27a4ef7f24e4de9bac681e1b18f4` |
| Gzip integrity | PASS |
| Annotation rows | 11,254,876 total; 5,087,789 exon rows |

Exact download URL:

```text
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.primary_assembly.annotation.gtf.gz
```

Using the primary-assembly FASTA and primary-assembly GTF from the same GENCODE
release avoids chromosome-name and coordinate mismatches. The FASTA contains
194 sequence names and the GTF uses 70 of them. Every GTF sequence name exists
in the FASTA; there are no annotation sequence names missing from the genome.
The 124 FASTA sequences without GTF rows are unannotated reference sequences,
not naming conflicts.

Like the FASTA, the compressed GTF is ignored by Git. The checksums, source,
and compatibility result above are the reproducible record kept in the
repository.

## HISAT2 index and first real technical pilot

The November 17 full-human index is complete at the ignored prefix
`references/hisat2_v50_pilot_nov17/genome`, built using HISAT2 **2.2.3**, two
threads and `--bmaxdivn 16`. The matching exon map has **698,612 junction keys**
and **1,229,403 target rows**; repeated keys retain their alternative targets
for ambiguity checks. The source FASTA/GTF SHA-256 values match the reference
record above. Eight index files, map and splice-site hints have output hashes
in `reference_manifest.json`.

The [bounded real pilot](../data/pilots/SRR36710387_prefix1m_nov17/README.md)
uses the first one million Cleaver-cleaned read pairs from one control,
`SRR36710387`. It reports 97.23% overall alignment and verifies CDH11 junction
counts and decomposition independently. Its unstranded alignment, 8-base
anchors and support thresholds are exploratory technical settings, not final
cohort approval. In particular, **73.7% of assignment-stage junction
observations are excluded as ambiguous** under the frozen map rule. Review
that loss of coverage before expanding; the policy was not relaxed here.

Exact reproduction commands and scope are in [real_pilot.md](real_pilot.md).

## HGNC symbol reference

The candidate-gene mapping in [`gene_sets/gene_sets.py`](../gene_sets/gene_sets.py)
uses the HGNC complete set for approved, previous, and alias symbols. It is
used together with the GENCODE v50 GTF above, which remains the source of the
versioned `gene_id`.

| Field | Recorded value |
|---|---|
| Provider | HUGO Gene Nomenclature Committee (HGNC) |
| Filename | `hgnc_complete_set.txt` |
| Local path | `references/hgnc/hgnc_complete_set.txt` |
| Download date | November 12, 2026 |
| Latest `date_modified` in file | 2026-11-09 |
| Size | 17,006,236 bytes; 45,388 lines including header |
| MD5 | `7c69fba9cf96520ccbeead1ada4ca1bd` |
| SHA-256 | `86cc190fcae1abcd77a13ddba9a34327150477f4e3c6bcadd7f445b26dd31246` |

Exact download URL:

```text
https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
```

HGNC replaces this file regularly and the URL keeps no versions, so the
checksum above identifies the copy actually used. Keep the same copy for the
whole extraction. If a newer copy is downloaded, record it here and re-run the
mapping.

## Frozen junction and graph conventions

These rules apply to every real GSE315857 sample. They keep the biological
meaning of a node or edge from changing between samples.

### Coordinates

- GENCODE GTF exon coordinates are **1-based and closed**: both `start` and
  `end` are included in the exon.
- A raw splice junction is represented as
  `(sequence, intron_start, intron_end, strand)`, where the two intron
  coordinates are also 1-based and closed.
- For a SAM/BAM alignment, `POS` is 1-based. A CIGAR `N` operation consumes
  the skipped reference interval; its first and last skipped bases become
  `intron_start` and `intron_end`.
- STAR `SJ.out.tab` uses the same 1-based, closed intron interval. The existing
  STAR conversion fixture remains useful for checking the convention even
  though HISAT2 is the selected aligner for the real pilot.
- Sequence names are copied exactly from the GENCODE files, including the
  `chr` prefix. No automatic renaming between `1` and `chr1` is allowed.

The coordinate rules follow the
[GENCODE GTF format](https://www.gencodegenes.org/pages/data_format.html), the
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf), and the
[STAR junction specification](https://github.com/alexdobin/STAR/blob/master/extras/doc-latex/STARmanual.tex).

### Graph nodes

**Rule since 2026-11-18 (`shared-splice-site`).** One graph node is one
**exon locus within a gene**: annotated exons of the same gene and strand that
share a start coordinate or an end coordinate are one node. A junction read
shows only its two splice sites, so exon variants that differ only at their
far ends (truncated, alternative-start or alternative-end transcripts) are the
same exon to the data. Exons that share no coordinate stay separate, so
skipped, included and mutually exclusive exons remain distinct nodes, and an
annotated exon-skipping pattern still forms the intended triangle.

A **retained-intron exon** (an exon that strictly contains an annotated intron
of its own gene) would otherwise glue the exons on either side of that intron
into one node. It therefore joins its start to its end only when that does not
merge two nodes that each already contain a normal exon.

The node identifier is:

```text
sequence:lowest_start-highest_end:strand
```

taken over the exons in the node. For example, exons `chr1:300-399:+` and
`chr1:300-420:+` (an alternative donor) form node `chr1:300-420:+`. Every node
contains at least one whole exon, and its lowest coordinate is a start site
that belongs to that node only, so identifiers are unique within a gene.
`exon_number` is not used because it can differ between transcripts.
`gene_id` is stored separately and preserves the versioned GENCODE identifier,
such as `ENSG00000123456.7`. Nodes come from the annotation alone, so they are
identical for every sample.

Why: under the earlier rule, 30% of annotated candidate-gene junctions and
74% of the pilot's junction observations matched several exon pairs and were
excluded. Under this rule 98.5% of candidate junctions map to one node pair.
The comparison is in [the node-rule review](#node-rule-review-2026-11-18).

**Earlier rule (2026-10-22 to 2026-11-17, `exon-interval`).** One node per
unique annotated exon interval. The November 17 technical pilot used it.
`build_junction_gene_map.py --node-rule exon-interval` and
`hisat2_pipeline.py prepare-reference --node-rule exon-interval` reproduce it
exactly.

### Directed edges and strand

An edge is directed in biological transcript order, from the upstream exon to
the downstream exon. Genomic coordinate order alone never determines the
direction.

| Strand | Source (`from`) | Destination (`to`) | Direction on the reference |
|---|---|---|---|
| `+` | Node with an exon whose `end = intron_start - 1` | Node with an exon whose `start = intron_end + 1` | Lower to higher coordinates |
| `-` | Node with an exon whose `start = intron_end + 1` | Node with an exon whose `end = intron_start - 1` | Higher to lower coordinates |

For plus-strand exons `chr1:100-199:+` and `chr1:300-399:+`, the intron is
`chr1:200-299:+` and the edge is `chr1:100-199:+` to
`chr1:300-399:+`. For an otherwise identical minus-strand pair, the edge runs
from `chr1:300-399:-` to `chr1:100-199:-`.

Junctions with unknown strand are not guessed into a direction. They are
excluded from the primary graph and counted in conversion QC.

### Annotation mapping and exclusions

The annotation map is built from exon pairs that are adjacent within at least
one GENCODE transcript, and lists every annotated intron. A junction is kept
only when its sequence, intron coordinates, and strand identify exactly one
`(gene_id, from_node, to_node)` mapping after duplicate mappings are
collapsed.

- No annotation match: record the junction as unannotated and exclude it.
- More than one gene (or, under the earlier rule, exon-pair) match: record it
  as ambiguous and exclude it.
- Both ends in one node (`same_node`): exclude it, because a self-loop has no
  graph meaning. This is rare: 40 of 711,654 GENCODE v50 map rows and none in
  the candidate genes.
- Exactly one match: keep it and aggregate supporting reads for that directed
  edge within the sample. Two junctions that join the same pair of nodes, such
  as alternative 5′ or 3′ splice sites, are **summed into one edge**; the
  choice between nearby sites is not measured.
- Primary counts use uniquely mapped junction-spanning reads. Multimapping
  support remains excluded unless a later reviewed rule explicitly changes
  that setting for every sample.

These exclusions favor a graph with clear biological node identities over
silently assigning a junction to the wrong exon or gene. The conversion QC
reports input-record filters separately from junction-assignment outcomes. Its
assignment denominator is the number of anchor-qualified junction observations
remaining after primary-alignment, unique-alignment, and known-strand filters.
Assigned, ambiguous, unannotated, and same-node observations are reported as
counts and fractions of that denominator; it is not a percentage of all reads.
Ambiguous mappings are also split into cross-gene and within-one-gene cases,
with the number of candidate targets recorded in the QC JSON. The converter
does not change its counting rules based on these summaries.

### Node-rule review (2026-11-18)

Compared on the full GENCODE v50 annotation, without reads, for the 661
candidate genes (18,052 annotated junctions). Splice-site hints for HISAT2 are
identical under both rules, so an existing index stays valid; only the map
file changes.

| Rule | Junctions mapped to one node pair | Same-node | Candidate genes with ≥ 1 graph cycle | Median cycle rank |
|---|---:|---:|---:|---:|
| Exact exon interval (earlier rule) | 69.6% | 0 | 298 | 0 |
| **Shared splice site, retained introns kept apart (adopted)** | **98.5%** | **0** | **536** | **4** |
| *Exploratory:* shared splice site, retained introns may merge | 93.1% | 976 | 526 | 3 |
| *Exploratory:* exact exon interval, GENCODE `basic` transcripts only | 77.9% of a set that drops 24% of junctions | 0 | 313 | 0 |

The first two rows are reproduced from the real maps:

```sh
python3 preprocessing/build_junction_gene_map.py --gtf <GENCODE v50 GTF> \
  --output maps/shared.tsv
python3 preprocessing/build_junction_gene_map.py --gtf <GENCODE v50 GTF> \
  --output maps/exon_interval.tsv --node-rule exon-interval
python3 scripts/node_rule_review.py --genes data/gene_sets/candidate_genes.tsv \
  --map exon-interval=maps/exon_interval.tsv --map shared-splice-site=maps/shared.tsv
```

The exploratory rows came from a one-off comparison script and were not
adopted. A graph-cycle count assumes every usable annotated junction is
observed; real graphs contain only observed junctions.

The remaining 1.5% match several genes and stay excluded. In the November 17
pilot prefix, the candidate genes had 19,983 assigned and 36,591 within-gene
ambiguous observations, so the adopted rule keeps about 100% instead of 35% of
same-gene candidate reads; candidate genes with at least 100 reads rise from
14 to 65. These are diagnostic estimates from the pilot's saved audit, not a
rerun.

### Interpretation boundary

Edge direction records transcript order. The Hodge cycle component measures
graph cycle-space structure after the edges are assembled; it does not require
a literal directed return path. Therefore a nonzero cycle component can
describe alternative paths such as exon skipping and does not by itself prove
circular RNA or a back-splice junction.

## D5 completion and remaining processing review

The genome, annotation, coordinate system, node identity, and strand direction
are frozen. The full human index, map and known splice-site hints were built
and checksum-recorded from this exact pair on November 17, completing D5.
The junction-ambiguity rule was settled on 2026-11-18 (see
[graph nodes](#graph-nodes)), and library strandedness on the same day (see
[below](#library-strandedness-2026-11-18)). Support/coverage eligibility
remains a processing decision; the technical pilot does not approve it.

## Library strandedness (2026-11-18)

All eight runs are **RF**: read 1 is antisense and read 2 sense to the
transcript, in 98.1–99.5% of reads matching nine reference mRNAs (both
groups; [results](../data/strandedness_2026-11-18/README.md)). The pilot's
HISAT2 alignment of SRR36710387 independently gave 99.7%. The cohort alignment
therefore uses `--strandedness RF` (HISAT2 `--rna-strandness RF`).

After alignment, each sample's spliced-alignment orientation is rechecked. A
sample whose read 1 is antisense in less than 95% of eligible spliced
alignments is held for review rather than analyzed. Eligible alignments are
primary, uniquely mapped (`NH:i:1`), spliced and carry an `XS` transcript
strand; a read is antisense when it aligned to the strand opposite to `XS`.
Since 2026-11-23 the junction counter records these counts in each sample's
`junction_qc.json` (`strand_orientation`), and the pipeline writes the pass or
hold call to `run_manifest.json` (`strand_check`) and prints a warning on a
hold. Samples aligned before then are checked from their saved BAM with
`preprocessing/strand_orientation.py` (see the
[HISAT2 guide](../preprocessing/HISAT2_GUIDE.md#check-strandedness-after-alignment)).

The November 17 pilot was aligned as `unstranded`; its results are kept as
recorded and not relabeled.
