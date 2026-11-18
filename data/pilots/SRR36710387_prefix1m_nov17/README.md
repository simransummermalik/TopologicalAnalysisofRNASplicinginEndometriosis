# Real-data technical pilot: CDH11

**Completed November 17, 2026.** Open [the report](verified/report.html).
This is one control sample's bounded technical pilot, not the full study.

## What actually ran

The first **1,000,000 synchronized cleaned read pairs** from `SRR36710387`
were aligned against the full verified GENCODE v50 GRCh38.p14 primary
reference using HISAT2 2.2.3. The existing Cleaver-cleaned source files matched
their recorded SHA-256 values. Cleaver sorted and indexed the resulting BAM.
The complete derivation and exploratory settings are in the [plan](plan.md)
and [reproduction guide](../../../docs/real_pilot.md).

| Measurement | Observed result |
|---|---:|
| Input biological samples | 1 control, 0 cases |
| Cleaned pairs used | 1,000,000 of 49,476,792 (about 2.02%) |
| HISAT2 overall alignment rate | 97.23% |
| Concordantly uniquely aligned pairs | 882,419 (88.24%) |
| Exported annotated junction rows across the prefix | 32,575 |
| Selected gene | CDH11, `ENSG00000140937.16` |
| Candidate-list label | `environment_only` (provided evidence, not donor exposure) |
| Selected graph at minimum support 3 | 6 nodes, 4 junctions |
| Retained read-junction support | 103 |
| Independent BAM recount | PASS: all 6 base-table edges and 106 observations agree |
| Rust/SVD comparison | PASS: largest exported-component difference about 4.6e-11 |
| Cycle rank / cycle fraction | 0 / effectively 0 |

The gene was selected by the rule written before the Hodge calculation:
coverage, graph size, then gene ID for ties. Only one candidate met this
small-graph rule. It was not selected for a large cycle fraction or disease
effect. At support 3, two low-support edges (1 and 2 observations) were removed,
leaving 103 of the 106 assigned gene observations.

## How to read the zero

The retained graph has two connected components and no independent graph loop:
4 edges minus 6 nodes plus 2 components equals cycle rank zero. Its signal
therefore fits entirely in the gradient space. Rust returns zero cycle
fraction; the independent SVD gives approximately 5.9e-31, numerical roundoff
around zero. The same conclusion holds at supports 1, 3, 5 and 10. This checks
the software on a real forest-shaped graph. It does not validate a nonzero
cycle signal in real data or establish the absence of alternative splicing,
back-splicing, circular RNA, disease effects or exposure effects.

## Important review before expanding

The full-prefix junction QC reports **432,599 ambiguous**, **152,907 assigned**,
and **1,583 unannotated** junction observations after the upstream filters.
Ambiguity accounts for **73.7% of these 587,089 observations**. These are
junction observations, not independent reads, read pairs or subjects.

The current frozen map requires exactly one full exon-pair target for a
junction. Shared boundaries across transcript/exon variants can therefore
remove a large fraction of otherwise aligned observations. These exclusions
were recorded and preserved, not reassigned. The team should review their
effect on graph coverage and comparability before approving the full cohort;
the frozen mapping policy was not changed to obtain this pilot result.

Also review library strandedness, final anchors/support/coverage thresholds,
and the deterministic-prefix sampling limit. There is no case sample here,
so no case–control effect, p-value or disease conclusion was calculated.

## Files and provenance

- [Verification JSON](verified/verification.json): source/output hashes, gene
  choice, exact commands, independent recount and SVD errors.
- [Junction table](verified/junctions.tsv) and
  [signed components](verified/results/junction_components.tsv).
- [Support sensitivity](verified/support_sensitivity.tsv): Rust and NumPy were
  rerun on the same gene at every listed support threshold.
- [Alignment summary](verified/alignment_summary.txt),
  [BAM counts](verified/alignment_flagstat.txt) and
  [junction exclusions](verified/junction_qc.json).
- The complete reference, subset and alignment manifests are copied into
  `verified/`. Large FASTQ, index, SAM and BAM files remain in ignored local
  folders. Only the small selected-gene review bundle belongs in Git.

Reference preparation reported 3,916 seconds (about 65 minutes); the bounded
alignment/conversion/BAM workflow reported 213 seconds. The selected six-node
Rust analysis/export took about 0.57 seconds including process startup. These
are single-run timings, not general performance guarantees. The optional
macOS `time -l` profiler could not access restricted kernel counters and made
its outer wrapper return nonzero; both actual pipeline manifests are complete,
BAM integrity passed, and independent verification completed. No peak-memory
claim is made from that profiler.
