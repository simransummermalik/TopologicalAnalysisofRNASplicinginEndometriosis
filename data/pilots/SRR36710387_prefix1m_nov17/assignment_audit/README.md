# Why junction observations were excluded

Diagnostic audit of the same one-million-pair control prefix, November 17,
2026. Open [report.html](report.html) or the
[assignment figure](../audit_figures/junction_assignment.pdf).

| Assignment outcome | Read-junction observations |
|---|---:|
| One gene and one full exon pair: retained | 152,907 |
| One gene, several possible exon pairs: excluded | 409,613 |
| Several possible genes: excluded | 3,352 |
| Several genes and exon pairs: excluded | 19,634 |
| No annotated target: excluded | 1,583 |
| Total after upstream filters | 587,089 |

All original filtering/assignment counters agree exactly with the saved
converter QC. Of 432,599 ambiguous observations, **94.69% identify one gene
but not one full exon pair**. The largest loss therefore concerns transcript
variants sharing splice boundaries, rather than mostly cross-gene assignment.
Gene identity alone is insufficient to assign the current full-interval exon
nodes, so none of these observations were added back into a graph.

The hypothetical same-gene upper bound of 562,520 observations is a diagnostic
sum, not a rescued count table. Changing to boundary nodes or choosing a
canonical exon would change the representation and requires a separate
method decision and validation. The frozen production policy remains intact.

## Strand diagnostic

Among primary, uniquely mapped spliced alignments passing the anchor and XS
checks, R1 has 264,362 antisense and 741 sense alignments; R2 has 278,423 sense
and 776 antisense alignments. Each alignment is counted once here, even if it
contains several junctions. These are approximately 99.72% RF-consistent
orientations for each mate. Confirm the library protocol and other samples
before freezing the production setting. See the
[source review](../../../../docs/dataset_review.md).

## Reproduce

From the repository root, with samtools available:

```sh
python3 preprocessing/ambiguity_audit.py \
  --bam preprocessing/runs/SRR36710387_pilot_nov17/alignment/aligned.bam \
  --gene-map references/hisat2_v50_pilot_nov17/junction_exon_map.tsv \
  --candidates data/gene_sets/candidate_genes.tsv \
  --sample-id SRR36710387 --min-anchor 8 \
  --expected-qc preprocessing/runs/SRR36710387_pilot_nov17/alignment/junction_qc.json \
  --output-dir preprocessing/runs/assignment_review_new
```

Choose a new output directory each time. The complete per-junction diagnostic
is about 27 MB and remains in the ignored run folder. This small review bundle
contains the summary, top 50 ambiguous keys, all 661 candidates' coverage,
HTML report and original audit manifest. The original manifest also describes
large local files that are intentionally not bundled. Per-gene *potential*
cross-gene support is non-additive and must not be summed across genes.

The [figure manifest](../audit_figures/figure_manifest.json) records source
table/script hashes and the SVG, PDF and PNG outputs. Figures are technical
diagnostics for one control prefix, with no disease comparison.
