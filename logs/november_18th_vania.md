# November 18, 2026 — Vania

**Status:** Pushed

**Commits:**

- `56c7330` — Priority 1 (D5, D9): node-mapping rule
- `8ed161c` — Priority 1 (D1, D9): library strandedness
- `787f8c0` — Junction threshold proposal (needs approval)
- `a0dc545` — Very minor README.md fixes

## Contribution

I settled two open input questions for the real-data pipeline (how junctions
map to graph nodes, and the library strandedness) and drafted the junction
threshold proposal. All three were done before any case/control result
existed.

### Node rule: exons that share a splice site are one node (D5, D9)

- **Problem.** Under the old rule (one node per exon interval), 74% of the
  pilot's junction observations and 30% of annotated candidate junctions
  matched several exon pairs. These differed only at exon ends a junction read
  cannot see, so the reads were discarded as ambiguous.
- **New rule** (`shared-splice-site`, now the default in
  `preprocessing/build_junction_gene_map.py`):
  - Within a gene and strand, exons that share a start or end coordinate are
    one node.
  - A retained-intron exon may not merge two nodes that each hold a normal
    exon.
  - Node IDs are `chrom:lowest-highest:strand`.
  - Junctions joining the same node pair are summed into one edge; junctions
    with both ends in one node are excluded and counted.
- **Effect:**

  | Rule | Candidate junctions mapped uniquely | Candidate genes with a cycle | Median cycle rank |
  |---|---|---|---|
  | exon-interval (old) | 69.6% | 298 | 0 |
  | shared-splice-site (new) | 98.5% | 536 | 4 |

  Same-gene candidate reads kept in the pilot rise from about 35% to about
  100%. Splice-site hints are unchanged, so the existing HISAT2 index stays
  valid; on GENCODE v50 the map has 536,254 nodes. The old rule stays
  available to reproduce the November 17 pilot.
- Updated the read counters, the ambiguity audit and the pipeline options;
  added `scripts/node_rule_review.py` and 6 tests. The decision is recorded in
  `docs/decisions.md` and `docs/data_methods.md`.

### Strandedness: all eight runs are RF (D1, D9)

- The written sources do not state the library orientation, so I checked the
  reads directly, without an aligner. I streamed the first 200,000 reads of
  each mate of every run from ENA and compared their 25-mers with both strands
  of nine RefSeq mRNAs of highly expressed genes
  (`preprocessing/strand_kmer_check.py`, with tests).
- Read 1 was antisense in 98.1–99.5% of matched reads in every run, in both
  groups: an RF (dUTP-style) library. The pilot alignment had independently
  given 99.7% for one run.
- The cohort alignment uses `--strandedness RF`; a sample below 95% read-1
  agreement after alignment is held for review. Results are in
  `data/strandedness_2026-11-18/`.

### Junction threshold proposal

Wrote `docs/threshold_proposal.md` for the six settings still open:

- anchor 8 bases;
- the read counter writes every junction (support 1);
- one shared edge set per gene: an edge needs ≥ 5 reads in ≥ 4 of 8 samples,
  chosen without the group labels;
- a sample needs ≥ 50 reads on a gene's shared graph;
- a gene needs cycle rank ≥ 1;
- a shared edge with no reads in a covered sample is a measured 0, and an
  uncovered sample is missing.

The shared edge set prevents a technical confound: the runs range from 49.5
to 80.0 million read pairs, and deeper samples would otherwise find more rare
junctions and more cycles. The proposal also declares a sensitivity grid of
six one-at-a-time runs.

## Current scope

The node rule and strandedness are decided. The thresholds were a proposal
awaiting sign-off.
