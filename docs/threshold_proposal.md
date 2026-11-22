# Proposal: junction support and coverage thresholds

**Status:** provisionally approved without changes on 2026-11-22 by Vania
Zeledon; confirmation by the rest of the team is pending (drafted 2026-11-18).
Recorded in [`decisions.md`](decisions.md). These are the primary settings for
the eight-sample run; if the team changes any value, record a new dated
decision before any case/control result is seen. The
[shared-graph implementation](shared_graphs.md) takes these values as explicit
settings. Complete inputs, storage/compute and final QC review are still needed
for the full eight-sample run.

**Must be decided before any case/control result is seen.** Every rule below
uses read counts only, never the case/control labels.

## What needs deciding

| # | Setting | Where it applies | Proposed value |
|---|---|---|---|
| 1 | Minimum anchor length | Read counter, at alignment | **8 bases** on each side of the junction |
| 2 | Minimum reads to write a junction | Read counter, at alignment | **1** (keep everything; filter later) |
| 3 | Which edges form a gene's graph | One shared edge set per gene, all samples | Edge has **≥ 5 reads in at least 4 of the 8 samples** |
| 4 | Sample coverage for a gene | Per sample and gene | **≥ 50 reads** in total on the gene's shared graph |
| 5 | Gene eligibility | Per gene | Shared graph has **at least one cycle** (cycle rank ≥ 1) |
| 6 | Zero versus missing | Per sample and gene | Shared edge with no reads in a covered sample = **measured 0**; sample below coverage = **missing** |

Already decided elsewhere: unique alignments only (`NH:i:1`), RF strandedness,
the shared-splice-site node rule, PCR duplicates kept, and a gene is tested
only with all 4 cases and 4 controls eligible.

## Why one shared edge set per gene (rule 3)

Cycle fraction depends on which edges exist. If each sample kept only its own
well-supported edges, a deeper sample would see more rare junctions, form more
cycles and get a different cycle fraction for technical reasons alone. The runs
range from 49.5 to 80.0 million read pairs, so this would be a real confound.

With one shared edge set, every sample of a gene has the same graph and only
the read proportions on its edges differ. A difference in cycle fraction then
reflects signal on a common topology. The edge set is chosen from all eight
samples together without the labels, so selection is unchanged when those
labels are permuted. Sequencing depth and unequal detection can still affect
the counts and require QC and sensitivity review.

"At least 5 reads in at least 4 of 8 samples" keeps an edge that is
reproducibly present in half the cohort. Five reads is enough that an edge
is unlikely to be a single misaligned read; requiring 4 samples stops one
unusual sample from adding edges for everyone.

## Why these values for the other rules

1. **Anchor 8.** A read must align at least 8 bases on each side of the
   junction. Only annotated junctions are counted, so a short anchor is
   unlikely to create a false junction; 8 is a common default (used in the
   pilot) and drops reads whose junction placement is uncertain.
2. **Write everything (1 read).** Filtering happens at rule 3 instead, so the
   sensitivity grid below can be rerun without realigning.
4. **50 reads per sample.** Cycle fraction is computed on read proportions.
   With fewer reads, one or two reads move the proportions a lot. Projected
   from the pilot's shallowest sample, about 365 of 661 candidate genes reach
   50 same-gene junction reads at full depth (401 at 20, 345 at 100). This is
   a rough estimate from one control at 2% depth.
5. **Cycle rank ≥ 1.** A graph with no cycle has cycle fraction exactly 0 in
   every sample. Such a gene carries no information and would only dilute the
   gene-set test. It is reported as "no cycle in shared graph", not tested.
6. **Zero versus missing.** This follows the project-wide rule: a measured zero
   is data, a missing measurement is never filled with zero.

## Sensitivity plan (declared now, reported whatever it shows)

The primary result uses the values above. The same analysis is repeated at:

- edge reads: 3, 5, 10 (samples fixed at 4 of 8);
- samples per edge: 2, 4, 6 of 8 (reads fixed at 5);
- sample coverage: 20, 50, 100 reads.

Each setting is varied one at a time, giving 6 extra runs. All are reported,
including settings where a gene drops out. A result that appears only at one
setting is described as such.

## Implementation available for review

- `scripts/shared_graph.py` builds each gene's shared edge set, writes covered
  samples' counts including measured zeros, and records sample/gene exclusions.
- The eligible statistics family requires cycle rank ≥ 1 and coverage in every
  manifest sample. Its original gene annotations are retained.
- `scripts/sensitivity.py --mode shared` runs the seven one-factor scenarios
  above, with optional explicitly configured Rust statistics. The earlier
  per-sample mode remains available to reproduce previous results.

The [synthetic walkthrough](../examples/shared_graph/README.md) exercises four
cases and four controls without human reads. A zero-count shared edge remains
part of the mathematical graph; a nonzero cycle component does not prove
observed use of that edge in that sample. See [interpretation changes](shared_graphs.md#interpretation-changes-to-explain).
Anchor sensitivity requires separately recounted junction exports and is not
evaluated by the count-table grid.

## Sign-off

- **Approved by:** Vania Zeledon (provisional); other team members to confirm:
  - [ ] Sydney Holmes
  - [ ] Summer Malik
  - [ ] Wynn Oo
- **Date:** 2026-11-22
- **Changes from this proposal, if any:** none
