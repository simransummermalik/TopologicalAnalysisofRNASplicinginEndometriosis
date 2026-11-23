# Primary cohort analysis settings

[`analysis.json`](analysis.json) is the saved configuration for the study's
primary analysis of GSE315857 (4 endometriosis cases, 4 controls). Every
value in it was decided before this file was written; this file only puts
those decisions in one place. Prepared 2026-11-23.

**Do not change a value here to get a different result.** A change needs a
new dated entry in [`docs/decisions.md`](../../docs/decisions.md) first, and
because a one-case-versus-one-control comparison has already been viewed, the
entry must say the change was made after case/control data were seen.

## Where each value comes from

| Setting | Value | Decision |
|---|---|---|
| Samples and labels | the 8 runs in `data/metadata/samples.tsv` | 2026-10-20 cohort; 2026-10-29 one sample per subject |
| Gene family | `data/gene_sets/candidate_genes.tsv` (661 genes), label `candidate_genes` | 2026-11-17 primary gene-set test |
| Shared edge set | edge ≥ 5 reads in ≥ 4 of 8 samples | [threshold proposal](../../docs/threshold_proposal.md), provisional 2026-11-22 |
| Sample coverage | ≥ 50 reads on the gene's shared graph | same |
| Gene eligibility | cycle rank ≥ 1 and covered in every sample (applied by the software) | same |
| Sensitivity grid | edge reads 3/5/10, samples per edge 2/4/6, coverage 20/50/100, one at a time | same |
| Test | two-sided, independent samples, all 4 cases and 4 controls required | 2026-11-17 inference rules |
| Primary statistic | mean over genes of the absolute case-minus-control cycle-fraction effect, 70 relabelings, α = 0.05 | 2026-11-17 inference rules |
| Node rule | `shared-splice-site` | 2026-11-18 |
| Strandedness | RF; samples below 95% read-1 agreement excluded | 2026-11-18 |
| Anchor, read-counter support, duplicates | 8 bases, 1 read, duplicates kept | threshold proposal; data methods |

Only the gene-set test is confirmatory. Gene-level results, the risk,
environment and overlap sets, and every sensitivity scenario are exploratory.

## Running it on the cluster

All paths are relative to this folder. The cohort and reference files stay on
the cluster and are ignored by Git.

1. Pull the repository on the cluster and build or install `splice-girl`.
2. Link the reference folder (the one holding `reference_manifest.json` and
   `junction_exon_map.tsv.manifest.json`) here:
   ```sh
   ln -s /path/to/reference analysis/cohort_primary/reference
   ```
3. After all eight runs finish and pass the strand check, combine them into
   this folder (see the
   [HISAT2 guide](../../preprocessing/HISAT2_GUIDE.md#combine-samples-that-were-run-as-separate-jobs)):
   ```sh
   python3 preprocessing/assemble_runs.py --samples runs.tsv \
     --output-dir analysis/cohort_primary/cohort \
     --strandedness RF --min-anchor 8 --min-unique 1
   ```
4. Check everything without writing results:
   ```sh
   splice-girl check --config analysis/cohort_primary/analysis.json
   ```
5. Run it:
   ```sh
   splice-girl run --config analysis/cohort_primary/analysis.json
   ```
   The primary gene-set result is
   `results/analysis/primary/statistics/set_test.tsv`; `results/report.html`
   shows everything.

## If the check fails

Fix the cause, do not edit the settings. In particular, if any sample was
excluded by the strand check, the check stops because the sample sheet lists
eight samples and the cohort has fewer. That is intended: with fewer than 4
per group the primary test cannot run, and how to proceed is a team decision
to record before any result is produced.

## Tested

On 2026-11-23 the configuration passed `check` and completed `run` against a
mock cohort: eight fake runs with the real sample IDs, the real sample sheet
and candidate list, and synthetic counts. That shows the files fit together;
it says nothing about the biology.
