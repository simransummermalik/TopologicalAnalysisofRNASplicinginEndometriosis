# Splice Girl project priorities

Updated October 22, 2026.

This is the short, ordered version of the full
[task guide](delegations_and_tasks/tasks.md). It covers the entire project.
Anyone can pick up any available task; there are no assigned leads.

## What matters most

Splice Girl already has working Hodge math, a Rust CLI, a synthetic demo,
comparison outputs, and a basic graph figure. The biggest remaining goal is to
connect trustworthy biological inputs to that working software and produce a
careful case-versus-control result.

## Priority 1 — Finish the real input pipeline

**Task IDs:** D1, D3, and D5--D12  
**Status:** Partly complete; this is the main blocker.

- [ ] Finish the dataset suitability review and raw-file checksum record.
- [ ] Build the HISAT2 index and splice-site hints from the verified reference.
- [ ] Run read QC and Cleaver cleaning for every sample.
- [ ] Align the cleaned paired reads with the frozen settings.
- [ ] Extract uniquely mapped splice-junction counts.
- [ ] Apply the agreed support, ambiguity, and coverage checks.
- [ ] Convert every eligible sample into Splice Girl's standard junction table.
- [ ] Confirm that the same gene uses comparable nodes across samples.
- [ ] Freeze a small integration subset before attempting a larger run.

**Finished when:** every eligible sample has a checked, standardized junction
table that Splice Girl can read without manual editing.

## Priority 2 — Finish the candidate-gene evidence

**Task IDs:** G1--G10 and E1--E11  
**Status:** Research work has not been formally completed.

- [ ] Freeze rules for genetic-risk evidence.
- [ ] Build the genetic-risk evidence table with stable gene IDs and sources.
- [ ] Freeze rules for environmental-response evidence.
- [ ] Build the environmental evidence table with stable gene IDs and sources.
- [ ] Review duplicate, outdated, and ambiguous gene names.
- [ ] Create risk-only, environment-only, and overlap gene lists.
- [ ] Record the exact source and reason for including every gene.

**Finished when:** every candidate gene can be traced to reviewed evidence and
placed into one reproducible category.

## Priority 3 — Prove one real gene works end to end

**Task IDs:** D/R6.1 and D/R6.2  
**Depends on:** Priorities 1 and 2.

- [ ] Choose one eligible gene using a written rule.
- [ ] Run its junction table through Splice Girl.
- [ ] Compare parsed counts with the source table.
- [ ] Inspect its graph and Hodge checks manually.
- [ ] Record excluded samples, junctions, and reasons.
- [ ] Review the pilot before expanding to more genes.

**Finished when:** one real gene moves from checked junction counts to a saved,
reproducible Splice Girl result.

## Priority 4 — Freeze and validate the statistics

**Task IDs:** S1--S8 and R11  
**Status:** Not started.

- [ ] Confirm the biological unit and whether samples are independent or paired.
- [ ] Freeze eligible samples, genes, thresholds, and missing-value rules.
- [ ] Define the main case-minus-control effect size.
- [ ] Specify the permutation test, random seed, and when inference is allowed.
- [ ] Test the calculations on tiny known examples.
- [ ] Define the multiple-testing correction and its gene family.
- [ ] Connect approved statistics fields to the saved CLI output.

**Finished when:** another person can reproduce the toy answers and explain
every statistics column without guessing.

## Priority 5 — Complete the case/control analysis

**Task IDs:** C1--C4, S9--S12, and R12  
**Depends on:** Priorities 3 and 4.

- [ ] Create one sample-by-gene result table.
- [ ] Check group labels, graph comparability, missingness, and exclusions.
- [ ] Calculate effect sizes and valid statistical results.
- [ ] Compare Hodge results with simpler baseline measurements.
- [ ] Run the approved sensitivity checks.
- [ ] Rank results with deterministic tie-breaking.
- [ ] Export statistics, baselines, rankings, and junction explanations.

**Finished when:** every reported difference is traceable to individual samples,
an approved analysis rule, and a versioned output file.

## Priority 6 — Finish and harden the software

**Task IDs:** R14--R17 and U1, U3, and U5  
**Status:** Several pieces are already in progress and can be done without real
data.

- [x] Finish tests for missing files, empty files, bad columns, invalid values,
  and unsafe output replacement.
- [x] Make every failure explain what happened and how to fix it.
- [ ] Record a runtime and peak-memory benchmark.
- [ ] Record the minimum supported Rust version.
- [ ] Test `cargo install --path .` from a clean checkout.
- [ ] Test `splice-girl --help`, `splice-girl --version`, and
  `splice-girl --start` after installation.
- [ ] Ask a beginner to run a clean copy and fix anything that confuses them.

**Finished when:** a new user can install, run, and troubleshoot the tool from
the written instructions alone.

## Priority 7 — Complete the figures and demo

**Task IDs:** V3--V11  
**Status:** The tiny demo and three-part Hodge SVG work.

- [x] Show the original, gradient, and cycle signals with identical graph layouts.
- [ ] Create one real-gene case/control figure after the pilot is approved.
- [ ] Create result figures only from versioned tables.
- [ ] Explain junctions, cycle fraction, and limitations in plain language.
- [ ] Choose the real demo gene using a written rule.
- [ ] Write speaking parts and handoff sentences.
- [ ] Rehearse with a timer and save an offline screenshot or recording backup.
- [ ] Audit every label and caption against its source output.

**Finished when:** the team can explain the whole project clearly and complete
the demo even without internet access.

## Priority 8 — Interpret the results carefully

**Task IDs:** I1--I3  
**Depends on:** Priorities 2 and 5.

- [ ] Join computational results to the evidence tables by stable gene ID.
- [ ] Compare risk-only, environment-only, and overlap categories.
- [ ] Link every discussed gene to its samples, junctions, and sources.
- [ ] State that cycle-space signal does not by itself prove circular RNA.
- [ ] State tissue, sample-size, alignment, coverage, and causation limits.
- [ ] Report a null result honestly if Hodge summaries do not beat baselines.

**Finished when:** every biological claim has visible computational evidence,
source provenance, and a matching limitation.

## Priority 9 — Assemble the final project

**Task IDs:** F1--F5

- [ ] Assemble the final methods from the frozen specifications.
- [ ] Assemble the final results from versioned outputs.
- [ ] Prepare the software release, license, version, and release notes.
- [ ] Compile and inspect the LaTeX documents.
- [ ] Run the complete acceptance checklist from a clean checkout.

**Finished when:** the report, presentation, software, figures, and reproducible
results agree with one another.

## Priority 10 — Keep the project organized throughout

**Task IDs:** P1.2, P1.3, and P1.5--P1.7

- [ ] Record course deadlines and internal checkpoints.
- [ ] Record where the main Overleaf source lives.
- [ ] Create `docs/decisions.md` and add dated decisions without deleting history.
- [ ] Add dependencies and target dates to the shared task board.
- [ ] Keep file-naming, large-file, and output-location rules in the decisions log.

**Finished when:** nobody needs chat history to understand a deadline, decision,
file location, or dependency.

## Four useful tasks that can happen at the same time

These are choices for four people, not assigned roles:

1. Continue the reference/index and standardized-junction work in Priority 1.
2. Continue the two evidence tables in Priority 2.
3. Write and test the statistical specification in Priority 4.
4. Finish the error audit and three-part graph view in Priorities 6 and 7.

## Already finished — do not redo

- [x] The core Hodge decomposition equations and mathematical tests.
- [x] The Rust parser, graph construction, normalization, and decomposition.
- [x] The one-sample and manifest-driven CLI workflow.
- [x] Saved result bundles and reproducibility manifests.
- [x] The synthetic acceptance tests and tiny demo.
- [x] A generated three-panel observed/gradient/cycle SVG view.
- [x] The beginner README, input contract, and math specification.
- [x] The eight-sample manifest and verified human reference pair.
- [x] Frozen exon-node, coordinate, strand, and annotation-mapping rules.

## Later, only if the required project is complete

- [ ] Extra website styling.
- [ ] Optional extensions X1--X4.
- [ ] Large-scale optimization beyond what the selected analysis needs.
- [ ] Extra biological claims that are not supported by the final results.
