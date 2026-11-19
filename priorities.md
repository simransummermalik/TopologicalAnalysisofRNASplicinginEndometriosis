# Splice Girl project priorities

Updated November 19, 2026.

This is the short, ordered version of the full
[task guide](delegations_and_tasks/tasks.md). It covers the entire project.
Anyone can pick up any available task; there are no assigned leads.

## Latest completed work — November 19

- [x] Prepare a local 0.2.0 release draft under CC BY 4.0, following the
  project's advice to allow commercial reuse with attribution. The published
  0.1.0 crate remains MIT-licensed. All four named contributors approved the
  0.2.0 license update on November 19; the registry publication is in progress.

- [x] Run the current packaged tool on all 243 candidate genes observed in
  the real control prefix; independently recount all 1,388 edges from BAM
  and check every gene against SVD. Six graphs have cycles and five have
  nonzero cycle components. This uses the historical pilot protocol, not
  the final cohort protocol. See the [real-data validation](data/pilots/SRR36710387_real_validation_nov18/README.md).
- [x] Add saved JSON configuration, read-only preflight and configured reports,
  including completed sample/batch/reference/map identity and count checks.
- [x] Publish 0.1.0 under MIT with a separate citation request, provide
  CITATION.cff and `cite`, and restrict the crate to software, fixtures and
  relevant guides. A local 0.2.0 CC BY 4.0 license update is now drafted.
- [x] Build, inspect, test and install the actual unpublished crate archive;
  exercise the installed guide outside the checkout. Independent human
  acceptance remains pending in [the reviewer exercise](docs/user_acceptance.md).
- [x] Make shared preparation and sensitivity available in the installed Rust
  CLI, including embedded example inputs, explicit statistics, provenance and
  binary-only checks with PATH empty. No Python or checkout is required for
  this junction-table analysis workflow; raw-read preprocessing stays separate.
- [x] Implement the [proposed shared-graph workflow](docs/shared_graphs.md):
  explicit edge/coverage thresholds, measured zeros only for covered samples,
  cycle eligibility, full-cohort gene-family handoff, and seven one-factor
  sensitivity runs with optional Rust statistics. An eight-sample synthetic
  example checks the whole path; biological threshold approval remains open.
- [x] Settle library strandedness: all eight runs are RF (read 1 antisense in
  98.1–99.5% of matched reads, both groups), checked by streaming 200,000
  reads per mate from ENA. The cohort alignment uses `--strandedness RF`. See
  the [check](data/strandedness_2026-11-18/README.md).
- [x] Settle the junction-ambiguity rule: annotated exons sharing a splice
  site are one node (`shared-splice-site`), with retained-intron exons kept
  from merging neighbors. 98.5% of candidate junctions now map uniquely
  (was 69.6%); the pilot's same-gene candidate reads kept rise from 35% to
  about 100%. Rebuild the exon map (seconds) before the cohort run; the HISAT2
  index stays valid. See the [decision](docs/decisions.md) and
  [review](docs/data_methods.md#node-rule-review-2026-11-18).
- [x] Add optional streaming BAM to single-sample and batch preprocessing,
  keeping the existing retained-SAM mode and unchanged count filters.
- [x] Check estimated output/scratch space before starting and recheck a
  configurable reserve between stages.
- [x] Supervise both sides of each pipe, stop children on failure/interruption,
  and retain stage logs plus incomplete outputs without publishing a complete
  cohort bundle.
- [x] Compare actual tools in both modes on synthetic paired reads; counts/QC
  agree. Force real sort spill with 12,000 synthetic pairs per sample.
- [x] Pass the complete project check: 100 Rust and 100 Python tests, strict
  Clippy, formatting, release build and actual-tool checks in both modes.

No human-read run was started. See the [storage guide](preprocessing/HISAT2_GUIDE.md#save-disk-space-with-streaming-bam)
and [validation evidence](preprocessing/validation/2026-11-18-storage-mode/README.md).

## Previously completed — November 17

- [x] Diagnose the real pilot's assignment loss with a tested
  [exon-pair/gene ambiguity and strand audit](data/pilots/SRR36710387_prefix1m_nov17/assignment_audit/README.md).
- [x] Benchmark Rust analysis/export on all 243 candidate genes observed in
  that real prefix; every group passes the internal math checks.
- [x] Export and inspect three plain scientific pilot figures as PDF/SVG/PNG.
- [x] Write a shared speaking plan with live commands, explanations and handoffs.
- [x] Assemble methods/results drafts with technical findings separated from
  biological comparisons that have not run.
- [x] Correct the BH/permutation-resolution explanation and add mathematical
  regression examples without changing the frozen four-versus-four policy.
- [x] Compile and inspect the available LaTeX roadmap section using a temporary
  standard wrapper; preserve the simple academic appearance.
- [x] Prepare an unpublished source-package workflow and release-review notes.
- [x] Pass the full offline check from an extracted source snapshot, verify
  its archive/file hashes, and run the isolated installed CLI and live demo.
- [x] Complete and document an [assistant-run beginner walkthrough](beginner_test/README.md):
  25 CLI checks, the 100-test Rust suite, copied outputs and a browser screenshot.
  Independent first-time human feedback remains open.
- [x] Pass 100 Rust and 88 Python tests, strict Clippy, formatting and release
  build; Rust tests/Clippy also pass on the minimum 1.85.0 compiler.

**Actual remaining dependencies:** the user confirmed there is no external
drive or server available. The laptop has about 9 GB free, and seven samples'
read pairs are not available locally. The full cohort has not been launched.
The [November 17 storage preflight](docs/cohort_execution.md) estimated about
768 GB for the retained-SAM mode keeping all raw/cleaned reads, SAMs and BAMs.
Streaming BAM now avoids SAM, but does not establish that all samples fit
this laptop. Final preprocessing choices, independent
evidence/beginner reviews, the Overleaf parent source, course dates and the
final release review also remain. These are not marked complete by software
tests or by writing a draft.

## What matters most

Splice Girl already has working Hodge math, a Rust CLI, a synthetic demo,
comparison outputs, and a basic graph figure. The biggest remaining goal is to
connect trustworthy biological inputs to that working software and produce a
careful case-versus-control result.

## Priority 1 — Finish the real input pipeline

**Task IDs:** D1, D3, and D5--D12  
**Status:** Partly complete; this is the main blocker.

- [ ] Resolve the [documented dataset review](docs/dataset_review.md), including
  conflicting library/read-length descriptions, and complete all raw checksums.
- [x] Build the full human HISAT2 index, exon map and splice-site hints from
  the verified reference (November 17).
- [x] Implement and test the Cleaver/HISAT2 adapter matching the current data
  methods, including a sequential batch command. The legacy STAR prototype
  remains separate. Full-cohort use still needs the checked human pilot.
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
**Status:** In progress. The rules, verified sources, dioxin-like chemical
list, evidence templates, and the tested `gene_sets/gene_sets.py` tooling
exist; see [`docs/gene_set_methods.md`](docs/gene_set_methods.md). The
environmental (CTD) evidence table is extracted, context-reviewed (E6),
fully mapped, deduplicated (E7), and audited (E8). The genetic-risk set is
extracted, reviewed, deduplicated, audited, and exported (G4–G10; 103 genes).
Overlap and categories are built (E9–E10) under the final publication rule
(decided 2026-11-17): 96 risk-only, 558 environment-only, 7 overlap, 661
candidates. Remaining: team sign-offs (G1 cold read, independent spot checks in
[`docs/spot_checks.md`](docs/spot_checks.md)).

- [ ] Freeze rules for genetic-risk evidence. *(Drafted and recorded in
  `docs/decisions.md`; waiting for a teammate's cold read, G1 step 7.)*
- [x] Build the genetic-risk evidence table with stable gene IDs and sources.
  *(Audited; drafted before the G1 cold read, so rerun if that changes a rule.)*
- [x] Freeze rules for environmental-response evidence.
- [x] Build the environmental evidence table with stable gene IDs and sources.
  *(Audited; final publication rule decided 2026-11-17.)*
- [x] Review duplicate, outdated, and ambiguous gene names.
  *(Environment E5 and risk G7 reviews recorded in `docs/gene_set_methods.md`.)*
- [x] Create risk-only, environment-only, and overlap gene lists.
  *(`candidate_genes.tsv`: 661 genes, 7 overlap.)*
- [x] Record the exact source and reason for including every gene.
  *(Every membership row lists its evidence IDs; every evidence row names its source record and any review reason.)*

**Finished when:** every candidate gene can be traced to reviewed evidence and
placed into one reproducible category.

## Priority 3 — Prove one real gene works end to end

**Task IDs:** D/R6.1 and D/R6.2  
**Depends on:** Priorities 1 and 2.

**Status:** The [bounded CDH11 technical pilot](data/pilots/SRR36710387_prefix1m_nov17/README.md)
passes on one million cleaned read pairs from one control. The full sample,
case/control comparison and team review remain; mapping ambiguity needs review.

- [x] Choose one pilot gene using a written coverage/size rule.
- [x] Run its junction table through Splice Girl.
- [x] Compare every count with an independent BAM recount.
- [x] Inspect its graph and compare Rust with an independent SVD projection.
- [x] Record the prefix scope, junction exclusions and sensitivity results.
- [ ] Review the pilot before expanding to more genes.

**Finished when:** one real gene moves from checked junction counts to a saved,
reproducible Splice Girl result.

## Priority 4 — Freeze and validate the statistics

**Task IDs:** S1--S8 and R11  
**Status:** In progress: `stats` connects exact permutation tests, BH
correction and a gene-set test to the CLI, saved tables, provenance, and
offline report. The inference rules are frozen (2026-11-17): two-sided, 4 per
group, one primary gene-set test on the candidate set; see
[statistics_spec.md](docs/statistics_spec.md#inference-rules-decided-2026-11-17).
Real-data eligibility and thresholds remain open.

- [x] Record the biological unit and unpaired cohort design in the statistics specification.
- [ ] Freeze eligible samples, genes, thresholds, and missing-value rules.
- [x] Define the main case-minus-control effect size.
- [x] Implement the exact permutation calculation; exhaustive enumeration uses no random seed.
- [x] Test permutation and BH correction calculations on tiny known examples.
- [x] Freeze when inference is allowed and the final test direction.
- [x] Define the multiple-testing correction and its gene family.
- [x] Connect explicitly configured statistics fields to saved CLI output and the report.

**Finished when:** another person can reproduce the toy answers and explain
every statistics column without guessing.

## Priority 5 — Complete the case/control analysis

**Task IDs:** C1--C4, S9--S12, and R12  
**Depends on:** Priorities 3 and 4.

- [ ] Create one sample-by-gene result table.
- [ ] Check group labels, graph comparability, missingness, and exclusions.
- [ ] Calculate effect sizes and valid statistical results.
- [x] Implement same-sample baseline exports alongside Hodge results.
- [x] Implement an explicit support/coverage sensitivity-grid runner.
- [x] Implement optional ranking with deterministic tie-breaking.
- [x] Export configured statistics, baselines, optional rankings, and signed
  junction components, retaining the original evidence-annotated gene family.
- [ ] Apply these tools to reviewed real inputs, compare methods, and interpret
  the approved sensitivity grid and rankings.

**Finished when:** every reported difference is traceable to individual samples,
an approved analysis rule, and a versioned output file.

## Priority 6 — Finish and harden the software

**Task IDs:** R14--R17 and U1, U3, and U5  
**Status:** The developer-run installation, compatibility, candidate-sized
synthetic benchmark, real-prefix observed-candidate benchmark, CI checks, and
CLI checks and the documented assistant walkthrough are complete; an independent
beginner review remains. Full-cohort
performance is not measured.

- [x] Finish tests for missing files, empty files, bad columns, invalid values,
  and unsafe output replacement.
- [x] Make every failure explain what happened and how to fix it.
- [x] Record runtime and peak-memory checks, including a repeatable
  candidate-sized synthetic benchmark and the real-prefix benchmark.
- [x] Record the minimum supported Rust version.
- [x] Test `cargo install --path .` from an isolated clean source copy.
- [x] Test `splice-girl --help`, `splice-girl --version`, and
  `splice-girl --start` after installation.
- [x] Provide one complete software check command, including an optional actual
  Cleaver/HISAT2 paired-read batch-to-report smoke test.
- [x] Run the beginner guide in a fresh local source copy, save each command's
  evidence and correct instruction gaps; see `beginner_test/RESULTS.md`.
- [ ] Ask a beginner to run a clean copy and fix anything that confuses them.

**Finished when:** a new user can install, run, and troubleshoot the tool from
the written instructions alone.

## Priority 7 — Complete the figures and demo

**Task IDs:** V3--V11  
**Status:** The tiny demo, three-part Hodge SVG, offline comparison report,
real-prefix technical figures and speaking draft are available. The final
real case/control figure and team rehearsal remain.

- [x] Show the original, gradient, and cycle signals with identical graph layouts.
- [x] Generate a plain offline comparison report with each sample's value,
  exclusions, provenance links, and Hodge graph figures from the same run.
- [ ] Create one real-gene case/control figure after the pilot is approved.
- [x] Create the available technical figures from saved, checksummed tables.
  Final biological result figures remain dependent on Priority 5.
- [x] Draft plain-language explanations of junctions, cycle fraction and limits
  in [presentation_notes.md](docs/presentation_notes.md); audience review remains.
- [x] Choose CDH11 for the bounded technical demonstration using a written rule.
  Final cohort-demo selection remains subject to real QC/results review.
- [x] Write shared speaking parts, live commands and handoff sentences.
- [ ] Rehearse with a timer and save an offline screenshot or recording backup.
- [x] Audit the current pilot figures' labels/captions against source outputs.
  Repeat for final case/control figures when those exist.

**Finished when:** the team can explain the whole project clearly and complete
the demo even without internet access.

## Priority 8 — Interpret the results carefully

**Task IDs:** I1--I3  
**Depends on:** Priorities 2 and 5.

- [ ] Join computational results to the evidence tables by stable gene ID.
- [ ] Compare risk-only, environment-only, and overlap categories.
- [ ] Link every discussed gene to its samples, junctions, and sources.
- [x] State that cycle-space signal does not by itself prove circular RNA.
- [x] Document current tissue, sample-size, alignment, coverage and causation
  limits in the methods/results drafts and speaking notes; apply them to final
  interpretations after the cohort run.
- [ ] Report a null result honestly if Hodge summaries do not beat baselines.

**Finished when:** every biological claim has visible computational evidence,
source provenance, and a matching limitation.

## Priority 9 — Assemble the final project

**Task IDs:** F1--F5

- [ ] Finalize the [assembled methods draft](docs/methods_draft.md) after real
  eligibility and preprocessing choices are resolved.
- [ ] Finalize the [technical results draft](docs/results_draft.md) after the
  full cohort analysis; do not substitute pilot measurements for disease results.
- [ ] Finalize the [prepared release review](docs/release_notes.md) after
  independent feedback. Version 0.1.0 is published on crates.io; version
  0.2.0 is prepared under CC BY 4.0 with contributor authorization, and is
  being published now.
- [ ] Compile the missing Overleaf parent. The available roadmap section is
  compiled, inspected and saved in `roadmap-preview.pdf`.
- [ ] Run the complete acceptance checklist from a clean checkout.

**Finished when:** the report, presentation, software, figures, and reproducible
results agree with one another.

## Priority 10 — Keep the project organized throughout

**Task IDs:** P1.2, P1.3, and P1.5--P1.7

- [ ] Record course deadlines and internal checkpoints.
- [ ] Record where the main Overleaf source lives. The source map now records
  that the parent `.tex` file is still in Overleaf and not in this checkout.
- [x] Create `docs/decisions.md` and add dated decisions without deleting history.
- [ ] Add dependencies and target dates to the shared task board.
- [x] Record file-naming, large-file and output-location rules in the decisions
  log and release checklist; maintain them as new outputs are added.

**Finished when:** nobody needs chat history to understand a deadline, decision,
file location, or dependency.

## Four useful tasks that can happen at the same time

These are choices for four people, not assigned roles:

1. Review the source-method discrepancy, pilot assignment loss and final data
   settings; the index and diagnostic software are already built.
2. Complete the independent evidence cold read/spot checks in Priority 2.
3. Have a beginner run the written demo and record where they get stuck.
4. Rehearse the shared speaking plan and make the offline recording backup.

## Already finished — do not redo

- [x] The core Hodge decomposition equations and mathematical tests.
- [x] The Rust parser, graph construction, normalization, and decomposition.
- [x] The one-sample and manifest-driven CLI workflow.
- [x] Saved result bundles and reproducibility manifests.
- [x] Compact input/manifest validation through the real graph and math path.
- [x] Sample-preserving descriptive comparison across every observed gene.
- [x] The synthetic acceptance tests and tiny demo.
- [x] A generated three-panel observed/gradient/cycle SVG view.
- [x] The beginner README, input contract, and math specification.
- [x] The eight-sample manifest and verified human reference pair.
- [x] Frozen exon-node, coordinate, strand, and annotation-mapping rules.
- [x] Documented that the current HISAT2 annotated-junction workflow cannot
  establish back-splice absence or circular RNA.

## Later, only if the required project is complete

- [ ] Extra website styling.
- [ ] Optional extensions X1--X4.
- [ ] Large-scale optimization beyond what the selected analysis needs.
- [ ] Extra biological claims that are not supported by the final results.
