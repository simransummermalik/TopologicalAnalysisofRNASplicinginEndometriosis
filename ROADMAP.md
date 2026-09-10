# Project Execution Roadmap

Topological Analysis of RNA Splicing in Endometriosis: Integrating Genetic Susceptibility, Environmental Response, and Hodge Theory

| Phase | What We Are Doing | Main Output | Depends On | Priority | Status |
|-------|-------------------|-------------|------------|----------|--------|
| 1 | Organize documents, owners, and interfaces | Source map, task assignments, draft schemas | Available proposal | MUST HAVE | Ready |
| 2 | Assemble genetics and environmental gene evidence | Risk, environment, overlap, and candidate TSVs | 1 | MUST HAVE | Not Started |
| 3 | Prepare endometriosis/control junction data | Metadata, junction tables, QC report | 1; dataset decision | MUST HAVE | Not Started |
| 4 | Specify and validate graph mathematics | Math specification, synthetic fixtures, expected answers | 1; graph contract | MUST HAVE | In Progress |
| 5 | Build Rust MVP | Synthetic input-to-output CLI | 1, 4; stable schema | MUST HAVE | In Progress |
| 6 | Integrate real splice data | Validated real-data tool run | 2, 3, 5 | MUST HAVE | Not Started |
| 7 | Produce disease/control comparisons | Per-sample values and gene differences | 6 | MUST HAVE | Not Started |
| 8 | Test differences and compare baselines | Effects, permutation results, FDR, baseline report | 7; valid sample design | MUST HAVE | Not Started |
| 9 | Interpret validated findings | Evidence-linked gene summaries and limitations | 2, 8 | MUST HAVE | Not Started |
| 10 | Finish reusable tool and documentation | Tested Rust CLI release and example | 5–8; starts alongside integration | MUST HAVE | Not Started |
| 11 | Rehearse a small live demonstration | Synthetic and real-gene demo; backup | 9, 10 | MUST HAVE | Not Started |
| 12 | Assemble final report and presentation | Figures, slides, report, reproducibility record | 8–11 | MUST HAVE | Not Started |
| 13 | Pursue proposal extensions if time permits | Separately documented extension | Core final milestone | STRETCH | Not Started |

Status vocabulary: **Not Started**, **Ready**, **In Progress**, **Blocked**, **Complete**. Ready means prerequisites permit starting; it does not mean work has been performed. A basic Rust implementation and two synthetic checks now exist; the remaining synthetic validation, gene selection, dataset acquisition, real-data analysis, statistics, ranking, and final tool work have not been performed.

## Authority and scope

The source is [hodgetheorynew.pdf](assets/hodgetheorynew.pdf), the complete 24-page proposal dated September 2, 2026. Repository inspection found no `.tex` source, no existing `sections/` folder, and no existing implementation. Consequently the PDF supplies the scientific commitments; exact Overleaf preamble compatibility must be checked when the source is restored. The separate `roadmap.tex` is an include-only section for `\input{roadmap}`. It does not replace or rewrite the proposal.

The approved question is whether genes associated with endometriosis differ in RNA-splicing topology between endometriosis and normal endometrial tissue. Genetic susceptibility and environmental response prioritize genes. Graph-level Hodge decomposition separates normalized junction signals into gradient and cycle-space components. The study is the first application of a reusable tool written primarily in Rust.

The detailed file contracts, task assignments, acceptance gates, and sequencing below operationalize that proposal. They are proposed implementation arrangements, not claims that the proposal already fixed those details. Multiple-testing correction and explicit baseline deliverables implement the additional roadmap requirements. No external evidence search or new scientific project selection is part of this task.

### Proposal-to-execution traceability

| Proposal sections (PDF pages) | Commitment to carry out | Execution location |
|---|---|---|
| 1–4, 23–24 (1–3, 16–17) | Research question; risk/environment/overlap groups; hypothesis remains testable | Biology workstreams; phases 2, 9 |
| 5–7, 26 (3–5, 18) | Junction measurements, splice graphs, public tissue data, candidate filtering | Data contract; phases 3, 6 |
| 8–12, 14–17 (6–12) | Incidence matrix, gradient and cycle projections, cycle-energy fraction, comparisons | Mathematics; phases 4, 7 |
| 13, 32 (8–10, 21–22) | Optional meaningful B2, exposure data, genome-wide and other-disease extensions | Phase 13 only |
| 18–22, 31 (12–16, 21) | Rust CLI, module structure, ranking, junction attribution, small MVP | Rust workstream; phases 5, 10 |
| 25 (17–18) | Linear, skipping, back-splicing, noise, weight-sweep validation | Phase 4 acceptance fixtures |
| 27–28 (19–20) | Fast one-gene demonstration; per-sample permutation comparison | Statistics; phases 8, 11 |
| 29–30, 33–36 (20, 22–24) | Success criteria, limitations, final computational biology product | Phases 9, 12; final gate |

## Who is doing what

Fill in names and update current tasks at each team meeting. One student can hold multiple roles; for four or five students, biology/environment and statistics/presentation can share owners. Every artifact still needs one accountable owner and a reviewer from the receiving workstream.

| Person / Role | Primary responsibility | Secondary responsibility | Current task to assign | Next task | Blocked by | Deliverable |
|---|---|---|---|---|---|---|
| ___ / Coordinator | Task board, source recovery, handoffs | Release checklist | Restore Overleaf source and assign owners | Record interface decisions | Source availability for original-document compilation | `docs/source_map.md`, task board |
| ___ / Biology Lead | Endometriosis genetic evidence | Interpretation | Extract proposal-aligned inclusion criteria | Assemble evidence rows | Criteria review | `data/gene_sets/endometriosis_risk.tsv` |
| ___ / Environment Lead | Exposure-response evidence | Candidate merging | Specify exposure scope and context fields | Assemble environmental evidence | Exposure scope | `data/gene_sets/environment_response.tsv` |
| ___ / Data Lead | Dataset, metadata, junction preparation | Demo data | Assess dataset suitability and access | Prepare comparable sample manifest | Dataset decision/access | Standardized junction TSVs and QC |
| ___ / Math Lead | Projection specification and validation | Numerical review | Fix orientation and synthetic expected answers | Review solver behavior | Graph contract | `docs/math_spec.md`, fixture expectations |
| ___ / Rust Lead | CLI and reusable analysis modules | Integration and packaging | Draft input/output interfaces | Implement parser and graph construction | Schema agreement | Rust CLI and tests |
| ___ / Stats Lead | Sample comparisons, FDR, baselines | Results review | Specify output fields and sample-design checks | Implement comparison validation | Real inference requires sample output | `docs/statistics_spec.md`, validated comparison outputs |
| ___ / Presentation Lead | Figures, narrative, live demo | Usability review | Sketch pipeline and demo storyboard | Prepare example run instructions | Results needed only for final figures | Figures, slides, demo script |

## Dependencies and simultaneous work

```text
Endometriosis genetics ----\
                           +--> Candidate gene sets ------------------\
Environmental response ---/                                           |
                                                                       v
RNA-seq / metadata --> Junction extraction --> Standard tables --> Real integration
                                                |                      ^
                                                v                      |
                                         Rust graph engine ------------+
                                                |                      |
Graph contract --> Synthetic math validation --> Hodge engine ----------+
                                                                       |
                                                                       v
                                                             Per-sample results
                                                                       |
                                                                       v
                                                        Case/control comparisons
                                                                       |
                                                                       v
                                                        Statistics + baselines
                                                                       |
                                                                       v
                                                       Biological interpretation
                                                                       |
                                            Reusable tool + demo + presentation
```

The graph contract connects data identifiers, orientation, matrix indexing, and output junction names. Freeze a reviewed first version before production integration; version subsequent changes. Gene lists filter and annotate the real analysis, but do not block a generic parser or synthetic engine.

- Genetics and environmental evidence collection proceed independently after their criteria are recorded; merging depends on normalized identifiers from both.
- Dataset suitability/access work proceeds alongside synthetic mathematics and parser design.
- Rust input validation and graph construction can use planned synthetic fixtures before real data arrive. Numerical implementation depends on the reviewed math contract.
- Statistics can specify interfaces and future validation cases early. Real tests depend on validated per-sample results and a valid sample design.
- Tool documentation and the pipeline figure can start early. Final interpretation and result figures wait for QC, statistics, and baseline comparisons.
- Phase 10 usability work runs alongside phases 6–9; phase 13 must not consume time needed for the core release.

## Shared deliverable contracts

All paths below are **future deliverables**, not files created by this documentation task. Record schemas in `docs/data_contract.md`, including required columns, allowed values, missing-value rules, and schema version. TSV is the initial interface proposed here; confirm it at phase 1. Keep raw evidence separate from deduplicated gene membership.

| Future artifact | Required content | Handoff check |
|---|---|---|
| `data/gene_sets/risk_evidence.tsv` | `gene_id`, `gene_symbol`, `original_id`, `evidence_type`, `source`, `source_version`, `mapping_method`, `notes` | One evidence assertion per row; ambiguous mappings retained for review |
| `data/gene_sets/endometriosis_risk.tsv` | `gene_id`, `gene_symbol`, `evidence_type`, `source`, `notes` | One normalized gene per row; multiple supporting sources traceable to evidence table |
| `data/gene_sets/environment_evidence.tsv` | Above evidence identifiers plus `exposure`, `response_type`, `species`, `tissue`, `context` | Exposure and experimental context explicit for each assertion |
| `data/gene_sets/environment_response.tsv` | `gene_id`, `gene_symbol`, `exposure`, `tissue`, `context`, `evidence_type`, `source`, `notes` | Gene membership deduplicated without losing evidence/context links |
| `data/gene_sets/overlap.tsv` | `gene_id`, `gene_symbol`, `risk_source`, `environment_source` | Intersection computed by stable ID |
| `data/gene_sets/candidate_genes.tsv` | `gene_id`, `gene_symbol`, `is_risk`, `is_environment`, `is_overlap` | Union of both lists; flags consistent; distinguish risk-only, environment-only, overlap when disjoint summaries are needed |
| `data/metadata/samples.tsv` | `sample_id`, `subject_id`, `condition`, `tissue`, `accession`, `library_type`, `batch`, `junction_path`, `genome_build`, `annotation_version` | Unique sample IDs; case/control provenance; unavailable metadata explicitly marked |
| `data/processed/junctions.tsv` | `sample_id`, `gene_id`, `from`, `to`, `count` | Nonnegative finite support; stable endpoints; sample/gene joins valid |
| `data/processed/nodes.tsv` | `gene_id`, `node_id`, `chromosome`, `start`, `end`, `strand`, `representation` | Document coordinate convention, exon/splice-site choice, reference version |
| `results/sample_metrics.tsv` | `run_id`, `sample_id`, `gene_id`, `total_support`, `n_vertices`, `n_edges`, `cycle_rank`, `gradient_energy`, `cycle_energy`, `signal_energy`, `cycle_fraction`, `qc_status` | One row per eligible sample/gene; excluded or undefined values recorded explicitly |
| `results/junction_components.tsv` | `run_id`, `sample_id`, `gene_id`, `from`, `to`, `normalized_usage`, `gradient_component`, `cycle_component` | Same edge ordering as projection; signed components preserved |
| `results/gene_comparisons.tsv` | `gene_id`, `n_case`, `n_control`, `mean_case`, `mean_control`, `delta_cycle_fraction`, `p_value`, `q_value`, `status`, category flags | Counts refer to biological samples; missing inference never presented as significant |
| `results/ranked_genes.tsv` | Comparison fields plus `rank`, `ranking_rule` | Deterministic order; distinguish effect ranking from evidence of significance |
| `results/run_manifest.json` | Input checksums, schema/reference versions, command, tool version, parameters, seed, exclusions, solver settings | Another teammate can reproduce the run |

The proposal's illustrative `gene,from,to,count` input lacks sample identity. Preserve those semantics while adding `sample_id` for individual-sample analysis, or document a file-to-sample mapping. The input adapter must not silently pool patients. Missing measurements must be distinguished from confirmed zero junction support.

## Assignable workstreams and semester checklist

Each unchecked item is future work. Complete an item only when its output and review evidence exist.

### 1. Endometriosis genetics

**Goal:** produce the reproducible genetic-risk/susceptibility gene set described in proposal sections 3 and 23. **Owner:** Biology Lead; **reviewer:** Environment or Data Lead. **Input:** approved scope and later literature/database evidence. **Output:** risk evidence and risk membership TSVs above.

- [ ] G1. Write inclusion/exclusion criteria distinguishing inherited susceptibility evidence from expression-only association in `docs/gene_set_methods.md`.
- [ ] G2. Record which literature/database sources will be used; the proposal names no specific source, so do not claim one was preselected.
- [ ] G3. Create the evidence-table template with source identifiers, versions/access dates, and mapping notes.
- [ ] G4. Extract candidate evidence from the first source; retain locus/variant-to-gene mapping evidence where relevant.
- [ ] G5. Extract the remaining agreed sources using the same criteria.
- [ ] G6. Normalize identifiers to an agreed stable gene namespace/reference release; keep symbols separately.
- [ ] G7. Flag unmapped and one-to-many mappings; resolve or explicitly exclude them with reasons.
- [ ] G8. Merge duplicate gene membership while preserving every supporting evidence row.
- [ ] G9. Have a second teammate audit provenance, mapping, and inclusion consistency.
- [ ] G10. Export `endometriosis_risk.tsv`; check required fields, uniqueness, and join compatibility.

**Done:** every included gene has traceable evidence and a stable identifier; duplicates and unresolved mappings have a documented disposition. **Handoff:** Environment Lead merges sets; Data/Rust Leads consume stable IDs.

### 2. Environmental response

**Goal:** build the proposal's environmental-response prioritization set. **Owner:** Environment Lead; **reviewer:** Biology Lead. **Output:** environment evidence/membership, overlap, and candidate TSVs.

- [ ] E1. Record the selected exposure or exposure class and inclusion criteria; this practical scope is not fixed in the PDF.
- [ ] E2. Create an evidence template with exposure, dose/duration when available, species, tissue, assay, and response context.
- [ ] E3. Record the agreed evidence sources and extract the first source.
- [ ] E4. Extract genes from remaining sources, distinguishing measured response from indirect association.
- [ ] E5. Normalize stable gene IDs using the genetics workstream's namespace and release.
- [ ] E6. Record tissue/context mismatches and mapping ambiguities for review rather than hiding them.
- [ ] E7. Merge membership duplicates while preserving exposure-specific evidence records.
- [ ] E8. Audit provenance, required fields, and consistency; export `environment_response.tsv`.
- [ ] E9. Compute the stable-ID intersection and union; export `overlap.tsv` and `candidate_genes.tsv` with membership flags.
- [ ] E10. Check overlap membership against both source sets and report set sizes, including risk-only and environment-only groups.
- [ ] E11. Add the interpretation statement to methods: without direct patient exposure data, these genes prioritize analysis and do not establish exposure causation.

**Done:** environmental scope and provenance are explicit; both teams can reproduce candidate membership. A small overlap does not prevent analyzing the broader risk and environment sets already allowed by the proposal.

### 3. RNA-seq / data

**Goal:** deliver comparable endometriosis/control splice-junction measurements in the tool's format. **Owner:** Data Lead; **reviewers:** Biology and Rust Leads. No datasets are obtained during this roadmap task.

- [ ] D1. Record dataset suitability: tissue identity, case/control definition, number of biological subjects, library design, metadata completeness, and available junction/BAM/FASTQ files.
- [ ] D2. Select and document the exact public dataset/accessions; assess access, storage, and processing needs before acquisition.
- [ ] D3. Obtain the selected data later; record sources and checksums without committing large raw files to the repository.
- [ ] D4. Build `samples.tsv`; verify diagnoses, subject IDs, paired samples, tissue comparability, batch, and exclusions.
- [ ] D5. Record genome build, annotation release, coordinate convention, strandedness, and stable gene mapping.
- [ ] D6. Use existing junction counts if suitable; otherwise derive them using established processing software and record versions/commands.
- [ ] D7. Document whether the assay and extraction support back-splice detection. Unavailable detection is not evidence of absence.
- [ ] D8. QC junction support, ambiguous mappings, low coverage, duplicated records, missing values, and sample outliers; export `qc/junction_qc.tsv` and `qc/data_qc.md`.
- [ ] D9. Apply reviewed count/coverage criteria consistently; preserve raw-to-filtered exclusion counts.
- [ ] D10. Map endpoints and genes into standard tables; resolve duplicated junctions under an explicit aggregation rule.
- [ ] D11. Verify sample/gene joins and graph comparability; distinguish unobserved measurements from measured zeros.
- [ ] D12. Export a small integration subset plus the full standardized candidate-gene input; have Rust Lead validate both against the schema.

**Done:** sample labels and tissue eligibility are reviewed, preprocessing is reproducible, and junction TSVs pass the agreed input contract. Missing essential tissue/design metadata blocks confirmatory interpretation, not synthetic tool development.

### 4. Mathematics

**Goal:** implement the equations already specified in proposal sections 9–16 and validate section 25's examples. **Owner:** Math Lead; **reviewer:** Rust Lead. **Outputs:** `docs/math_spec.md`, future `tests/fixtures/`, and `tests/expected/`.

- [ ] M1. Specify deterministic vertex/edge indexing and biological orientation. For each source-to-destination edge, B1 has -1 at the source and +1 at the destination.
- [ ] M2. Specify normalization: the proposal's initial example is `w_e = count_e / total_gene_support` for each sample/gene. Confirm the initial choice; preserve raw total support for QC.
- [ ] M3. Define zero-total, empty-graph, disconnected-graph, duplicate-edge, and unsupported self-loop handling. Zero signal has an undefined energy fraction and must receive an explicit status.
- [ ] M4. Specify the least-squares problem minimizing `||F - B1^T phi||^2`; record how the solver handles nonunique potentials on disconnected components. The projected edge signal must remain well-defined.
- [ ] M5. Specify `F_grad = B1^T phi` and projection onto `ker(B1)` for `F_cycle`. An exact unweighted orthogonal least-squares residual lies in this kernel; numerically verify that membership before accepting it as the cycle component.
- [ ] M6. Specify `C = ||F_cycle||^2 / ||F||^2` for nonzero signals; preserve signed edge components even though the input counts are nonnegative.
- [ ] M7. Choose and record absolute/relative tolerances, solver convergence criteria, precision, and failure behavior using small known examples; never silently accept failed projection checks.
- [ ] M8. Verify reconstruction, `B1 F_cycle` near zero, gradient/cycle inner product near zero, and energy additivity with scale-aware tolerance. Check `0 <= C <= 1` up to roundoff.
- [ ] M9. Verify incidence dimensions/signs and graph cycle rank `m - n + k`, with `k` connected components, using hand-worked fixtures.
- [ ] M10. Prepare the linear chain and a branching tree: cycle rank zero and cycle component zero within tolerance.
- [ ] M11. Prepare the proposal's exon-skipping graph and compute its expected projection independently. It can have an undirected cycle-space component without a directed cycle; do not require zero merely because the directed graph is acyclic.
- [ ] M12. Prepare the directed three-edge back-splice cycle. An equal-weight circulation has zero gradient component and cycle fraction one.
- [ ] M13. Prepare fixed-seed noise fixtures and report how perturbations affect the measurement; distinguish weak extra edges from changes to weights on fixed support.
- [ ] M14. Run the proposal's planned back-splice weight sweep `0.01, 0.05, 0.10, 0.20, 0.40` later, with all other weights and normalization recorded; compare against independent expected projections rather than assuming a universal monotonic response.
- [ ] M15. Add zero-support, disconnected, and ordering-invariance cases; reversing an algebraic orientation must also reverse the corresponding signal coordinate.
- [ ] M16. Review reference expectations before porting them into Rust regression tests.

**Done:** expected answers and tolerances are documented and the eventual implementation passes all accepted identities. Cycle-space energy is a graph signal summary; it does not by itself establish a directed biological cycle or a circular RNA molecule. This makes the proposal's residual and interpretation cautions operational without changing its decomposition.

### 5. Rust tool — major software workstream

**Goal:** a reusable Rust CLI that takes standardized junction data through validated graph decomposition, comparison, ranking, and export. **Owner:** Rust Lead; **reviewers:** Math, Data, and Stats Leads. Follow the proposal's module responsibilities; keep input adapters separate from the graph/mathematics core so another compatible dataset can use the tool.

| Module / future file | Assignable implementation task | Acceptance evidence |
|---|---|---|
| `src/main.rs` | Define argument parsing, help, subcommands, exit codes, configuration | Help and malformed-command checks; no hard-coded local paths |
| `src/parser.rs` | Stream/read TSV; validate headers, types, identifiers, counts and metadata joins | Valid small fixture accepted; malformed row reports line/field and nonzero exit |
| `src/genes.rs` | Group by sample and gene; attach candidate categories | Shuffled input gives equivalent groups; samples never pooled implicitly |
| `src/graph.rs` | Build deterministic directed graphs and endpoint maps; apply reviewed duplicate/support rules | Expected vertices/edges and stable mapping to junction IDs |
| `src/matrix.rs` | Construct B1 and expose component/rank metadata | Source/destination signs and dimensions match reference fixtures |
| `src/hodge.rs` | Normalize, solve projection, verify identities, calculate energies/fraction | Synthetic expectations pass; failed/undefined cases have explicit status |
| `src/compare.rs` | Read eligible per-sample outputs, compute group differences, run specified comparisons | Known fixture effects; valid sample design and missingness behavior |
| `src/ranking.rs` | Rank by documented effects/evidence and select explanatory junctions | Deterministic ties; signed direction retained; no causal attribution |
| `src/output.rs` | Export TSVs, manifest, and concise human-readable summary | Schema-valid files and round-trip checks |

- [ ] R1. Record the settled `splice-girl` command name and agree on the supported input contract, output schema, and stable module interfaces.
- [ ] R2. Set up the Rust package later with documented compiler/toolchain and dependency versions.
- [ ] R3. Evaluate and select the numerical crate using singular/disconnected synthetic cases; prioritize correct projection for the MVP before scale optimizations.
- [ ] R4. Build CLI shell and error reporting.
- [ ] R5. Implement parser validation and sample/gene grouping.
- [ ] R6. Implement graph construction, raw support retention, normalization, and B1 generation.
- [ ] R7. Implement the graph-level projection and numerical acceptance checks.
- [ ] R8. Connect a single-sample analysis command to exported metrics and edge components.
- [ ] R9. Pass all synthetic acceptance fixtures end-to-end before real-data interpretation.
- [ ] R10. Implement one-gene case/control analysis and multi-gene comparison with subject metadata.
- [ ] R11. Implement the agreed statistics interface and multiple-testing outputs; record seeds and parameters.
- [ ] R12. Implement gene filtering/ranking and a documented junction contribution rule based on exported components. Decide aggregation and direction reporting before interpretation.
- [ ] R13. Implement exports and run manifests; avoid embedding machine-specific paths in examples.
- [ ] R14. Handle missing files, bad IDs, zero support, solver failure, insufficient samples, and existing output paths predictably; avoid silent overwrites.
- [ ] R15. Measure runtime/memory on the small and candidate-gene inputs; add sparse operations or gene-level parallelism as SHOULD HAVE improvements if scale requires them, preserving reproducibility.
- [ ] R16. Add user installation/input/tutorial documentation, license decision, versioned release, and a small redistributable example.
- [ ] R17. Have another teammate build from a clean checkout and reproduce the example without developer help; run formatting, lint, unit and integration checks appropriate to the implementation.

**Target commands:** the software is named Splice Girl and uses the command `splice-girl`. The subcommands and IDs below remain planned interfaces beyond the current basic synthetic prototype. The CLI may preserve the proposal's `cargo run -- compare ...` development interface.

```sh
splice-girl analyze sample.tsv
splice-girl gene GENE_ID --case endometriosis/ --control controls/
splice-girl compare --case endometriosis/ --control controls/
splice-girl rank --genes candidate_genes.tsv --case endometriosis/ --control controls/
```

`analyze` validates input and returns per-sample/per-gene metrics. `gene` exposes graph and junction components for one gene and compares eligible samples. `compare` produces sample-based group differences and statistical outputs. `rank` applies candidate membership and the recorded ranking rule. Directory inputs require explicit sample-to-metadata mappings. A successful command reports analyzed/excluded counts, output locations, and meaningful statuses; an input/numerical failure is actionable.

**Done:** a teammate can supply a new conforming junction table, obtain reproducible metrics/comparisons/ranks and trace them back to junctions using documented commands. A notebook or precomputed report alone does not satisfy this workstream.

### 6. Statistics and baseline comparison

**Goal:** quantify sample-level differences and assess whether the Hodge summary adds useful information beyond simpler measurements. **Owner:** Stats Lead; **reviewers:** Math and Biology Leads. **Outputs:** `docs/statistics_spec.md`, comparisons, `results/baseline_comparisons.tsv`, and `reports/statistics.md`.

- [ ] S1. Identify the biological unit and whether samples are independent, repeated, or paired; record any restrictions needed for valid label permutations.
- [ ] S2. Freeze eligible samples/genes, graph-support comparison policy, coverage rules, and missing-value handling before inspecting case/control rankings.
- [ ] S3. Receive one valid `C_g,s` per eligible sample/gene; verify uniqueness and exclusions. Do not treat junctions as independent patients.
- [ ] S4. Define the primary effect as case mean minus control mean cycle fraction, consistent with the proposal's difference; record group sizes and variability.
- [ ] S5. Specify permutation count, random seed, sidedness, exchangeability restrictions, and p-value calculation. If an implementation uses sampled permutations, use a documented finite-permutation correction rather than reporting a zero p-value.
- [ ] S6. Validate later with tiny known-label examples, no-difference fixtures, missing-data cases, and paired/restricted cases if applicable.
- [ ] S7. Run valid per-gene tests later; if biological replication or design is inadequate, export descriptive effects with inference marked unavailable.
- [ ] S8. Apply a recorded FDR procedure (proposed implementation: Benjamini–Hochberg) to the defined family of eligible tested genes; keep raw p-values, adjusted q-values, and untested statuses.
- [ ] S9. Specify baseline outputs from the same samples and filtering: junction-usage differences and graph cycle rank. Add back-splice support only when detection is supported. Expression comparisons are SHOULD HAVE if compatible expression measurements exist.
- [ ] S10. Compare baseline and cycle-energy effects/rankings; use the fixed-topology weight sweep to distinguish signal-weight information from cycle rank.
- [ ] S11. Check sensitivity to coverage/filtering and graph-support choices within documented implementation settings; retain all attempted settings, avoiding significance-driven selection.
- [ ] S12. Export final tables and explain whether the method offers additional information. A null result or no advantage over baselines remains a valid outcome.

**Done:** sample identity, effect direction, testing family, permutation design, baseline definitions, and exclusions are reproducible. Small p-values alone are not a biological explanation; interpretation follows QC and baseline review.

### 7. Visualization and presentation

**Goal:** explain the biology-to-tool pipeline to students with basic biology/computer science. **Owner:** Presentation Lead; **reviewers:** all workstream owners. **Outputs:** `figures/`, `docs/demo.md`, `presentation/`, and final report.

- [ ] V1. Draw the approved pipeline linking evidence, candidates, RNA junctions, Rust graphs, decomposition, comparison, and reusable output.
- [ ] V2. Draw a readable synthetic splice graph with endpoint names, junction support, and orientation matching the actual fixture.
- [ ] V3. Show original, gradient, and cycle components with consistent edge ordering and clear signed-component legends.
- [ ] V4. Prepare a one-gene case/control display showing individual samples, group summaries, support, and exclusions.
- [ ] V5. Prepare top-result figures only after statistics are validated; show effects and uncertainty/variation with p/q annotations where valid.
- [ ] V6. Link each selected gene's interpretation to risk/environment provenance; state observational and circular-RNA limitations.
- [ ] V7. Assign the biology introduction, data explanation, math explanation, CLI operator, statistics interpretation, and limitations segments.
- [ ] V8. Select one eligible real gene after analysis using a recorded rule that favors clear, supported output; a null-result gene is acceptable.
- [ ] V9. Package a tiny demo subset, exact commands, expected outputs, and a measured runtime target of seconds.
- [ ] V10. Rehearse on the presentation machine and prepare saved output, graph images, and a short recording as backup.
- [ ] V11. Verify that every figure matches a versioned run and that illustrative numbers are never labeled as experimental results.

**Done:** another teammate can run and explain the example, and the presentation makes the reusable tool and the limits of the evidence clear.

## Detailed execution phases

Timing is relative because no semester deadlines are recorded in the proposal. Map the opening, midpoint, and final gates to actual course dates at the first meeting. Dependencies govern execution; the phase numbering does not prohibit parallel work.

### Phase 1 — Repository/document organization

- **Goal:** make existing commitments and interfaces easy to assign and review.
- **Who:** Coordinator with all leads.
- **Exact tasks:** recover original Overleaf sources; record proposal/source locations; fill role slots; map course dates; agree first schemas and artifact ownership; create a task board referencing G/E/D/M/R/S/V IDs; record implementation decisions and review handoffs.
- **Inputs:** proposal PDF, this roadmap, later recovered Overleaf sources.
- **Outputs:** `docs/source_map.md`, `docs/data_contract.md`, `docs/decisions.md`, assigned task board.
- **Dependencies:** none for organizing; source recovery required to verify original-document compilation.
- **Definition of Done:** every workstream has an owner/reviewer, paths and schema versions are recorded, midpoint/final dates are assigned.
- **Possible blocker:** missing `.tex` source or team availability.
- **Fallback:** use the supplied PDF authority and standalone section until sources return; combine roles while keeping reviewers explicit.
- **Priority:** MUST HAVE. **Checkpoint:** interface and ownership review.

### Phase 2 — Gene-set work

- **Goal:** supply reproducible biological filters and interpretation provenance.
- **Who:** Biology and Environment Leads, Data reviewer.
- **Exact tasks:** execute G1–G10 and E1–E11; review source/mapping choices; export normalized membership, overlap, and candidate flags.
- **Inputs:** proposal scope; later selected evidence sources; shared ID namespace.
- **Outputs:** all `data/gene_sets/` artifacts; `docs/gene_set_methods.md`.
- **Dependencies:** phase 1 interfaces; two evidence collections can run simultaneously.
- **Definition of Done:** schema/provenance/duplicate audits pass and all category joins reproduce.
- **Possible blocker:** ambiguous mappings or tiny overlap.
- **Fallback:** document unresolved exclusions; retain separate risk and environmental sets as allowed in the proposal.
- **Priority:** MUST HAVE. **Checkpoint:** candidate-set handoff before real-data filtering.

### Phase 3 — Dataset preparation

- **Goal:** produce usable disease/control junction input.
- **Who:** Data Lead, Biology reviewer, Rust interface reviewer.
- **Exact tasks:** execute D1–D12; verify tissue/design suitability before downloading; organize metadata; obtain/derive counts; QC; export standard tables.
- **Inputs:** approved comparison, chosen public dataset, schema and reference conventions.
- **Outputs:** sample manifest, junction/node TSVs, preprocessing record, QC report.
- **Dependencies:** phase 1; phase 2 needed for final candidate subset, not initial dataset preparation.
- **Definition of Done:** acquisition and transformations are reproducible; sample labels, comparability, coverage, and input validity reviewed.
- **Possible blocker:** unavailable junction counts or unusable tissue metadata.
- **Fallback:** derive counts from available BAM/FASTQ using established software; if access or suitability fails, document a replacement public dataset meeting the same comparison criteria. Continue synthetic development meanwhile.
- **Priority:** MUST HAVE. **Checkpoint:** data-readiness review before phase 6.

### Phase 4 — Synthetic graph and mathematical validation

- **Goal:** establish known answers before interpreting real biology.
- **Who:** Math Lead with Rust reviewer.
- **Exact tasks:** execute M1–M16; prepare linear/tree, exon-skipping, back-splice, noise, and changing-weight cases; specify numerical acceptance and independent expected values.
- **Inputs:** proposal equations and graph contract.
- **Outputs:** math specification, fixture inputs, expected components/metrics, validation record.
- **Dependencies:** phase 1 contract; independent of gene selection and real-data acquisition.
- **Definition of Done:** every required fixture has reviewed expectations; identities and failure criteria are explicit and verified in the later implementation.
- **Possible blocker:** singular solver behavior or confusion between directed cycles and graph cycle space.
- **Fallback:** work through tiny incidence matrices and rank-aware least-squares examples; retain the exact graph-level model and defer scale optimization.
- **Priority:** MUST HAVE. **Checkpoint:** mathematical acceptance gate before real-data claims.

### Phase 5 — Rust MVP implementation

- **Goal:** run standardized synthetic data through the core tool.
- **Who:** Rust Lead with Math and Data reviewers.
- **Exact tasks:** implement R1–R9 and the initial module interfaces; validate parser errors, grouping, graph/B1 construction, projection, metrics, and exports.
- **Inputs:** stable schema, math specification, expected synthetic fixtures.
- **Outputs:** Rust package, working analysis command, tests, synthetic metrics/components.
- **Dependencies:** phases 1 and 4; parser work can start while fixture expectations are reviewed.
- **Definition of Done:** a documented command reproduces a reviewed synthetic decomposition with no hard-coded gene/sample assumptions and passes its checks.
- **Possible blocker:** numerical crate or sparse integration complexity.
- **Fallback:** use a validated small-graph solver appropriate to the same equations; defer parallelism and large-scale optimization.
- **Priority:** MUST HAVE; performance improvements SHOULD HAVE. **Checkpoint:** midpoint demonstration.

### Phase 6 — Real splice-data integration

- **Goal:** connect real candidate-gene data to the validated engine.
- **Who:** Data and Rust Leads, Math reviewer.
- **Exact tasks:** select an eligible small integration subset later; validate IDs and metadata; run sample-level decomposition; inspect edge mappings and QC statuses; expand to the candidate set only after the small run passes.
- **Inputs:** phase 2 gene files, phase 3 data, phase 5 tool.
- **Outputs:** integration log, accepted inputs, `sample_metrics.tsv`, junction components, run manifest.
- **Dependencies:** phases 2–5 and synthetic gate.
- **Definition of Done:** real rows trace back to source junctions; numerical checks pass; failures/exclusions are counted and reproducible.
- **Possible blocker:** annotation mismatch, low coverage, or runtime.
- **Fallback:** harmonize identifiers/reference versions; use a documented smaller eligible candidate subset under the proposal's MVP. Do not substitute synthetic results for real evidence.
- **Priority:** MUST HAVE. **Checkpoint:** reviewed one-gene real-data run.

### Phase 7 — Endometriosis/control analysis

- **Goal:** generate sample-preserving comparisons and descriptive ranks.
- **Who:** Rust and Stats Leads, Data reviewer.
- **Exact tasks:** verify eligible biological samples; calculate each gene's per-sample cycle fraction; summarize case/control groups; compute signed differences; produce component/junction summaries and deterministic descriptive ranks.
- **Inputs:** validated real outputs, metadata, frozen comparison/support rules.
- **Outputs:** sample and comparison tables, provisional effect rankings with inference status explicit.
- **Dependencies:** phase 6; S1–S4 contract.
- **Definition of Done:** every comparison reports sample counts, direction, missingness, and support; no implicit patient pooling or significance claims from rank alone.
- **Possible blocker:** insufficient valid samples or inconsistent graph support.
- **Fallback:** report descriptive comparisons and exclusions; resolve support policy consistently before inference; record unmet inferential requirements.
- **Priority:** MUST HAVE. **Checkpoint:** comparison-table review.

### Phase 8 — Statistics and baseline comparison

- **Goal:** evaluate differences and the added information of the method.
- **Who:** Stats Lead with Math/Data reviewers.
- **Exact tasks:** execute S5–S12; validate permutation implementation; apply design-appropriate permutations and FDR; compute simpler baselines on the same eligible inputs; document sensitivity and null outcomes.
- **Inputs:** phase 7 sample outputs, sample design, frozen statistical specification.
- **Outputs:** p/q/effect tables, baseline comparisons, statistical methods/results record.
- **Dependencies:** phase 7 and adequate replication/design for inference.
- **Definition of Done:** seeds, testing family, sample restrictions and baseline definitions reproduce; infeasible inference is explicitly unavailable.
- **Possible blocker:** nonexchangeable labels, low replication, or no benefit over baselines.
- **Fallback:** restrict permutations when justified by design; otherwise give descriptive results and an unmet-inference note. Report no added value honestly.
- **Priority:** MUST HAVE; extra compatible expression baseline SHOULD HAVE. **Checkpoint:** statistical release gate before interpretation.

### Phase 9 — Result interpretation

- **Goal:** connect validated patterns to candidate-gene evidence.
- **Who:** Biology and Environment Leads with Stats reviewer.
- **Exact tasks:** join validated results to provenance; inspect supported junction contributions; compare risk/environment/overlap categories descriptively; write evidence-linked summaries; document alignment, coverage, causality, and method limitations.
- **Inputs:** phases 2 and 8, graph/component views.
- **Outputs:** `reports/interpretation.md`, figure captions, gene evidence summaries.
- **Dependencies:** final statistical and QC review.
- **Definition of Done:** all claims link to actual outputs and source evidence; environmental prioritization and circular-RNA cautions are explicit; unsupported hypotheses remain hypotheses.
- **Possible blocker:** weak, null, or biologically ambiguous findings.
- **Fallback:** explain the supported null/ambiguous result and baseline comparison; use a clear eligible gene to demonstrate functionality without claiming a discovery.
- **Priority:** MUST HAVE. **Checkpoint:** biology/statistics joint review.

### Phase 10 — Tool cleanup and usability

- **Goal:** deliver reusable software beyond this study.
- **Who:** Rust Lead; independent teammate as user reviewer.
- **Exact tasks:** finish R10–R17; document installation, schemas, commands, limitations, parameters and examples; add robust errors/manifests; test clean-checkout execution; measure runtime and produce a versioned release.
- **Inputs:** validated engine, integrated comparisons, real/synthetic example data.
- **Outputs:** documented Rust CLI, test suite, example package, release notes and version.
- **Dependencies:** phase 5 onward; final acceptance after phases 6–8.
- **Definition of Done:** a teammate runs a new conforming input through analysis/comparison/ranking and reproduces shipped examples without manual code edits.
- **Possible blocker:** platform/dependency problems or scale.
- **Fallback:** document a tested environment and provide reproducible source/build instructions; demonstrate a smaller real subset. Performance optimization cannot replace correctness.
- **Priority:** MUST HAVE; sparse/parallel scale improvements SHOULD HAVE. **Checkpoint:** release candidate acceptance.

### Phase 11 — Live demo

- **Goal:** show a real reusable command completing in seconds.
- **Who:** Presentation Lead and CLI operator; Math/Biology speakers.
- **Exact tasks:** execute V7–V10; package small inputs; record commands and expected outputs; rehearse timing; verify backup on the presentation machine.
- **Inputs:** release candidate, validated synthetic and real-gene examples, reviewed interpretation.
- **Outputs:** `docs/demo.md`, `examples/demo/`, saved outputs, backup images/recording.
- **Dependencies:** phases 9–10; storyboarding starts earlier.
- **Definition of Done:** live commands finish in seconds during rehearsal and a backup can explain the same output without network access.
- **Possible blocker:** machine failure or unexpected latency.
- **Fallback:** use saved outputs and recording labeled as precomputed; do not run the full dataset live.
- **Priority:** MUST HAVE. **Checkpoint:** timed rehearsal.

### Phase 12 — Final report and presentation

- **Goal:** deliver the study, reusable tool, and clear evidence of validation.
- **Who:** all leads; Presentation Lead assembles, Coordinator checks completeness.
- **Exact tasks:** finish V1–V11; assemble methods/results/limitations, source citations, pipeline and graph figures; link software version and commands; review all captions against outputs; compile any edited LaTeX and resolve formatting warnings.
- **Inputs:** phases 8–11 outputs and preserved proposal.
- **Outputs:** final report, slides, figures, tool release, reproducibility record, tested Overleaf documentation when sources are available.
- **Dependencies:** validated results and demo readiness.
- **Definition of Done:** team can distinguish demonstrated functionality, real findings, illustrative examples, and unmet goals; submission opens and figures/tables are readable.
- **Possible blocker:** late results or unavailable original LaTeX.
- **Fallback:** present the validated scope and explicit limitations; use the separately tested roadmap section while retaining the original-source integration task.
- **Priority:** MUST HAVE. **Checkpoint:** final acceptance review.

### Phase 13 — Research extensions

- **Goal:** execute only the optional extensions already listed in the proposal.
- **Who:** relevant workstream owners after core release.
- **Exact tasks:** if time and evidence permit, document patient-level exposure comparison; define biologically meaningful faces before introducing B2; investigate genome-wide scaling or a compatible other-disease application. Record prerequisites and separate outputs before work begins.
- **Inputs:** completed core milestone; suitable exposure data or justified higher-order structure, as applicable.
- **Outputs:** `extensions/` documentation and separately validated extension artifacts.
- **Dependencies:** core tool, statistical checks, and final presentation readiness.
- **Definition of Done:** any attempted extension has its own validated inputs, method conditions and limitations; no unsupported curl/harmonic interpretation is attached to the graph-only MVP.
- **Possible blocker:** absent exposure metadata, unjustified faces, or insufficient time.
- **Fallback:** defer the extension and retain the approved graph-level tool.
- **Priority:** STRETCH. **Checkpoint:** optional go/no-go after core acceptance.

## Milestones and immediate assignment queue

### Start together at the next meeting

1. Coordinator: recover the Overleaf source, assign names to roles, and place midpoint/final dates on the task board (phase 1).
2. Biology Lead: draft G1–G3, the genetic-evidence criteria and table contract.
3. Environment Lead: draft E1–E2, the exposure scope and context/provenance contract.
4. Data Lead: start D1–D2, the dataset suitability/access review; document exact accessions before acquisition.
5. Math Lead: start M1–M3 and hand-work the linear and three-edge-cycle expectations.
6. Rust Lead: draft R1 interfaces with Data/Math; implementation begins later under phase 5.
7. Stats Lead: start S1–S2 and the sample-output contract; Presentation Lead starts V1 and the storyboard.

Items 2–7 can proceed simultaneously, with a short joint review to settle shared identifiers, graph representation, and schemas. These are assignments for the group after this documentation task, not actions performed now.

### Midpoint gate

- [ ] Gene-set criteria and environmental scope recorded; candidate files with provenance exist.
- [ ] Dataset selected and accessible; metadata and a small standardized junction subset prepared; any remaining full-data preparation has an owner and date.
- [ ] Junction, sample, node, and output schemas stable and versioned.
- [ ] Synthetic fixtures and independent expected projections prepared.
- [ ] Mathematical identities verified on required synthetic cases with documented tolerances.
- [ ] Rust parser, sample/gene grouping, graph construction, and B1 creation work.
- [ ] One documented synthetic command runs normalization, decomposition, metrics, and export end-to-end.
- [ ] Team demonstrates the run and reviews blockers for real-data integration.

The midpoint artifact is a tested synthetic Rust MVP plus usable biological/data inputs and a stable contract. Final biological results are not required at midpoint. If data access is blocked, record that midpoint data requirement as unmet; synthetic progress does not make it complete.

### Final gate

- [ ] Real endometriosis/control data integrated with reproducible metadata and QC.
- [ ] Per-sample graph decomposition, group comparison, candidate filtering and ranking work.
- [ ] Effect sizes, valid permutation tests, FDR, and baseline comparisons are complete where design supports inference; any unmet requirement is named explicitly.
- [ ] Supported findings, including null outcomes, are interpreted with provenance and limitations.
- [ ] Reusable Rust CLI accepts new conforming data and exports metrics, junction components, comparisons, ranks, and a run manifest.
- [ ] Installation, input specification, example, numerical checks, error handling, and release version are documented and independently exercised.
- [ ] One synthetic example and one real gene can run live in seconds; backup is ready.
- [ ] Pipeline, graph, sample comparison, and final result figures match validated outputs.
- [ ] Report, slides, references, and available LaTeX documentation are checked.

If scope must shrink, apply proposal section 31: fewer eligible genes and a small real case/control comparison using the same normalization, B1, projections, ranking, and Rust demonstration. A synthetic-only tool is useful progress but does not complete the planned real-data deliverable. Missing replication cannot be repaired by pooling junctions or manufacturing inference.

## Concrete live demo

Prepare the full analysis beforehand; live commands use a tiny local package. The following sequence is a rehearsal target, not a measured runtime claim.

| Step | Speaker/action | What the audience sees | Acceptance |
|---|---|---|---|
| 1 | Math speaker introduces tiny chain and cycle | Named nodes, directed junctions, known counts | No unexplained notation |
| 2 | Operator runs `splice-girl analyze examples/demo/synthetic.tsv` | Gradient/cycle metrics and exported components | Command finishes in seconds; matches expected fixture |
| 3 | Math speaker explains decomposition | Which signal is explained by a vertex potential and which lies in cycle space | Clarify this is a mathematical projection |
| 4 | Operator runs `splice-girl gene GENE_ID --case examples/demo/case/ --control examples/demo/control/` | One real gene, nodes/junctions, sample fractions, group difference | Stable local data; command finishes in seconds |
| 5 | Biology/Stats speakers interpret | Supported junction contributions, sample variation, relevant evidence and inference status | No circular-RNA or exposure-causation claim from graph energy |
| 6 | Operator shows precomputed larger ranking and help | Reusable input/output interface and full-study context | Precomputed results labeled with run version |

Freeze the actual gene ID and paths only after integration and review. Use the same tool version and small files during rehearsal and presentation. Record observed timings in `docs/demo.md`. Backup: saved terminal output, matching figures and a short recording from the validated run, all available offline. If no eligible real gene is ready, demonstrate the synthetic tool and explicitly state that the real-gene milestone is unmet.

## Open Implementation Decisions

These are execution details left open by the proposal or necessary to make its interfaces reproducible. The research question, endometriosis focus, Rust deliverable, graph-level decomposition and optional status of B2 are settled.

| Decision | Owner | Resolve by | Record / boundary |
|---|---|---|---|
| Recover main `.tex`, preamble and exact source layout | Coordinator | Original-document integration | `docs/source_map.md`; do not reconstruct/rewrite proposal silently |
| Course dates and role assignments | Coordinator/team | Phase 1 review | Task board |
| Exact genetic evidence sources and mapping criteria | Biology | Before G4 | `gene_set_methods.md`; preserve susceptibility focus |
| Exact exposure/class and evidence context | Environment | Before E3 | Same methods file; prioritization without patient exposure data |
| Exact public dataset/accessions and eligible tissue/subjects | Data/Biology | Before acquisition | Data methods and manifest |
| Stable ID release, exon vs splice-site nodes, coordinates and duplicate rules | Data/Math/Rust | Before graph integration | Data contract; proposal permits either node representation |
| Minimum support/coverage and missing-vs-zero policy | Data/Stats | Before real comparisons | QC methods and run parameters |
| Initial normalization | Math/Data | Before fixture acceptance | Proposal suggests within-gene total normalization; confirm and document zero-total handling |
| Comparable graph support across samples | Math/Data/Stats | Before phase 6 | Decide shared support or explicitly justified varying support; avoid label-dependent filtering |
| Numerical crate, solver, precision and tolerances | Rust/Math | Before M/R acceptance | Math specification; equations unchanged |
| Exact subcommand flags, output and manifest formats | Rust/Data/Stats | Before public interface stabilizes | The `splice-girl` command name is settled; command details above remain conceptual |
| Permutation design/count/seed/sidedness, FDR family | Stats | Before real testing | Statistics specification; design must support exchangeability |
| Baseline definitions, ranking/tie rules and junction contribution aggregation | Stats/Math/Rust | Before phase 7–8 outputs | Distinguish effects, significance and explanatory components |
| Example redistribution, software license, release platform | Coordinator/Rust | Before packaging | Document source terms and tested build environment |

## Practical risks and fallbacks

| Problem / trigger | Owner | Fallback and required record |
|---|---|---|
| Original LaTeX missing | Coordinator | Use PDF authority; compile include-only section in a basic wrapper; integration remains unverified until source returns |
| Overlap set tiny | Biology/Environment | Analyze risk, environment and overlap membership as already allowed; report actual group sizes |
| Gene IDs/reference builds disagree | Data | Reconcile versions and mapping; retain ambiguous exclusions and evidence trail |
| Junction data unavailable | Data | Derive from available BAM/FASTQ using established tools; record cost/versions; use another suitable public dataset if necessary |
| Tissue or batch prevents defensible case/control inference | Data/Stats | Review eligible samples/design; report limitation and defer unsupported inference; do not relabel unsuitable controls |
| Back-splice detection unsupported | Data/Math | State detection limits; use available junction topology without interpreting absent back-splices as biological absence |
| Low support/noisy edges | Data/Math | Apply documented thresholds and sensitivity checks; record exclusions and denominator failures |
| Solver residual fails kernel/orthogonality checks | Math/Rust | Fail the calculation explicitly, debug small fixtures and solver rank handling before real interpretation |
| Cycle fraction adds little beyond cycle rank/junction usage | Stats | Report comparable baseline performance honestly; retain functioning tool and validated result |
| Too few biological samples | Stats | Descriptive effects only; p/q unavailable; inferential goal remains unmet |
| Full dataset too large | Rust/Data | Use candidate filtering and small validated subset; profile before sparse/parallel improvements |
| Full Hodge extension lacks meaningful faces | Math | Defer B2; deliver graph-level MVP |
| Live command fails or is slow | Presentation/Rust | Use precomputed output/recording labeled clearly; demonstrate one gene, never full processing |
| No compelling significant gene | Biology/Presentation | Demonstrate a clear QC-passing real gene and explain null findings without implying discovery |

## Documentation validation

`roadmap.tex` contains a short section with standard sectioning, a simple margin-aware table, and verbatim command examples. It declares no document class, document environment, custom fonts, or package requirements. Include it using `\input{roadmap}` from the original main document once that source is available. If the team later moves it into a `sections/` folder, update the input path to `\input{sections/roadmap}`.

The accompanying `roadmap-preview.pdf` is a standalone rendering in a temporary standard article wrapper. It is a preview of the new section, not a rebuilt proposal. Original Overleaf integration cannot be certified without the missing `.tex` source.

Validation completed with Tectonic 0.16.9 using cached resources and a temporary 11-point `article` wrapper with standard one-inch margins. The three-page preview compiled successfully; the log contained no overfull/underfull boxes or LaTeX warnings, and all rendered pages were visually checked. The include-only section contains no document wrapper commands. The initial sandboxed compiler startup failed; the same cached-resource build succeeded outside the sandbox. No proposal source or existing project file was modified.
