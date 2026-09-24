# Topological Analysis of RNA Splicing in Endometriosis

**Integrating Genetic Susceptibility, Environmental Response, and Hodge Theory**  
BINF 2111 research project, 2026

This project asks whether genes associated with endometriosis show differences in RNA-splicing topology between endometriosis and normal endometrial tissue. Genetic susceptibility and environmental-response evidence will prioritize genes. The study will represent splice junctions as weighted graphs and use graph-level Hodge decomposition to examine how their junction-usage patterns differ between samples.

The final computational deliverable is a reusable command-line tool named **Splice Girl**, written primarily in Rust. Its command name is `splice-girl`, and the endometriosis study will be its first application.

## Current status

The scientific proposal and execution roadmap are available. A basic Rust prototype now performs the graph-level decomposition on two synthetic fixtures. Gene-set construction, dataset preparation, real-data integration, statistics, ranking, and research analysis have not started. Most commands and output paths in the roadmap remain future deliverables.

The approved proposal remains the scientific source of truth. The roadmap translates it into tasks, dependencies, and completion checks.

## Start here

| Document | Purpose |
|---|---|
| [Project proposal](assets/hodgetheorynew.pdf) | Full research question, biological motivation, mathematical method, validation strategy, and limitations |
| [Execution roadmap](ROADMAP.md) | Detailed team responsibilities, 13 execution phases, task checklists, file contracts, milestones, and fallbacks |
| [Beginner task guide](delegations_and_tasks/tasks.md) | Very detailed step-by-step guide for every workstream task |
| [LaTeX roadmap section](roadmap.tex) | Short roadmap section for inclusion in the existing Overleaf proposal |
| [Roadmap PDF preview](roadmap-preview.pdf) | Compiled standalone preview of the LaTeX section |
| [Basic math specification](docs/math_spec.md) | Orientation, normalization, decomposition equations, validation identities, and first expected examples |
| [Input data contract](docs/data_contract.md) | Required TSV columns, validation rules, examples, and current prototype limits |
| [Contribution logs](logs/README.md) | Dated record of Summer's and Sydney's repository contributions |
| [Splice Girl startup art](assets/splice_girl_banner.txt) | Detailed ASCII figure printed when the CLI starts |

## Run the basic synthetic prototype

The current command accepts a tab-separated file with exactly these columns:

```text
sample_id	gene_id	from	to	count
```

Run the chain and cycle examples:

```sh
cargo run -- tests/fixtures/linear.tsv
cargo run -- tests/fixtures/cycle.tsv
```

To start Splice Girl with its interactive menu, run:

```sh
cargo run -- --start
```

The menu currently contains one action, **Analyze a junction TSV**, plus quit. Choose `1`, enter a path such as `tests/fixtures/cycle.tsv`, and press Enter to return to the menu. The direct file-path command remains available for scripts and quick checks.

Run the automated checks:

```sh
cargo test
```

The chain should report a cycle fraction near zero. The equal-weight directed cycle should report a cycle fraction near one. Both runs also check reconstruction, cycle-space membership, and orthogonality. The prototype currently rejects zero-total groups, duplicate directed edges, and self-loops with explicit errors.

These `cargo run` commands build and run the `splice-girl` executable. After a release build, the equivalent direct command begins with `./target/release/splice-girl`.

Splice Girl prints its detailed ASCII startup figure before each run. The figure is stored in `assets/splice_girl_banner.txt` and is shown in a soft terminal color when output is connected to a terminal.

## How the project connects

```text
Genetic-risk evidence + Environmental-response evidence
                          |
                          v
                    Candidate genes
                          |
Endometriosis/control RNA-seq --> Standardized junction tables
                                          |
                                          v
                                 Rust splice-graph engine
                                          |
                                          v
                              Gradient/cycle-space decomposition
                                          |
                                          v
                              Per-sample measurements
                                          |
                                          v
                         Case/control statistics and baselines
                                          |
                                          v
                         Gene ranking, interpretation, and demo
```

In the graph, vertices represent exons or splice sites, edges represent observed splice junctions, and normalized read support describes junction usage. The decomposition separates the edge signal into a gradient component explained by values assigned to vertices and a component in the graph's cycle space. Synthetic examples will validate the calculations before real-data interpretation.

Environmental-response evidence is initially a gene-prioritization component. Without patient-level exposure measurements, it cannot establish that an exposure caused a splicing change. Likewise, a cycle-space component does not establish the existence of circular RNA. Full simplicial Hodge decomposition remains an optional extension.

## Team workstreams

| Workstream | Main responsibility |
|---|---|
| Endometriosis genetics | Build a stable-ID susceptibility gene set with traceable evidence |
| Environmental response | Record exposure and tissue context; build response-gene and overlap sets |
| RNA-seq / data | Verify sample comparability, prepare junction counts, and deliver quality-checked standard inputs |
| Mathematics | Specify graph orientation and projections; validate synthetic examples and numerical identities |
| Rust tool | Implement input validation, graph construction, decomposition, comparison, ranking, and export |
| Statistics | Compare biological samples, quantify effects, run valid permutation tests, correct multiple testing, and assess baselines |
| Visualization / presentation | Prepare understandable figures, reviewed interpretations, and a small live demonstration |

Names and assignments belong in the roadmap's **Who is doing what** table. Team members may share roles; each deliverable should have an owner and reviewer.

## Next steps for the group

1. Assign workstream owners and reviewers, and record course milestone dates.
2. Recover the original Overleaf source files; only the proposal PDF is currently present.
3. Document genetic-evidence criteria and the environmental exposure scope.
4. Review public dataset suitability, tissue comparability, metadata, and access requirements.
5. Agree on identifiers, junction-table fields, graph representation, and numerical validation expectations.
6. Assign the roadmap's first tasks and update their statuses as work begins.

Genetics and environmental evidence preparation can proceed in parallel. Dataset review can run alongside synthetic mathematics and Rust interface design. Real statistical analysis depends on validated per-sample tool outputs; final biological interpretation follows statistical and baseline review.

## Milestones

**Midpoint:** candidate-gene files with provenance, a selected accessible dataset and prepared small input subset, stable schemas, reviewed synthetic examples, and a Rust MVP that runs one synthetic decomposition from input through export. Final biological results are not required at midpoint.

**Final:** a documented reusable Rust CLI that processes real junction data, calculates per-sample components, compares case/control groups, ranks genes, identifies relevant junctions, and exports reproducible outputs. The study should include validated statistics where the sample design supports inference, baseline comparisons, interpreted findings, final figures, and a real-gene demonstration that runs in seconds.

The live demonstration will use tiny local inputs. Larger analyses will be precomputed, with saved output and a recording available as backup. Null findings or no advantage over simpler baselines remain valid outcomes and should be reported honestly.

## Overleaf documentation

Add `roadmap.tex` alongside the original main `.tex` file and include it where the execution roadmap belongs:

```latex
\input{roadmap}
```

The section inherits the parent document's formatting and contains no document class or document environment. It uses standard LaTeX sectioning, a simple table, and verbatim command examples.

The supplied three-page preview compiled successfully with Tectonic in a temporary standard article wrapper and was checked for formatting issues. Integration with the original Overleaf document remains unverified until its source files are restored. The preview is not a rebuilt version of the proposal.

See [ROADMAP.md](ROADMAP.md) for open implementation decisions, detailed acceptance criteria, and practical fallbacks.
