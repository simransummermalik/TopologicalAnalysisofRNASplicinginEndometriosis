# Topological Analysis of RNA Splicing in Endometriosis

**Integrating Genetic Susceptibility, Environmental Response, and Hodge Theory**  
BINF 2111 research project

This project asks whether genes associated with endometriosis show differences in RNA-splicing topology between endometriosis and normal endometrial tissue. Genetic susceptibility and environmental-response evidence will prioritize genes. The study will represent splice junctions as weighted graphs and use graph-level Hodge decomposition to examine how their junction-usage patterns differ between samples.

The final computational deliverable is a reusable command-line tool written primarily in Rust. The endometriosis study will be its first application.

## Current status

The scientific proposal and execution roadmap are available. Software implementation, gene-set construction, dataset preparation, and research analysis have not started in this repository. There is no runnable CLI or installation procedure yet. Commands and output paths in the roadmap describe future deliverables.

The approved proposal remains the scientific source of truth. The roadmap translates it into assignments, dependencies, and completion checks.

## Start here

| Document | Purpose |
|---|---|
| [Project proposal](assets/hodgetheorynew.pdf) | Full research question, biological motivation, mathematical method, validation strategy, and limitations |
| [Execution roadmap](ROADMAP.md) | Detailed team responsibilities, 13 execution phases, task checklists, file contracts, milestones, and fallbacks |
| [LaTeX roadmap section](roadmap.tex) | Short roadmap section for inclusion in the existing Overleaf proposal |
| [Roadmap PDF preview](roadmap-preview.pdf) | Compiled standalone preview of the LaTeX section |

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
