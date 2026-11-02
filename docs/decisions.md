# Splice Girl decisions log

This is the dated record of project-wide choices. Add a new dated entry when a
decision changes; do not rewrite an old entry. Open questions stay open until
the team records an owner and a reasoned choice.

## Confirmed decisions

| Date | Decision | Reason and evidence |
|---|---|---|
| 2026-09-02 | The scientific question is whether splice-graph organization differs between endometriosis and control endometrial samples. | Proposal in [`assets/hodgetheorynew.pdf`](../assets/hodgetheorynew.pdf). |
| 2026-10-20 | Use four control and four endometriosis samples from GSE315857 for the first study cohort. | The exact accessions and labels are frozen in [`data/metadata/samples.tsv`](../data/metadata/samples.tsv). |
| 2026-10-22 | Use the GENCODE v50 GRCh38.p14 primary assembly and its matching comprehensive GTF. | The reference pair and checksum record are documented in [`docs/data_methods.md`](data_methods.md). |
| 2026-10-22 | Represent an exon as a node and an observed splice junction as a directed edge; orient edges in transcript direction and use 1-based closed coordinates. | Frozen graph conventions in [`docs/data_methods.md`](data_methods.md). |
| 2026-10-22 | Keep raw reads, genome references, indexes, BAM files, and preprocessing runs outside Git. | They are large and reproducibility is preserved through paths and checksums in run records; see [`.gitignore`](../.gitignore). |
| 2026-10-28 | Keep descriptive case-minus-control effects separate from inferential statistics and gene ranking. | The current Rust output is sample-preserving and does not claim significance; see [`docs/statistics_spec.md`](statistics_spec.md). |
| 2026-10-29 | Treat each subject/sample as the biological unit in the fixed eight-sample cohort. | The statistics specification records an unpaired 4-versus-4 design and its exact permutation space. |
| 2026-11-02 | The minimum supported Rust version is 1.85.0, and the package has no external Rust dependencies. | `Cargo.toml`, the isolated installation check, and [`docs/software_validation.md`](software_validation.md). |

## Open decisions that still need a team answer

| Question | Why it matters | What must be recorded before implementation |
|---|---|---|
| Which exact support and coverage thresholds define an eligible junction? | Thresholds change which graphs and samples enter the analysis. | Threshold values, rationale, and a small sensitivity plan. |
| Which splice-aware aligner is final for production junction extraction? | The helper script currently plans STAR, while the frozen data methods name HISAT2; Bowtie2 is currently documented only for optional PhiX cleanup. | One chosen aligner, version, index command, splice-site options, and a checked pilot before D6 is marked done. |
| Will the final permutation test be one-sided or two-sided? | The direction of a cycle-fraction change is not known in advance. | The chosen alternative and the person approving it. |
| Which reviewed candidate-gene table is final? | Ranking and interpretation must use a frozen, source-linked set. | Stable gene IDs, evidence sources, and inclusion rules. |
| What are the course midpoint and final deadlines? | The task board cannot assign realistic due dates without them. | Dates, owner, and the required artifact at each checkpoint. |
| Where is the parent Overleaf proposal source? | The local roadmap section cannot be compiled alone. | Main `.tex` filename, folder, and compile command. |

An open decision is not a failure. It is safer to leave an item visible than to
hide an assumption inside code or a final result.
