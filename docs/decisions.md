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
| 2026-11-12 | Candidate-gene evidence uses the versioned GENCODE v50 `gene_id` as its stable ID, matching the junction tables. The NCBI GRCh38.p14 records in `human-reference-genome/` are not the canonical reference. | Reaffirms the 2026-10-22 reference decision so evidence tables join to Splice Girl outputs; see [`docs/gene_set_methods.md`](gene_set_methods.md). |
| 2026-11-12 | Genetic-risk evidence is limited to endometriosis associations at genome-wide significance (p < 5×10⁻⁸); the gene is the one reported or mapped by the source, and loci without a single reported gene are marked ambiguous rather than guessed. Expression-only evidence never counts as risk. | G1 criteria in [`docs/gene_set_methods.md`](gene_set_methods.md). |
| 2026-11-12 | Approved risk-evidence sources: the GWAS Catalog (endometriosis trait) and the largest endometriosis GWAS meta-analysis. | G2 source list in [`docs/gene_set_methods.md`](gene_set_methods.md); exact release and citation are recorded at extraction. |
| 2026-11-12 | The environmental-response scope is dioxin-like aryl-hydrocarbon-receptor (AhR) ligands: TCDD and other dioxin-like compounds. | E1 scope in [`docs/gene_set_methods.md`](gene_set_methods.md). Narrow class with established endometriosis literature. |
| 2026-11-12 | The "largest endometriosis GWAS meta-analysis" risk source is Koller, Cormand et al., *Nature Genetics* (2026), doi:10.1038/s41588-026-02582-2 (~1.4 million women, 105,869 cases, 80 risk loci). It supersedes Rahmioglu et al. 2023 (doi:10.1038/s41588-023-01323-z; 42 loci). Before extraction, confirm from the full text that the threshold is p < 5×10⁻⁸ and that a lead-variant table with reported or mapped genes exists. | Larger and newer; matches the "largest" rule. Comparison recorded in [`docs/gene_set_methods.md`](gene_set_methods.md). |
| 2026-11-12 | From Koller et al. (*Nat Genet* 58(5):1051–1061), extract only the cross-ancestry endometriosis combined-definition meta-analysis (105,869 cases; 80 loci at P < 5×10⁻⁸; Supplementary Table 4). Ancestry-specific, clinical-only, self-reported-only, and adenomyosis analyses are excluded. | The paper's headline analysis; its threshold matches G1. Full-text checks recorded in [`docs/gene_set_methods.md`](gene_set_methods.md). |
| 2026-11-12 | The dioxin-like definition stays anchored to the WHO 2005 TEF table (29 congeners). | Verified against Van den Berg et al. 2006, Table 1. The 2022 WHO update's congener table could not be accessed; its abstract describes a reevaluation of the 2005 TEFs. The project uses the congener list, not TEF values. |
| 2026-11-12 | Approved environmental-response source: CTD (Comparative Toxicogenomics Database) curated chemical–gene interactions, limited to TCDD plus a fixed list of dioxin-like compounds written before extraction, *Homo sapiens* only, and measured-response interaction types only. Rodent data may enter later only as a labeled indirect tier by a new dated decision. | Literature-curated, versioned, and traceable to primary publications; see E3 source in [`docs/gene_set_methods.md`](gene_set_methods.md). |
| 2026-11-12 | **Provisional:** require at least 2 supporting publications per gene–chemical pair in the environmental set. | Keeps a heavily studied exposure manageable; confirm or revise after E10 category counts, with a new dated entry. |
| 2026-11-12 | Expression-only genes named in the group research notes go into a separate, clearly labeled literature-context list used only for interpretation, never as risk or environment evidence. | [`research/Research for Group BINF 2111 Project .pdf`](../research/) reports disease-versus-control expression, not inherited or exposure evidence. |

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
