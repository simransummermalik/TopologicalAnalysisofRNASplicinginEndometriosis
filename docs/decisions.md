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
| 2026-11-17 | E6 context rules for CTD environmental evidence: exclude co-treatment rows (1,028), because a mixture response cannot be credited to the dioxin-like compound alone. Keep, labeled, rows where another agent alters the dioxin's effect (`dioxin_effect_modified`) and rows where the dioxin alters another agent's effect (`dioxin_modifies_other`). Count a row only when its own gene is the measured one. | E6 review in [`docs/gene_set_methods.md`](gene_set_methods.md); the context choice changes the ≥2-publication gene count only from 1,693 to 1,640. |
| 2026-11-17 | **Final environmental publication rule (replaces the provisional 2026-11-12 rule):** a gene–chemical pair needs at least 2 distinct publications, and at least 1 of them must not be a large screen. A large screen is a publication cited by 100 or more mapped evidence rows in the CTD extraction (17 publications in the 2026-10-29 release). | With the provisional rule, 64% of qualifying pairs rested only on genome-wide screens. The final rule keeps replication, removes screen-only support, and keeps every gene with endometrial evidence the provisional rule kept. Result: 565 environmental genes, 661 candidates, 7 overlap. Options and counts are in [`docs/gene_set_methods.md`](gene_set_methods.md); chosen before any real splicing results existed. |
| 2026-11-17 | GWAS Catalog risk rows: exclude cross-trait (pleiotropy) analyses (75 rows) and MTAG analyses (92 rows). Exclude Koller et al. analyses other than the cross-ancestry combined definition, GCST90841381 (139 rows), which applies the 2026-11-12 Koller decision to the Catalog copy. Keep endometriosis subtype analyses from other papers (stage 3 or 4; deep disease). Excluded rows stay in the table with a reason. | G1 requires an endometriosis association in a GWAS or meta-analysis. Pleiotropy tests a shared signal, and MTAG borrows power from other traits. Recorded in the G4–G5 extraction record in [`docs/gene_set_methods.md`](gene_set_methods.md). |
| 2026-11-17 | The Koller source is read from Supplementary Table 4, first section only (80 cross-ancestry loci), using the authors' "Nearest gene" column as the reported gene. A ":"-joined value names several genes and is treated as ambiguous. | The supplement became available from the publisher's static content server. The column is the source's own gene report, which G1 allows. |
| 2026-11-17 | G7 risk ambiguity rules. When the authors list several genes that all overlap the lead SNP, keep every listed gene, one row each, flagged `multi-gene locus`. Translate retired clone names through GENCODE v19 stable IDs. Exclude Catalog copies of Koller loci in favor of the Supplementary Table 4 row that carries the authors' gene. Use a single author-reported gene when the Catalog records one; otherwise exclude flanking-gene and no-gene rows. | G1 forbids choosing a gene by proximity. These rules only use what the source itself reports. See the G7 review in [`docs/gene_set_methods.md`](gene_set_methods.md). |
| 2026-11-17 | CTD rows get a paper-level tissue tag (`cited_paper:endometrium_uterus`, `cited_paper:other_female_reproductive`, or `cited_paper:other_or_unstated`) from PubMed MeSH headings and titles, with reviewed overrides. Tissue is never a reason to exclude a row. | CTD records no tissue; E1 does not restrict it. 98% of kept rows come from non-reproductive tissue, so the tag keeps that limitation visible. |
| 2026-11-17 | Statistical inference rules: every test is two-sided; a gene is testable only with all 4 cases and 4 controls eligible; the single confirmatory test is one gene-set test on `candidate_genes.tsv` at α = 0.05, using the mean over genes of the absolute case-minus-control cycle-fraction effect with all 70 relabelings applied to every gene at once. Risk, environmental, and overlap set tests and gene-level BH across testable candidates are exploratory. | With 4 versus 4 the smallest two-sided p is 2/70 ≈ 0.029: gene-level BH over hundreds of genes cannot pass, and any correction over 2 or more tests needs p ≤ 0.025. The absolute effect keeps genes changing in opposite directions from cancelling. Chosen before any real splicing results existed; see [inference rules](statistics_spec.md#inference-rules-decided-2026-11-17). |
| 2026-11-18 | Junction-to-node rule (replaces the node part of the 2026-10-22 convention): within a gene and strand, annotated exons that share a start or end coordinate are one node (`shared-splice-site`); a retained-intron exon may not merge two nodes that each hold a normal exon; node ID is `chrom:lowest-highest:strand`. Junctions joining the same node pair (alternative 5′/3′ sites) are summed into one edge. Junctions with both ends in one node are excluded and counted. The earlier `exon-interval` rule stays available to reproduce the 2026-11-17 pilot. | Under the earlier rule, 74% of the pilot's junction observations and 30% of annotated candidate junctions matched several exon pairs that differ only at ends a junction read cannot see. The new rule maps 98.5% of candidate junctions uniquely, keeps about 100% instead of 35% of same-gene candidate reads in the pilot, and raises candidate genes with any annotated cycle from 298 to 536. Splice-site hints are unchanged, so the existing HISAT2 index stays valid. Chosen before any case/control result existed; see the [node-rule review](data_methods.md#node-rule-review-2026-11-18). |
| 2026-11-18 | Library strandedness for GSE315857 is RF; the cohort alignment uses `--strandedness RF`. After alignment, a sample whose read 1 is antisense in less than 95% of eligible spliced alignments is held for review. | Read k-mers against nine reference mRNAs show read 1 antisense in 98.1–99.5% of matched reads for all eight runs, in both groups; the pilot alignment independently gave 99.7%. The written sources do not state the orientation. See the [strandedness check](../data/strandedness_2026-11-18/README.md). |

## Open decisions that still need a team answer

**Mathematical clarification, November 17:** the statistical-policy row above
retains its historical rationale, but two claims in that rationale were too
broad. BH discoveries are possible with balanced p-values of 2/70 when enough
genes share small p-values; Bonferroni/Holm and BH have different thresholds.
An unequal 3-versus-4 absolute-effect test can attain 1/35, not necessarily 2/35.
The added tests in `src/fdr.rs` and `src/permutation.rs` demonstrate both facts.
The chosen primary candidate-set test, two-sided alternative, four-per-group
eligibility and exploratory status of other tests are unchanged. The four-per-
group requirement is a consistent complete-cohort policy, not a claim that
significance is impossible for every smaller design.

Implementation note, 2026-11-17: `preprocessing/hisat2_pipeline.py` now implements
the HISAT2 path already named in the data methods; a tiny actual paired-read
Cleaver/HISAT2-to-Rust-report test passes. The STAR prototype is retained as
legacy code. This resolves the implementation mismatch without treating a
synthetic smoke test as approval of the full human pilot or its thresholds.

| Question | Why it matters | What must be recorded before implementation |
|---|---|---|
| Which exact support and coverage thresholds define an eligible junction? | Thresholds change which graphs and samples enter the analysis. | Threshold values, rationale, and a small sensitivity plan. |
| Does the real HISAT2 pilot pass the production-method review? | The implemented adapter now matches the data methods; full-cohort suitability remains unverified. Bowtie2 alone is not the splice-junction extractor. | Reviewed human alignment/QC, frozen tool versions and thresholds, and a checked real gene before D6 is marked done. |
| Which reviewed candidate-gene table is final? | Ranking and interpretation must use a frozen, source-linked set. | Stable gene IDs, evidence sources, and inclusion rules. |
| What are the course midpoint and final deadlines? | The task board cannot assign realistic due dates without them. | Dates, owner, and the required artifact at each checkpoint. |
| Where is the parent Overleaf proposal source? | The local roadmap section cannot be compiled alone. | Main `.tex` filename, folder, and compile command. |

An open decision is not a failure. It is safer to leave an item visible than to
hide an assumption inside code or a final result.

## File and output handling — November 17 clarification

Use the stable `SRR...` run IDs across read sheets, alignment outputs and
junction tables; retain the `case`/`control` labels from the study manifest.
Create a new named run folder for each execution rather than replacing prior
outputs. Large reads, indexes and alignments stay in ignored local/external
storage. Small review tables, plots and their checksum manifests may be kept
under `data/pilots/<sample_scope_date>/`; generated software review archives
stay under ignored `target/`. Preserve original source checksums and distinguish
full samples from deterministic prefixes. The user confirmed that there is
currently no external drive/server, so no full-cohort run is started on the
nearly full laptop. This records output handling, not a new scientific method.
