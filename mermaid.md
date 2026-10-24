# Complete Splice Girl project map

This diagram shows how the biological question, gene research, RNA-seq
preparation, Rust software, Hodge mathematics, statistics, interpretation, and
final outputs connect. It also shows which layers are already working and which
ones still depend on real-data processing or reviewed research.

**Status date:** October 24, 2026

| Color | Meaning |
|---|---|
| Green | Implemented, verified, or frozen |
| Yellow | Active work or partly complete |
| Gray | Waiting for an earlier dependency |
| Blue | Planned later work |
| Purple | Optional or alternate path |
| Red | Interpretation boundary that must remain visible |

```mermaid
flowchart TB
    classDef complete fill:#e8f3e8,stroke:#2f6b3b,color:#153a1d,stroke-width:1.5px
    classDef active fill:#fff3cd,stroke:#8a6d1d,color:#4d3b08,stroke-width:1.5px
    classDef waiting fill:#eeeeee,stroke:#777777,color:#333333,stroke-width:1.2px
    classDef planned fill:#e8f0fa,stroke:#406a9b,color:#173451,stroke-width:1.2px
    classDef optional fill:#f1e8f7,stroke:#76508d,color:#3d2350,stroke-width:1.2px
    classDef caution fill:#fbe9e9,stroke:#9a3d3d,color:#5c1f1f,stroke-width:1.5px

    subgraph L0["LAYER 0 - Question, scope, and frozen rules"]
        direction LR
        Q["Research question:<br/>Does splice-graph structure differ between<br/>endometriosis and control endometrium?"]
        TISSUE["Comparison boundary:<br/>same tissue and proliferative phase"]
        CONTRACT["Input contract v0.5:<br/>sample_id, gene_id, from, to, count"]
        CONVENTIONS["Frozen graph rules:<br/>1-based closed coordinates;<br/>full exon intervals; transcript direction"]
        SAMPLE_RULE["Samples remain separate:<br/>one biological sample is never replaced<br/>by pooled junctions"]
        REPRO_RULE["Reproducibility rules:<br/>versions, commands, checksums,<br/>settings, exclusions, run ID"]
        Q --> TISSUE
        TISSUE --> SAMPLE_RULE
        CONTRACT --> CONVENTIONS
        CONVENTIONS --> REPRO_RULE
    end

    subgraph L1["LAYER 1 - Candidate-gene evidence"]
        direction LR
        G_RULES["Define genetic-risk<br/>inclusion and exclusion rules"]
        G_SOURCES["Review genetic sources:<br/>variants, loci, inherited susceptibility"]
        G_ROWS["Genetic evidence table:<br/>stable ID, symbol, evidence, source"]
        E_RULES["Define environmental-response<br/>scope and evidence rules"]
        E_SOURCES["Review environmental sources:<br/>response context and exposure evidence"]
        E_ROWS["Environmental evidence table:<br/>stable ID, symbol, context, source"]
        ID_MAP["Normalize gene identifiers<br/>and audit ambiguous mappings"]
        RISK_ONLY["Risk-only genes"]
        ENV_ONLY["Environment-only genes"]
        OVERLAP["Risk + environment overlap"]
        CANDIDATES["Reviewed candidate-gene set<br/>with provenance for every gene"]

        G_RULES --> G_SOURCES --> G_ROWS --> ID_MAP
        E_RULES --> E_SOURCES --> E_ROWS --> ID_MAP
        ID_MAP --> RISK_ONLY --> CANDIDATES
        ID_MAP --> ENV_ONLY --> CANDIDATES
        ID_MAP --> OVERLAP --> CANDIDATES
    end

    subgraph L2["LAYER 2 - Study samples and metadata"]
        direction LR
        GEO["Public study:<br/>GEO GSE315857"]
        CTRL["4 control samples:<br/>SRR36710390-SRR36710387"]
        CASE["4 endometriosis samples:<br/>SRR36710386-SRR36710383"]
        SRA["NCBI SRA paired-end runs"]
        FASTQ["16 raw FASTQ files:<br/>forward + reverse for 8 samples"]
        META["Sample metadata:<br/>condition, tissue, subject, accession"]
        MANIFEST["data/metadata/samples.tsv:<br/>one traceable row per sample"]
        SUITABILITY["Suitability review:<br/>tissue, phase, pairing, replication,<br/>batch information, exclusions"]

        GEO --> CTRL --> SRA
        GEO --> CASE --> SRA
        SRA --> FASTQ
        GEO --> META --> MANIFEST --> SUITABILITY
    end

    subgraph L3["LAYER 3 - Exact human reference"]
        direction LR
        FASTA["GENCODE v50 GRCh38.p14<br/>primary-assembly FASTA"]
        GTF["Matching GENCODE v50<br/>comprehensive primary-assembly GTF"]
        REF_CHECK["Checksum and sequence-name<br/>compatibility checks"]
        HISAT_INDEX["HISAT2 genome index"]
        SPLICE_HINTS["Known splice-site hints<br/>from the matching GTF"]
        EXON_MAP["Annotation map:<br/>intron coordinates -> unique gene<br/>and adjacent exon intervals"]

        FASTA --> REF_CHECK
        GTF --> REF_CHECK
        REF_CHECK --> HISAT_INDEX
        GTF --> SPLICE_HINTS
        GTF --> EXON_MAP
    end

    subgraph L4["LAYER 4 - RNA-seq preprocessing, one sample at a time"]
        direction TB
        RAW_QC["FastQC on raw read pairs"]
        MULTIQC["MultiQC summary across samples"]
        CLEAVER["Cleaver-X / fastp:<br/>adapter trimming, filtering, JSON QC"]
        PHIX_DECISION{"PhiX contamination<br/>supported by experiment or QC?"}
        PHIX["Optional Bowtie2 PhiX removal"]
        CLEAN_READS["Clean paired FASTQ reads"]
        HISAT["HISAT2 splice-aware alignment<br/>with one frozen command"]
        BAM["Coordinate-sorted SAM/BAM<br/>kept per sample"]
        ALIGN_QC["Alignment QC:<br/>input reads, mapped reads,<br/>unique mapping, concordant pairs"]
        JUNCTIONS["Extract junction-spanning reads<br/>with unique-mapping support"]
        ANNOTATE["Map each junction to one unique<br/>gene and adjacent exon pair"]
        FILTER["Exclude unannotated, ambiguous,<br/>unknown-strand, and low-support rows"]
        STANDARD["Standard Splice Girl TSV:<br/>sample_id, gene_id, from, to, count"]
        CONVERT_QC["Conversion QC:<br/>kept and excluded junction counts"]
        STAR["Legacy STAR prototype<br/>for development only"]

        FASTQ --> RAW_QC --> MULTIQC
        FASTQ --> CLEAVER --> PHIX_DECISION
        PHIX_DECISION -- "yes" --> PHIX --> CLEAN_READS
        PHIX_DECISION -- "no" --> CLEAN_READS
        CLEAN_READS --> HISAT --> BAM --> ALIGN_QC --> JUNCTIONS
        JUNCTIONS --> ANNOTATE --> FILTER --> STANDARD
        FILTER --> CONVERT_QC
        CLEAN_READS -. "legacy prototype" .-> STAR -.-> JUNCTIONS
    end

    subgraph L5["LAYER 5 - Splice Girl input and safety checks"]
        direction LR
        CLI["splice-girl CLI:<br/>analyze, gene, visualize, --start"]
        TSV_PARSER["Junction parser validates:<br/>file, header, fields, IDs, counts"]
        MANIFEST_PARSER["Manifest parser validates:<br/>sample_id and case/control label"]
        ERROR_STOP["Safe failure:<br/>filename, line, column, bad value,<br/>likely correction, nonzero exit"]
        GROUP["Group rows by<br/>(sample_id, gene_id)"]
        NO_POOL["Preserve each sample;<br/>never pool patients"]
        OUTPUT_GUARD["Refuse to overwrite an existing<br/>result folder or SVG"]

        CLI --> TSV_PARSER
        CLI --> MANIFEST_PARSER
        TSV_PARSER -- "invalid" --> ERROR_STOP
        MANIFEST_PARSER -- "invalid" --> ERROR_STOP
        TSV_PARSER -- "valid" --> GROUP --> NO_POOL
        CLI --> OUTPUT_GUARD
    end

    subgraph L6["LAYER 6 - Deterministic splice-graph construction"]
        direction LR
        GRAPH["Build one directed graph<br/>for one sample and one gene"]
        NODES["Nodes:<br/>unique annotated exon intervals"]
        EDGES["Edges:<br/>observed directed splice junctions"]
        COUNTS["Raw edge weights:<br/>junction-spanning read counts"]
        DUPLICATES["Reject duplicate directed edges,<br/>self-loops, and zero-total groups"]
        ORDER["Sort node and edge IDs<br/>for deterministic order"]
        NORM["Normalize within sample/gene:<br/>F_e = count_e / total support"]
        GRAPH_OUT["Validated weighted graph<br/>and normalized edge signal F"]

        GROUP --> GRAPH
        GRAPH --> NODES
        GRAPH --> EDGES
        GRAPH --> COUNTS
        NODES --> ORDER
        EDGES --> ORDER
        COUNTS --> DUPLICATES
        DUPLICATES --> NORM
        ORDER --> GRAPH_OUT
        NORM --> GRAPH_OUT
    end

    subgraph L7["LAYER 7 - Hodge decomposition on each graph"]
        direction LR
        INCIDENCE["Incidence matrix B1:<br/>-1 at source, +1 at destination"]
        POTENTIAL["Fit node potential phi:<br/>minimize ||F - B1^T phi||^2"]
        SOLVER["Dense solver:<br/>one pinned node per component,<br/>partial pivoting, f64"]
        GRAD["Gradient component:<br/>F_grad = B1^T phi"]
        CYCLE["Cycle component:<br/>F_cycle = F - F_grad"]
        ENERGY["Energies:<br/>||F||^2, ||F_grad||^2, ||F_cycle||^2"]
        FRACTION["Cycle fraction:<br/>||F_cycle||^2 / ||F||^2"]
        RECON["Check 1:<br/>F approximately F_grad + F_cycle"]
        KERNEL["Check 2:<br/>B1 F_cycle approximately 0"]
        ORTHO["Check 3:<br/>F_grad dot F_cycle approximately 0"]
        MATH_GATE{"All mathematical checks<br/>pass at tolerance 1e-9?"}
        GROUP_RESULT["Per-sample, per-gene result:<br/>metrics + signed edge components"]
        MATH_STOP["Stop and report the<br/>sample/gene math failure"]

        GRAPH_OUT --> INCIDENCE --> POTENTIAL --> SOLVER --> GRAD
        GRAPH_OUT --> CYCLE
        GRAD --> CYCLE
        GRAPH_OUT --> ENERGY
        GRAD --> ENERGY
        CYCLE --> ENERGY --> FRACTION
        GRAPH_OUT --> RECON
        GRAD --> RECON
        CYCLE --> RECON
        INCIDENCE --> KERNEL
        CYCLE --> KERNEL
        GRAD --> ORTHO
        CYCLE --> ORTHO
        RECON --> MATH_GATE
        KERNEL --> MATH_GATE
        ORTHO --> MATH_GATE
        MATH_GATE -- "yes" --> GROUP_RESULT
        MATH_GATE -- "no" --> MATH_STOP
    end

    subgraph L8["LAYER 8 - Synthetic validation and software testing"]
        direction LR
        LINEAR["Linear chain:<br/>expected cycle fraction 0"]
        BRANCH["Branching tree:<br/>expected cycle fraction 0"]
        SKIP["Exon-skipping triangle:<br/>expected cycle fraction 1/9"]
        PURE_CYCLE["Equal directed cycle:<br/>expected cycle fraction 1"]
        NOISE["Fixed-noise example:<br/>reviewed exact expectation"]
        SWEEP["Back-splice weight sweep:<br/>five reviewed answers"]
        EDGE_CASES["Zero support, duplicate edge,<br/>self-loop, disconnected, row order"]
        ERROR_TESTS["CLI error suite:<br/>nine real command failures"]
        ACCEPTANCE["Automated acceptance suite:<br/>math, parser, CLI, output, SVG"]

        LINEAR --> ACCEPTANCE
        BRANCH --> ACCEPTANCE
        SKIP --> ACCEPTANCE
        PURE_CYCLE --> ACCEPTANCE
        NOISE --> ACCEPTANCE
        SWEEP --> ACCEPTANCE
        EDGE_CASES --> ACCEPTANCE
        ERROR_TESTS --> ACCEPTANCE
    end

    subgraph L9["LAYER 9 - Sample-preserving case/control comparison"]
        direction LR
        JOIN_MANIFEST["Join per-sample gene results<br/>to the case/control manifest"]
        INCLUDE["Include:<br/>sample exists in both inputs<br/>and mathematical checks pass"]
        EXCLUDE["Exclude with reason:<br/>missing metadata, missing gene data,<br/>or failed checks"]
        GROUP_GATE{"At least one eligible case<br/>and one eligible control?"}
        SAMPLE_TABLE["One row per sample/gene:<br/>condition, cycle fraction,<br/>status, exclusion reason"]
        CASE_VALUES["Endometriosis sample values"]
        CONTROL_VALUES["Control sample values"]

        GROUP_RESULT --> JOIN_MANIFEST
        MANIFEST_PARSER --> JOIN_MANIFEST
        JOIN_MANIFEST --> INCLUDE --> GROUP_GATE
        JOIN_MANIFEST --> EXCLUDE --> SAMPLE_TABLE
        GROUP_GATE -- "yes" --> SAMPLE_TABLE
        GROUP_GATE -- "no" --> ERROR_STOP
        SAMPLE_TABLE --> CASE_VALUES
        SAMPLE_TABLE --> CONTROL_VALUES
    end

    subgraph L10["LAYER 10 - Statistical design, baselines, and ranking"]
        direction LR
        STAT_SPEC["Freeze statistics specification:<br/>biological unit, eligibility,<br/>missing values, graph comparability"]
        EFFECT["Primary effect:<br/>case summary - control summary"]
        PERMUTE["Permutation test only when<br/>sample design permits label exchange"]
        SEED["Saved random seed and<br/>exact or sampled permutations"]
        FDR["Multiple-testing correction:<br/>defined gene family + q-values"]
        BASE_JUNCTION["Baseline 1:<br/>junction-usage differences"]
        BASE_RANK["Baseline 2:<br/>graph cycle rank"]
        BASE_EXPR["Optional baseline 3:<br/>expression, if compatible"]
        SENSITIVITY["Sensitivity checks:<br/>approved alternate thresholds"]
        RANK["Deterministic ranking:<br/>effect, inference status,<br/>candidate category, junction explanation"]
        STAT_REPORT["Statistics report:<br/>effects, sample counts, p/q status,<br/>baselines, settings, limitations"]

        CASE_VALUES --> EFFECT
        CONTROL_VALUES --> EFFECT
        STAT_SPEC --> EFFECT
        STAT_SPEC --> PERMUTE --> SEED --> FDR
        EFFECT --> PERMUTE
        SAMPLE_TABLE --> BASE_JUNCTION
        SAMPLE_TABLE --> BASE_RANK
        SAMPLE_TABLE --> BASE_EXPR
        EFFECT --> SENSITIVITY
        FDR --> RANK
        EFFECT --> RANK
        BASE_JUNCTION --> RANK
        BASE_RANK --> RANK
        BASE_EXPR --> RANK
        SENSITIVITY --> RANK --> STAT_REPORT
    end

    subgraph L11["LAYER 11 - Evidence-linked biological interpretation"]
        direction LR
        PROVENANCE_JOIN["Join results to candidate-gene<br/>provenance by stable gene ID"]
        CATEGORY_COMPARE["Compare risk-only,<br/>environment-only, and overlap groups"]
        JUNCTION_EXPLAIN["Trace each result back to<br/>samples, graph edges, and junction values"]
        NULL_OK["A null or no-advantage result<br/>is still a valid result"]
        INTERPRET["Reviewed gene summaries<br/>with evidence and limitations"]
        CIRC_LIMIT["Cycle-space signal<br/>does not prove circular RNA"]
        CAUSE_LIMIT["Group difference<br/>does not prove disease causation"]
        EXPOSURE_LIMIT["Environmental-response gene<br/>does not prove exposure causation"]
        STUDY_LIMIT["Report tissue, sample-size,<br/>batch, alignment, and coverage limits"]

        RANK --> PROVENANCE_JOIN
        CANDIDATES --> PROVENANCE_JOIN
        PROVENANCE_JOIN --> CATEGORY_COMPARE --> INTERPRET
        STAT_REPORT --> JUNCTION_EXPLAIN --> INTERPRET
        NULL_OK --> INTERPRET
        CIRC_LIMIT --> INTERPRET
        CAUSE_LIMIT --> INTERPRET
        EXPOSURE_LIMIT --> INTERPRET
        STUDY_LIMIT --> INTERPRET
    end

    subgraph L12["LAYER 12 - Reproducible outputs and audience-facing products"]
        direction LR
        METRICS_TSV["sample_metrics.tsv:<br/>one row per sample/gene"]
        COMPONENTS_TSV["junction_components.tsv:<br/>observed, gradient, cycle per edge"]
        COMPARISONS_TSV["comparisons.tsv:<br/>included and excluded samples"]
        RUN_JSON["run_manifest.json:<br/>command, version, checksums,<br/>settings, exclusions, failures"]
        THREE_PANEL["Generated academic SVG:<br/>observed | gradient | cycle"]
        DEMO["Tiny synthetic CLI demo<br/>with known answers"]
        README_OUT["Beginner README,<br/>dictionary, troubleshooting"]
        REPORT["Final methods + results<br/>with versioned tables"]
        SLIDES["Presentation figures,<br/>speaking plan, live demo, backup"]
        RELEASE["Reusable Splice Girl release:<br/>source, tests, examples, version"]

        GROUP_RESULT --> METRICS_TSV
        GROUP_RESULT --> COMPONENTS_TSV
        SAMPLE_TABLE --> COMPARISONS_TSV
        REPRO_RULE --> RUN_JSON
        METRICS_TSV --> RUN_JSON
        COMPONENTS_TSV --> RUN_JSON
        COMPARISONS_TSV --> RUN_JSON
        GROUP_RESULT --> THREE_PANEL
        ACCEPTANCE --> DEMO
        ERROR_TESTS --> README_OUT
        STAT_REPORT --> REPORT
        INTERPRET --> REPORT
        THREE_PANEL --> SLIDES
        INTERPRET --> SLIDES
        DEMO --> SLIDES
        ACCEPTANCE --> RELEASE
        README_OUT --> RELEASE
    end

    %% Cross-layer dependencies that join separate workstreams.
    SAMPLE_RULE --> MANIFEST
    MANIFEST --> MANIFEST_PARSER
    SUITABILITY --> RAW_QC
    HISAT_INDEX --> HISAT
    SPLICE_HINTS --> HISAT
    EXON_MAP --> ANNOTATE
    CONVENTIONS --> EXON_MAP
    CONTRACT --> STANDARD
    STANDARD --> CLI
    CONVERT_QC --> CLI
    CANDIDATES --> STAT_SPEC
    TISSUE --> STAT_SPEC
    ACCEPTANCE --> RELEASE

    class CONTRACT,CONVENTIONS,SAMPLE_RULE,REPRO_RULE,MANIFEST,FASTA,GTF,REF_CHECK,EXON_MAP,TSV_PARSER,MANIFEST_PARSER,ERROR_STOP,GROUP,NO_POOL,OUTPUT_GUARD,GRAPH,NODES,EDGES,COUNTS,DUPLICATES,ORDER,NORM,GRAPH_OUT,INCIDENCE,POTENTIAL,SOLVER,GRAD,CYCLE,ENERGY,FRACTION,RECON,KERNEL,ORTHO,MATH_GATE,GROUP_RESULT,MATH_STOP,LINEAR,BRANCH,SKIP,PURE_CYCLE,NOISE,SWEEP,EDGE_CASES,ERROR_TESTS,ACCEPTANCE,JOIN_MANIFEST,INCLUDE,EXCLUDE,GROUP_GATE,SAMPLE_TABLE,CASE_VALUES,CONTROL_VALUES,METRICS_TSV,COMPONENTS_TSV,COMPARISONS_TSV,RUN_JSON,THREE_PANEL,DEMO,README_OUT complete
    class GEO,CTRL,CASE,SRA,FASTQ,META,SUITABILITY,RAW_QC,MULTIQC,CLEAVER,CLEAN_READS,HISAT_INDEX,SPLICE_HINTS active
    class G_RULES,G_SOURCES,G_ROWS,E_RULES,E_SOURCES,E_ROWS,ID_MAP,RISK_ONLY,ENV_ONLY,OVERLAP,CANDIDATES,HISAT,BAM,ALIGN_QC,JUNCTIONS,ANNOTATE,FILTER,STANDARD,CONVERT_QC,STAT_SPEC,EFFECT,PERMUTE,SEED,FDR,BASE_JUNCTION,BASE_RANK,SENSITIVITY,RANK,STAT_REPORT,PROVENANCE_JOIN,CATEGORY_COMPARE,JUNCTION_EXPLAIN,NULL_OK,INTERPRET,REPORT,SLIDES,RELEASE waiting
    class BASE_EXPR planned
    class PHIX_DECISION,PHIX,STAR optional
    class CIRC_LIMIT,CAUSE_LIMIT,EXPOSURE_LIMIT,STUDY_LIMIT caution
```

## The main handoffs

1. **Research to analysis:** reviewed gene evidence becomes the candidate set
   used for filtering, annotation, and later interpretation.
2. **Raw reads to Splice Girl:** each sample becomes one checked five-column
   junction table; raw FASTQ, BAM, and reference files do not enter the Rust
   graph engine directly.
3. **Graph to mathematics:** deterministic nodes, edges, and normalized values
   become the incidence matrix and Hodge decomposition.
4. **Mathematics to comparison:** every sample keeps its identity and produces
   one gene-level value before case/control comparison.
5. **Statistics to interpretation:** effects and any valid inference are joined
   back to gene provenance and individual junction explanations.
6. **Every layer to reproducibility:** versions, commands, checksums, settings,
   exclusions, and failures are preserved with the output.

## Detailed documents

- [Project priorities](priorities.md)
- [Full task guide](delegations_and_tasks/tasks.md)
- [Preprocessing methods](docs/data_methods.md)
- [Input and CLI contract](docs/data_contract.md)
- [Hodge mathematics](docs/math_spec.md)
- [Runnable preprocessing notes](preprocessing/README.md)
- [Synthetic demonstration](examples/demo/README.md)
