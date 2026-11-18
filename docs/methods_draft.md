# Methods draft

Draft assembled November 17, 2026. Completed steps are described in the past
tense; planned cohort analysis is explicitly conditional. This does not replace
the approved proposal or the controlling specifications linked below.

## Study inputs and evidence

The selected public bulk RNA-seq cohort, GSE315857, contains four endometriosis
and four control samples of proliferative-phase eutopic endometrium. Accessions,
condition labels and internal sample/subject IDs are retained in the sample
manifest. The [dataset review](dataset_review.md) records unresolved library
methods discrepancies, incomplete batch information and the need to confirm
donor independence. Patient exposure was not measured in the available manifest.

Genetic-risk and environmental-response evidence were normalized to versioned
GENCODE identifiers. The candidate table retains the source-derived risk,
environment and overlap labels; its current categories are 96 risk-only,
558 environment-only and seven overlap genes. Inclusion rules, evidence
provenance and remaining independent checks are in
[gene_set_methods.md](gene_set_methods.md) and [spot_checks.md](spot_checks.md).

## Reference, processing and technical validation

The verified GENCODE v50 primary-assembly GRCh38.p14 FASTA and comprehensive
GTF were used to construct the HISAT2 index, annotated splice-site hints and
exon-pair map. Exact source/output hashes and versions are retained. The
implemented adapter uses the lab's Cleaver executable for cleaning and BAM
operations, HISAT2 for splice-aware alignment, and the checked junction
converter. Its outputs retain sample IDs and stage/QC manifests.

One control, SRR36710387, had an existing Cleaver-cleaned paired dataset.
The first one million synchronized cleaned pairs were selected by a recorded
deterministic-prefix rule and aligned using HISAT2 2.2.3. This technical pilot
used exploratory unstranded alignment, eight-base anchors and unique primary
read alignments. Junctions were assigned only when the annotation map
identified one exact exon pair; ambiguous and unannotated observations were
excluded and counted. No alternative pair was selected to rescue support.
The [data methods](data_methods.md), [pilot guide](real_pilot.md) and
[assignment audit](../data/pilots/SRR36710387_prefix1m_nov17/assignment_audit/README.md)
give the coordinate conventions, commands and exclusions. These exploratory
settings have not been approved as final settings for all samples.

A coverage/size rule, recorded before inspecting Hodge values, selected CDH11
for a readable real-input example. Counts were independently reconstructed
from the indexed BAM. The Rust result was compared edge by edge against a
NumPy SVD projection at support thresholds 1, 3, 5 and 10. This is independent
verification of one selected graph's calculation; it does not validate every
observed gene's biology or demonstrate a nonzero real-data cycle signal.

## Graph and Hodge calculation

Each graph contains one sample and one gene. Nodes identify annotated exon
intervals; directed edges follow transcript order. Each edge count is divided
by the total support for that sample/gene to obtain relative junction usage.
The incidence matrix has minus one at an edge's source and plus one at its
target. Least squares projects the edge signal onto node-potential differences.
Subtracting this gradient projection gives the cycle-space component.

The cycle fraction is the squared length of that component divided by the
squared length of the input signal. Reconstruction, zero divergence of the
residual and orthogonality are checked numerically. Zero-total groups are
rejected rather than assigned a fabricated fraction. The exact equations,
tolerances, synthetic expected answers and graph-cycle interpretation are in
[math_spec.md](math_spec.md). The Rust implementation has no external crates.

## Planned cohort comparison

After the remaining inputs and preprocessing rules are reviewed, comparable
per-sample graphs will be analyzed for the fixed candidate family. The
reported effect is mean case cycle fraction minus mean control cycle fraction.
Missing or excluded sample/gene values remain missing. The frozen primary
test uses the candidate-set mean absolute gene effect and shared exhaustive
label assignments; with four independent samples per group there are 70
assignments. Gene-level exact tests and BH-adjusted values are exploratory.
Validity depends on the design's exchangeability assumptions, not merely the
existence of eight filenames.

The software exports same-sample support/topology baselines, optional rankings
and an explicitly configured support/coverage sensitivity grid. Applying
these to the full cohort, recording eligibility/exclusions and interpreting
the resulting biological effects are pending. The authoritative statistical
policy is [statistics_spec.md](statistics_spec.md).

## Reproducibility and limits

The software retains input hashes, parameters, status fields, signed junction
components, result tables and offline graph reports. See
[software_validation.md](software_validation.md) for compiler/test records and
the actual-prefix benchmark. Large reads and alignment files stay outside Git.
Graph cycle space alone does not identify circular RNA; evidence membership
does not establish exposure, causation or a gene–environment interaction.
