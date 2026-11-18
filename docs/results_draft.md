# Results draft — technical progress only

Updated November 17, 2026. This is usable report text about completed work;
the full eight-sample biological results are still pending.

## Software and evidence inputs

Splice Girl implements input validation, graph construction, within-gene
normalization, graph Hodge decomposition, descriptive comparisons, explicit
exact permutation tests, candidate-set tests, BH correction, baseline exports,
support/coverage sensitivity runs and offline reports. The candidate table
contains 661 versioned gene IDs: 96 risk-only, 558 environment-only and seven
overlap genes. Those labels come from evidence curation and do not describe
the sequenced donors' exposure histories. Independent evidence sign-offs
remain listed in [spot_checks.md](spot_checks.md).

## Human-reference technical pilot

The first one million synchronized cleaned pairs from control SRR36710387
were aligned against the verified GENCODE v50 human primary assembly. HISAT2
reported 97.23% overall alignment. The annotation converter exported 32,575
junction rows. A predeclared coverage/size rule selected CDH11; at minimum
support three its graph retained six nodes, four edges and 103 read-junction
observations. Independent BAM recounting matched the six base edges and 106
observations, and the Rust projection agreed with an independent SVD to
approximately 4.6e-11 per exported component. See the
[versioned pilot bundle](../data/pilots/SRR36710387_prefix1m_nov17/README.md).

The selected graph had two connected components and cycle rank zero.
Its cycle fraction was effectively zero at support thresholds 1, 3, 5 and 10.
This is an expected mathematical result for this graph and a successful
technical check. It does not demonstrate a nonzero cycle component in real
data or establish absence of alternative splicing or circular RNA.

## Assignment quality and execution scale

Of 587,089 assignment-stage observations, 152,907 were uniquely assigned,
432,599 were ambiguous, and 1,583 were unannotated. The detailed audit found
that 409,613 of the ambiguous observations identified a single gene but
multiple full exon pairs. These exclusions were retained under the frozen
mapping policy. Library orientation and assignment coverage require review
before expanding. The [audit figures](../data/pilots/SRR36710387_prefix1m_nov17/audit_figures/README.md)
show these diagnostics without presenting them as biological differences.

Across all 243 candidate genes observed in this prefix, the Rust analysis and
table export processed 1,388 junction rows in a median of 11.2 ms over five
runs, with a maximum measured child-process peak of about 4.0 MiB. All groups
passed internal math checks on every run. The largest observed graph had 66
nodes. This measures analysis/export on one prefix, not alignment, report
rendering or full-cohort runtime. The
[benchmark record](../data/pilots/SRR36710387_prefix1m_nov17/candidate_benchmark/benchmark.json)
retains per-run measurements and hashes. The 418 unobserved candidates were
not assigned zero-valued results.

## Biological comparisons still to report

No case/control effect, disease-associated gene, significant gene set, or
advantage over simpler baselines is claimed from this one-control pilot.
The final results need eight reviewed sample inputs, comparable gene graphs,
eligibility/exclusion tables, the configured primary candidate-set test,
exploratory gene-level results, same-sample baselines and sensitivity analysis.
Tissue/cell composition, the small cohort, incomplete batch information and
the lack of measured exposures must accompany any later interpretation.
