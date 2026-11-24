# GSE315857 suitability review

Technical review, November 17, 2026. **Closed 2026-11-24:** the dataset is
accepted for the study; see [final disposition](#final-disposition-2026-11-24).
Per-run checks still apply to every sample as it is processed.

| Question | Evidence and assessment |
|---|---|
| What is compared? | Four endometriosis and four control samples. The [sample manifest](../data/metadata/samples.tsv) links each run to its GEO record and condition. |
| Are the tissues comparable? | Both groups use eutopic endometrium in the proliferative phase. This avoids comparing a lesion directly with a different tissue. It does not control all cell-composition or clinical differences. |
| What are the controls? | The recorded definition is laparoscopically negative for endometriosis. |
| Are these paired subjects? | The current design is independent groups. Internal subject IDs identify the eight samples; they are not published patient identifiers or proof of matching. Confirm there are no repeated donors before final inference. |
| What files exist? | Saved raw FastQC/MultiQC summaries for all eight runs; one complete cleaned control locally; one million-pair aligned prefix from that control. Other FASTQ pairs still need an accessible storage location. |
| How much sequence? | The saved QC tables contain 545,210,220 pairs across the eight runs. Both mates agree for every run; median read length is 100 bases. These are raw counts, before cleaning. |
| Is library information settled? | Orientation yes: all eight runs are RF ([2026-11-18 check](../data/strandedness_2026-11-18/README.md)). The kit description still differs between sources; see below. |
| Are batch/covariates complete? | Sequencing batch is not reported in the local manifest. Eight samples cannot support a large covariate model. Check donor independence and clinical metadata before interpreting diagnosis effects. |
| Can this test exposure effects? | No patient exposure measurements are present in the study manifest. Environmental evidence prioritizes genes; it does not label these donors as exposed. |
| Can this test the proposed graph method? | Yes (accepted 2026-11-24; see the [final disposition](#final-disposition-2026-11-24)): paired bulk RNA-seq, matched tissue/phase, accessible run identifiers and an approved human-reference technical pilot. Each sample must still pass its own per-run checks. |

## A methods discrepancy to resolve

The [GEO record for the pilot sample](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9438615)
describes rRNA depletion and fastp/HISAT2/StringTie processing. The
[associated YWHAZ publication](https://academic.oup.com/reproduction/article/171/4/xaag038/8541717)
describes KAPA mRNA library preparation, 50-base paired reads, and STAR/edgeR
processing. The saved raw QC summaries have a 100-base median read length.
These descriptions do not agree sufficiently to infer the original library
orientation or reproduce the authors' processing by assumption. Splice Girl's
own chosen reference, pipeline and observed file lengths are recorded separately.
Do not substitute one source's protocol for another without resolving this.

## What the alignment itself tells us

The [assignment audit](../data/pilots/SRR36710387_prefix1m_nov17/assignment_audit/README.md)
counts each eligible spliced alignment once for orientation. Read 1 is
antisense to the inferred transcript in 99.7205% of these alignments; read 2
is sense in 99.7221%. This is consistent with an RF library as defined by the
[HISAT2 manual](https://daehwankimlab.github.io/hisat2/manual/).
It is technical evidence from one prefix, not an independently confirmed
library protocol for all eight samples. The exploratory unstranded pilot
has not been relabeled or silently rerun as a final RF analysis.

## What remains before accepting all eight

1. ~~Resolve library orientation and inspect it across samples.~~ Done
   2026-11-18: the reads of all eight runs show RF orientation
   ([check](../data/strandedness_2026-11-18/README.md)). The written sources
   still disagree on kit and read length; that affects only how the authors'
   own processing is described, not this project's settings.
2. ~~Review the exact-exon-pair ambiguity loss and retain one documented rule
   across the cohort.~~ Done 2026-11-18: the
   [shared-splice-site node rule](data_methods.md#graph-nodes).
3. ~~Set anchor, support, gene coverage and missing-value rules before
   comparing case/control effects; check graph comparability after applying
   them.~~ Done: thresholds approved 2026-11-22 and confirmed 2026-11-24
   ([proposal](threshold_proposal.md)). One shared graph per gene makes every
   covered sample's graph identical, so graph comparability holds by
   construction.
4. Verify raw/cleaned checksums, donor identity, QC and sample exclusions for
   every run. Resolved as a procedure; see below. A technical prefix is not a
   substitute for a full biological sample.

## Final disposition (2026-11-24)

**GSE315857 is accepted for the study's primary analysis**, with the points
below recorded as settled, handled per run, or accepted limits.

| Point | Disposition |
|---|---|
| Design and tissue | Settled: 4 endometriosis and 4 laparoscopically negative controls, eutopic endometrium, proliferative phase, analyzed as independent groups. |
| Library orientation | Settled: RF in all eight runs by direct measurement (2026-11-18), rechecked after alignment for every run (95% rule). |
| Kit and read-length conflict between the GEO record and the paper | **Not resolvable from the sources, and not needed.** It concerns how the original authors processed their data. This project's orientation was measured from the reads and the read length (median 100 bases) from the files themselves; no setting relies on either written description. The methods state the disagreement. |
| Raw checksums | Handled per run: `hisat2_pipeline.py` records the SHA-256 of every input read file and output in each run's manifest, and `assemble_runs.py` refuses a run whose recorded outputs changed. Complete for each sample when its run completes. |
| QC and sample exclusions | Handled per run: raw FastQC/MultiQC exists for all eight runs; cleaning QC, alignment summaries and the strand check are written by each run; a sample held by the strand check is excluded with its reason and the analysis records it. |
| Repeated donors | **Accepted assumption.** The eight runs are eight distinct GEO samples (GSM records) with no shared subject identifier and no paired design described. Public records cannot prove that no donor contributed twice; the analysis treats each run as one subject (decision 2026-10-29), and the methods state this as an assumption. |
| Batch, clinical covariates, exposure | Accepted limits: sequencing batch is not reported, eight samples cannot support covariate models, and no exposure was measured. These are stated as limits, not corrected for. |
| Thresholds and graph comparability | Settled (item 3 above). |

What this does **not** cover: whether each individual run passes its own
checks. That is decided when the run is processed, by the recorded checks
above, not by this review.
