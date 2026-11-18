# GSE315857 suitability review

Technical review, November 17, 2026. This records what the available sources
support; it is not final approval of the cohort or library settings.

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
| Can this test the proposed graph method? | Provisionally yes: paired bulk RNA-seq, matched tissue/phase, accessible run identifiers and a successful human-reference technical pilot. Full eligibility still depends on cleaning, alignment, assignment coverage and graph comparability for every sample. |

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
3. Set anchor, support, gene coverage and missing-value rules before comparing
   case/control effects; check graph comparability after applying them.
4. Verify raw/cleaned checksums, donor identity, QC and sample exclusions for
   every run. A technical prefix is not a substitute for a full biological sample.

The dataset is a reasonable candidate, but these unresolved points prevent a
final biological result from being claimed today.
