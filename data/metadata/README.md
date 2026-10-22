# Sample metadata

The current Splice Girl analysis uses the eight bulk RNA-seq samples in
[`samples.tsv`](samples.tsv) from GEO study
[`GSE315857`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315857).
The table is the single sample manifest for the current analysis. Older
candidate-dataset tables remain in this folder for provenance, but they are not
inputs to the GSE315857 comparison.

## What one row means

Each row is one RNA-seq library from one biological subject. `sample_id` uses
the SRA run accession because that same value will be carried into future
junction rows. `condition` is deliberately limited to `case` or `control`, the
two values accepted by Splice Girl.

The publication states that the samples came from eight different women: four
with confirmed endometriosis and four laparoscopically negative controls. It
does not publish patient identifiers. The `subject_id` values in this manifest
are therefore internal labels, not original clinical IDs. They keep subject
identity separate from the sequencing-run ID without inventing patient data.

All eight samples are proliferative-phase eutopic endometrium, so the tissue is
matched for the planned case-control comparison. Tissue eligibility is marked
`yes`, while final analysis eligibility remains pending until preprocessing,
alignment, and junction QC are complete.

## Processing status in the table

- Raw FastQC and MultiQC are complete for all eight runs. The group reports are
  stored under [`data/galaxy_data/october_20`](../galaxy_data/october_20).
- `SRR36710387` is the completed Cleaver pilot. Its cleaned reads passed the
  recorded post-cleaning checks and are ready for alignment.
- Cleaver is not yet recorded as run for the other seven samples.
- No project alignment or junction table exists yet. These fields say
  `not_started`, `not_generated`, or `pending_D6` instead of pretending the
  missing work is complete.
- The project genome build and annotation remain `pending_D5` until the team
  freezes one matching reference pair. The source study reports hg38, but that
  does not replace the project reference decision.

The GEO record and paper disagree on some library-preparation and source
processing details: GEO describes an rRNA-depletion workflow and HISAT2, while
the paper describes KAPA mRNA preparation and STAR. The manifest therefore
uses the conservative label `bulk_RNA-seq` and records only project-generated
processing status. This discrepancy should be resolved in the final data
methods rather than silently choosing one description.

## Source evidence

- [GEO GSE315857](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315857)
- [Published study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13058273/)
- [Original accession mapping](../../endometriosis-splicing/samples.tsv)
- [Cleaver pilot record](../galaxy_data/october_20/cleaver_pilot/README.md)

## Checks performed

The manifest is valid when it has eight unique `sample_id` values, four
`control` rows, four `case` rows, unique GEO/SRA accessions, and no blank cells.
The Rust manifest parser accepts it because required `sample_id` and
`condition` columns are present; the extra provenance columns are allowed.
