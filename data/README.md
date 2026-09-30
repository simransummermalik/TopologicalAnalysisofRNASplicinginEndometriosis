# Data directory

The first public dataset selected for review is **GEO GSE179640**:

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179640>

It contains human bulk RNA-seq samples from control endometrium, eutopic
endometrium from people with endometriosis, and ectopic lesions. The primary
comparison marked in the manifest is:

- 5 control endometrium samples;
- 7 eutopic endometrium samples from endometriosis cases.

The 12 ectopic lesion samples are listed for context and future secondary
comparisons. They are not part of the first control-versus-eutopic comparison.

The GEO record identifies the linked sequencing study as BioProject
`PRJNA744463` / SRA study `SRP327335`.

## Files here

- `geo/GSE179640_series_matrix.txt.gz` is the small GEO metadata snapshot.
- `metadata/GSE179640_bulk_samples.tsv` is the cleaned 24-sample manifest.

The GEO record also links to per-sample featureCounts files and a much larger
raw archive. The featureCounts files are gene-level counts, not the junction
table required by Splice Girl. The raw archive is about 0.8 GB, so it is not
copied into this repository. Junction counts still need to be derived from
appropriate alignments before the data can be passed to the Rust tool.
One possible later extraction route is the public recount3 exon-exon junction
resource for SRA projects; that matrix has not been downloaded here yet.

The manifest preserves GEO, BioSample, SRA, and processed-file identifiers so
the source data can be retrieved later. Do not treat the `featurecounts_url`
column as a Splice Girl input path.

The downloaded metadata snapshot has SHA-256 checksum
`50caaef6ffcfa5d6e11597202945793f98ea678d159c02d0b2973077b66564c4`.
