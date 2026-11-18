# Real-data technical pilot plan — November 17, 2026

This checks one real control sample through alignment, junction extraction,
Rust analysis and an independent numerical calculation. It is not the final
cohort analysis and will not calculate a disease effect or p-value.

- Sample: GSE315857 control `SRR36710387`, using the existing Cleaver-cleaned
  R1/R2 files. Their SHA-256 values must match the October 22 pilot record.
- Input scope: the **first 1,000,000 synchronized read pairs** (2,000,000 reads)
  of 49,476,792 cleaned pairs. This deterministic prefix may have ordering bias;
  it is not a random sample or a new biological replicate.
- Reference: the full verified GENCODE v50 GRCh38.p14 primary assembly and
  matching GTF. Do not restrict the reference to the selected gene; uniqueness
  must be assessed against the complete reference.
- Index: HISAT2 2.2.3, two threads, suffix-sort bucket divisor 16 to limit build
  memory. The temporary expanded FASTA copy can be removed after successful
  indexing; the original compressed reference is preserved.
- Alignment: the existing HISAT2 adapter with known splice-site hints and
  `--dta-cufflinks`. Use an explicitly **unstranded exploratory alignment**;
  the final library-strandedness decision remains subject to protocol review.
- Base junction table: primary uniquely mapped read alignments (`NH:i:1`),
  known transcript strand (`XS`), at least 8 aligned bases on each junction
  side, minimum junction support 1, exact unambiguous annotated exon-pair map.
  Duplicate flags remain retained; counts are reads, not deduplicated fragments.
- Pick one gene from the current `data/gene_sets/candidate_genes.tsv` using
  this rule: retain edges supported by at least 3 reads; require at least 3
  retained edges, at most 12 observed exon nodes, and total retained support
  at least 100. Select the gene with greatest retained support, then lexical
  `gene_id` for ties. Do not use cycle fraction, p-values or case/control
  differences for this choice. If none qualify, record that failure.
- Independently recount the chosen gene's eligible junctions from its BAM
  region and compare every count to the converter output.
- Compare the Rust gradient/cycle values to NumPy's SVD-based least-squares
  projection, and check reconstruction, divergence and orthogonality.
- Recalculate this same gene at supports 1, 3, 5 and 10, keeping all attempted
  settings visible. This is technical sensitivity, not an approved study grid.

The small selected-gene table, verification results and figure may be retained
here for review. FASTQs, the full junction table, reference/index, SAM and BAM
remain in ignored local folders. No full-sample or case/control completion
will be inferred from this technical pilot.
