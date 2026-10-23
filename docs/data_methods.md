# Data methods

This file records the reference and preprocessing choices that must stay the
same across every GSE315857 sample. Large reference files are stored locally in
the ignored `references/` directory; this document is the reviewable record
that belongs in Git.

## Human genome reference

Splice Girl uses the **GENCODE release 50 GRCh38.p14 primary assembly**. GENCODE
lists release 50 as the current GRCh38 release and provides a primary-assembly
FASTA whose sequence names match its primary-assembly annotation.

| Field | Recorded value |
|---|---|
| Provider | GENCODE / EMBL-EBI |
| GENCODE release | 50 |
| Assembly | GRCh38.p14 |
| Scope | Primary assembly: chromosomes and scaffolds; alternate loci, haplotypes, and assembly patches excluded |
| Filename | `GRCh38.primary_assembly.genome.fa.gz` |
| Local path | `references/gencode_v50/GRCh38.primary_assembly.genome.fa.gz` |
| Download date | October 22, 2026 |
| Compressed size | 845,635,028 bytes |
| Uncompressed size | 3,151,417,447 bytes |
| Official and observed MD5 | `da1a11258be075cfa7af718162c894e7` |
| Observed SHA-256 | `b760d18dbb651dd14dfc290083371b3ef3bff122d43a9cefb13ca4ecf38f05ca` |
| Gzip integrity | PASS |

Exact download URL:

```text
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh38.primary_assembly.genome.fa.gz
```

The compressed FASTA remains compressed to conserve local disk space. It is
ignored by Git and must not be pushed. Decompression and HISAT2 index creation
will be done only when enough space is available for the 3.15 GB FASTA and the
index at the same time.

## Matching annotation

Splice Girl uses the **GENCODE release 50 comprehensive primary-assembly GTF**.
It was downloaded from the same GENCODE release as the genome FASTA and checked
before being recorded in the sample manifest.

| Field | Recorded value |
|---|---|
| Annotation | GENCODE v50 comprehensive primary-assembly GTF |
| Filename | `gencode.v50.primary_assembly.annotation.gtf.gz` |
| Local path | `references/gencode_v50/gencode.v50.primary_assembly.annotation.gtf.gz` |
| Download date | October 22, 2026 |
| Compressed size | 124,650,284 bytes |
| Official and observed MD5 | `289b91e5e95e8b0450d223246f10a12e` |
| Observed SHA-256 | `89bbad69a8c89fee5fadec0a2f14a098752d27a4ef7f24e4de9bac681e1b18f4` |
| Gzip integrity | PASS |
| Annotation rows | 11,254,876 total; 5,087,789 exon rows |

Exact download URL:

```text
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.primary_assembly.annotation.gtf.gz
```

Using the primary-assembly FASTA and primary-assembly GTF from the same GENCODE
release avoids chromosome-name and coordinate mismatches. The FASTA contains
194 sequence names and the GTF uses 70 of them. Every GTF sequence name exists
in the FASTA; there are no annotation sequence names missing from the genome.
The 124 FASTA sequences without GTF rows are unannotated reference sequences,
not naming conflicts.

Like the FASTA, the compressed GTF is ignored by Git. The checksums, source,
and compatibility result above are the reproducible record kept in the
repository.

## Frozen junction and graph conventions

These rules apply to every real GSE315857 sample. They keep the biological
meaning of a node or edge from changing between samples.

### Coordinates

- GENCODE GTF exon coordinates are **1-based and closed**: both `start` and
  `end` are included in the exon.
- A raw splice junction is represented as
  `(sequence, intron_start, intron_end, strand)`, where the two intron
  coordinates are also 1-based and closed.
- For a SAM/BAM alignment, `POS` is 1-based. A CIGAR `N` operation consumes
  the skipped reference interval; its first and last skipped bases become
  `intron_start` and `intron_end`.
- STAR `SJ.out.tab` uses the same 1-based, closed intron interval. The existing
  STAR conversion fixture remains useful for checking the convention even
  though HISAT2 is the selected aligner for the real pilot.
- Sequence names are copied exactly from the GENCODE files, including the
  `chr` prefix. No automatic renaming between `1` and `chr1` is allowed.

The coordinate rules follow the
[GENCODE GTF format](https://www.gencodegenes.org/pages/data_format.html), the
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf), and the
[STAR junction specification](https://github.com/alexdobin/STAR/blob/master/extras/doc-latex/STARmanual.tex).

### Graph nodes

One graph node is one **unique annotated exon interval within a gene**, not a
single splice-site coordinate. Exact duplicate exon intervals appearing in
multiple transcripts of the same gene collapse to one node. This convention
connects consecutive junctions through their shared exon and makes an
annotated exon-skipping pattern form the intended triangle.

The node identifier is:

```text
sequence:exon_start-exon_end:strand
```

For example, `chr1:100-199:+` is a plus-strand exon covering bases 100 through
199. `exon_number` is not used because it can differ between transcripts.
`gene_id` is stored separately and preserves the versioned GENCODE identifier,
such as `ENSG00000123456.7`.

### Directed edges and strand

An edge is directed in biological transcript order, from the upstream exon to
the downstream exon. Genomic coordinate order alone never determines the
direction.

| Strand | Source (`from`) | Destination (`to`) | Direction on the reference |
|---|---|---|---|
| `+` | Exon whose `end = intron_start - 1` | Exon whose `start = intron_end + 1` | Lower to higher coordinates |
| `-` | Exon whose `start = intron_end + 1` | Exon whose `end = intron_start - 1` | Higher to lower coordinates |

For plus-strand exons `chr1:100-199:+` and `chr1:300-399:+`, the intron is
`chr1:200-299:+` and the edge is `chr1:100-199:+` to
`chr1:300-399:+`. For an otherwise identical minus-strand pair, the edge runs
from `chr1:300-399:-` to `chr1:100-199:-`.

Junctions with unknown strand are not guessed into a direction. They are
excluded from the primary graph and counted in conversion QC.

### Annotation mapping and exclusions

The annotation map is built from exon pairs that are adjacent within at least
one GENCODE transcript. A junction is kept only when its sequence, intron
coordinates, and strand identify exactly one `(gene_id, from_node, to_node)`
mapping after duplicate mappings are collapsed.

- No annotation match: record the junction as unannotated and exclude it.
- More than one gene or exon-pair match: record it as ambiguous and exclude it.
- Exactly one match: keep it and aggregate supporting reads for that directed
  edge within the sample.
- Primary counts use uniquely mapped junction-spanning reads. Multimapping
  support remains excluded unless a later reviewed rule explicitly changes
  that setting for every sample.

These exclusions favor a smaller graph with clear biological node identities
over silently assigning a junction to the wrong exon or gene. The conversion
QC must report the number of input, unannotated, ambiguous, filtered, and
written junctions.

### Interpretation boundary

Edge direction records transcript order. The Hodge cycle component measures
graph cycle-space structure after the edges are assembled; it does not require
a literal directed return path. Therefore a nonzero cycle component can
describe alternative paths such as exon skipping and does not by itself prove
circular RNA or a back-splice junction.

## Remaining D5 work

The genome, annotation, coordinate system, node identity, and strand direction
are now frozen. D5 remains in progress only until the HISAT2 index and known
splice-site hints are built and verified from this exact FASTA/GTF pair.
