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

## Remaining D5 decisions

D5 remains in progress until all of the following are recorded and checked:

1. Confirm the final exon-versus-splice-site node convention.
2. Record the conversion from 1-based GTF coordinates to junction endpoints.
3. Record strand handling.
4. Build and verify the HISAT2 index and splice-site hints from this exact
   FASTA/GTF pair.
