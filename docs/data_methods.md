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

The matching annotation has been selected but has not yet been downloaded or
verified:

| Field | Planned value |
|---|---|
| Annotation | GENCODE v50 comprehensive primary-assembly GTF |
| Filename | `gencode.v50.primary_assembly.annotation.gtf.gz` |
| Official MD5 | `289b91e5e95e8b0450d223246f10a12e` |
| Status | Pending download and integrity check |

Exact planned URL:

```text
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.primary_assembly.annotation.gtf.gz
```

Using the primary-assembly FASTA and primary-assembly GTF from the same GENCODE
release avoids chromosome-name and coordinate mismatches.

## Remaining D5 decisions

D5 remains in progress until all of the following are recorded and checked:

1. Download and verify the matching GTF.
2. Confirm the final exon-versus-splice-site node convention.
3. Record the conversion from 1-based GTF coordinates to junction endpoints.
4. Record strand handling.
5. Build and verify the HISAT2 index and splice-site hints from this exact
   FASTA/GTF pair.
