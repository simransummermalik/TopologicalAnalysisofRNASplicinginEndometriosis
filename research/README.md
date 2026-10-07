# E-MTAB-15117 endometriosis RNA-seq dataset

This folder records a second public endometriosis RNA-seq dataset identified as
**E-MTAB-15117**. The original source sheet is preserved as
[research/dataorg.csv](dataorg.csv). Despite its filename, that file is
tab-separated MAGE-TAB sample metadata.

The source sheet describes 48 unique female human samples, all paired-end
PolyA-selected RNA-seq libraries with a nominal insert length of 300 bp. It
contains direct ENA FASTQ links and one ENA run accession per sample.

## Sample breakdown

| Condition | Samples | Tissue or site |
|---|---:|---|
| Control | 15 | Endometrium |
| Deep-infiltrating endometriosis | 16 | Body proper / lesion site |
| Ovarian endometriosis | 11 | Ovary |
| Peritoneal endometriosis | 6 | Peritoneal cavity |

The disease label for one source row (DE-14) contains a missing separator.
The reader preserves that original value in original_disease and maps it to
the normalized deep_infiltrating_endometriosis condition.

## Clean manifest

Run this from the repository root:

```sh
python3 research/read_dataorg.py
```

The script validates that each sample has exactly one R1 and one R2 FASTQ,
checks that paired rows agree on their metadata, normalizes condition labels,
and writes:

```text
data/metadata/E-MTAB-15117_samples.tsv
```

To validate without writing the output:

```sh
python3 research/read_dataorg.py --check-only
```

The generated manifest has one row per sample and includes the sample, tissue,
condition, ENA experiment and run IDs, and both FASTQ URLs. Keep the original
sheet unchanged so the normalized file can always be traced back to its source.

The first two-sample technical pilot is recorded in
[data/metadata/E-MTAB-15117_pilot.tsv](../data/metadata/E-MTAB-15117_pilot.tsv):
`B-01` is one normal-endometrium control and `DE-01` is one
deep-infiltrating-endometriosis lesion. The pilot is deliberately small enough
to test the workflow before downloading the rest of the study.

Before running the pilot, record the human genome build and annotation release
chosen by the research group, then build the matching STAR index and junction
gene map. Those files are not in this repository yet; the FASTQ download should
wait until that reference choice is recorded so reads are aligned consistently.

## Important comparison note

The control samples are normal endometrium. The disease samples are lesions
from body proper, ovary, or peritoneal cavity. This makes the dataset useful
for lesion-subtype and exploratory comparisons, but a direct disease-versus-
control result may also reflect tissue location. The research group should
choose the primary comparison before full-scale processing.

For the first technical pilot, use one control and one lesion sample. Do not
download all 48 samples or commit FASTQ/BAM files to Git until the pipeline has
passed that pilot.
