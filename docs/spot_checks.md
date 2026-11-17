# Independent spot checks: candidate-gene evidence

**For:** someone who did not build the evidence tables. The builder's own
audits (E8, G9) already passed; this check exists so that a second person
looks at the original sources.

**Time:** about 45–60 minutes for all ten checks.

**Status:** not yet done. Fill in each **Result** line, then sign off at the
bottom.

## How to check

1. Open the links for each row and compare what the source says with the
   **Claim** line.
2. On the **Result** line, write `pass` or `problem`, then a one-line note.
   For a `problem`, say what the source actually shows.
3. Do not edit the data files. If something is wrong, it is fixed through the
   review tables and recorded in [`gene_set_methods.md`](gene_set_methods.md).

The rows were chosen to cover the places most likely to hide a mistake:
genome-wide screens, genes not named in an abstract, weak citations,
multi-gene loci, and genes assigned from an author report. Each table also
has one row expected to pass easily, as a control.

## Environmental evidence (CTD): 5 rows

Each row comes from
[`environment_evidence.tsv`](../data/gene_sets/environment_evidence.tsv).
CTD's page for a gene's chemical interactions is linked as **CTD**.

### E-1. MMP2 · PCB 126 · cited paper on endometriosis

- **Row:** `CTD:136637:expression` (gene `ENSG00000087245.14`)
- **Claim:** 3,4,5,3',4'-pentachlorobiphenyl (PCB 126) increases MMP2 mRNA in
  human cells. At least one cited paper is about endometrium.
- **Sources:**
  - [CTD MMP2 interactions](https://ctdbase.org/detail.go?type=gene&acc=4313&view=ixn)
    (search the page for "pentachlorobiphenyl")
  - [PMID 27663891](https://pubmed.ncbi.nlm.nih.gov/27663891/) (endometriosis)
  - [PMID 29947894](https://pubmed.ncbi.nlm.nih.gov/29947894/)
  - [PMID 35618242](https://pubmed.ncbi.nlm.nih.gov/35618242/)
- **Look for:** in PMID 27663891's full text or figures, whether MMP2
  expression was measured after PCB 126 (CB126) in human endometrial cells,
  and in which direction.
- **Result:**

### E-2. ADM · TCDD · genome-wide screens

- **Row:** `CTD:2803579:expression` (gene `ENSG00000148926.12`)
- **Claim:** TCDD increases ADM mRNA in human cells. Two of the three cited
  papers are large screens; PMID 26238291 is the targeted paper that lets it
  pass the final rule.
- **Sources:**
  - [CTD ADM interactions](https://ctdbase.org/detail.go?type=gene&acc=133&view=ixn)
    (search for "Tetrachlorodibenzodioxin")
  - [PMID 26238291](https://pubmed.ncbi.nlm.nih.gov/26238291/)
  - [PMID 20106945](https://pubmed.ncbi.nlm.nih.gov/20106945/) and
    [PMID 21632981](https://pubmed.ncbi.nlm.nih.gov/21632981/) (screens)
- **Look for:** in PMID 26238291's supplementary gene lists, whether ADM is
  listed as changed by TCDD in human liver cells, and in which direction.
- **Result:**

### E-3. CTNNB1 · TCDD · gene not named in the abstract

- **Row:** `CTD:2815733:protein` (gene `ENSG00000168036.21`)
- **Claim:** TCDD increases β-catenin (CTNNB1) protein, from a single paper on
  hepatocyte dedifferentiation.
- **Sources:**
  - [CTD CTNNB1 interactions](https://ctdbase.org/detail.go?type=gene&acc=1499&view=ixn)
  - [PMID 22902829](https://pubmed.ncbi.nlm.nih.gov/22902829/)
- **Look for:** a figure or result showing β-catenin protein after TCDD in
  primary human hepatocytes, and in which direction.
- **Result:**

### E-4. CYP1A1 · TCDD · weak citation

- **Row:** `CTD:2816616:expression` (gene `ENSG00000140465.17`), context
  `dioxin_effect_modified`
- **Claim:** oxygen deficiency (hypoxia) reduces TCDD's induction of CYP1A1
  mRNA. The audit flagged PMID 20610613 as mainly a pathway model.
- **Sources:**
  - [CTD CYP1A1 interactions](https://ctdbase.org/detail.go?type=gene&acc=1543&view=ixn)
    (search for "Oxygen deficiency")
  - [PMID 19154183](https://pubmed.ncbi.nlm.nih.gov/19154183/)
  - [PMID 20610613](https://pubmed.ncbi.nlm.nih.gov/20610613/)
- **Look for:** whether either paper *measured* CYP1A1 mRNA under TCDD with
  and without hypoxia, rather than modeling it.
- **Result:**

### E-5. BRCA1 · TCDD · control

- **Row:** `CTD:2810658:activity` (gene `ENSG00000012048.28`), context
  `dioxin_effect_modified`
- **Claim:** TCDD (through AHR) increases BRCA1 transcriptional activity, and
  alpha-naphthoflavone blocks this.
- **Sources:**
  - [CTD BRCA1 interactions](https://ctdbase.org/detail.go?type=gene&acc=672&view=ixn)
  - [PMID 18259752](https://pubmed.ncbi.nlm.nih.gov/18259752/)
- **Look for:** the abstract statement about dioxin and alpha-naphthoflavone,
  and that the context label "another agent alters the dioxin's effect" fits.
- **Result:**

## Genetic-risk evidence: 5 genes

Each gene comes from
[`endometriosis_risk.tsv`](../data/gene_sets/endometriosis_risk.tsv).

**Supplementary Table 4** is the "Supplementary Table 4" sheet of
[Koller et al. 2026 Supplementary Tables (xlsx, 26 MB)](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-026-02582-2/MediaObjects/41588_2026_2582_MOESM4_ESM.xlsx).
Its positions are on GRCh37. Catalog positions are on GRCh38. To see where a
gene is, open its Ensembl page and check the location.

### R-1. ETAA1 · assigned from an author report, 217 kb away

- **Rows:** `GWASCAT:21308389` (gene `ENSG00000143971.11`)
- **Claim:** In GCST004549 (PMID 28537267), rs6546324 is associated with
  endometriosis at p = 3×10⁻⁸ and the authors report ETAA1. The Catalog itself
  maps two lncRNAs (LINC01812, LINC02831).
- **Sources:**
  - [Catalog variant rs6546324](https://www.ebi.ac.uk/gwas/variants/rs6546324)
  - [Catalog study GCST004549](https://www.ebi.ac.uk/gwas/studies/GCST004549)
  - [PMID 28537267](https://pubmed.ncbi.nlm.nih.gov/28537267/)
  - [ETAA1 in Ensembl](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000143971)
- **Look for:** the "Reported gene(s)" field for this association, and
  whether the paper names ETAA1 for this locus.
- **Result:**

### R-2. TFPI · split multi-gene locus

- **Row:** `KOLLER2026_ST4:rs8176498/2` (gene `ENSG00000003436.18`)
- **Claim:** Supplementary Table 4 row 16 lists lead SNP rs8176498
  (chr2:188353919, GRCh37), p = 8.67×10⁻¹³, nearest gene
  `AC007319.1:TFPI`. The SNP lies inside both genes, so both are kept.
  AC007319.1 is CALCRL-AS1.
- **Sources:**
  - Supplementary Table 4, sheet row 16
  - [TFPI in Ensembl](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000003436)
  - [rs8176498 in Ensembl](https://www.ensembl.org/Homo_sapiens/Variation/Explore?v=rs8176498)
- **Look for:** the row's values, and whether rs8176498 falls within TFPI.
- **Result:**

### R-3. FSHB · nearest gene to an intergenic variant

- **Row:** `KOLLER2026_ST4:rs11031005` (gene `ENSG00000131808.12`)
- **Claim:** Supplementary Table 4 row 58 lists rs11031005 (chr11:30226356,
  GRCh37), p = 1.07×10⁻⁶¹, nearest gene FSHB, about 26 kb from the gene.
- **Sources:**
  - Supplementary Table 4, sheet row 58
  - [FSHB in Ensembl](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000131808)
- **Look for:** the row's values.
- **Result:**

### R-4. CD109 · the two sources use different lead SNPs

- **Rows:** `GWASCAT:213218824` (rs2917893, GCST90693144, p = 1×10⁻¹¹) and
  `KOLLER2026_ST4:rs2917868` (row 33, chr6:74610073 GRCh37, p = 9.28×10⁻²¹)
- **Claim:** both rows map to CD109 (`ENSG00000156535.16`). rs2917893 is
  inside CD109; rs2917868 is about 72 kb away.
- **Sources:**
  - [Catalog variant rs2917893](https://www.ebi.ac.uk/gwas/variants/rs2917893)
  - Supplementary Table 4, sheet row 33
  - [CD109 in Ensembl](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000156535)
- **Look for:** the Catalog's mapped gene for rs2917893, and the nearest gene in
  row 33.
- **Result:**

### R-5. PAPPA2 · control

- **Rows:** `GWASCAT:218952499` and `KOLLER2026_ST4:rs10913186` (gene
  `ENSG00000116183.12`)
- **Claim:** rs10913186 is a Koller cross-ancestry locus (p ≈ 4.6×10⁻¹⁰),
  inside PAPPA2. The Catalog and Supplementary Table 4 row 6 agree.
- **Sources:**
  - [Catalog variant rs10913186](https://www.ebi.ac.uk/gwas/variants/rs10913186)
  - Supplementary Table 4, sheet row 6
- **Look for:** both sources give PAPPA2 and similar p-values.
- **Result:**

## Optional: reproductive-tissue tags

This check takes about 10 minutes. In
[`ctd_publication_tissue.tsv`](../data/gene_sets/ctd_publication_tissue.tsv),
read the 30 papers tagged `endometrium_uterus` or `other_female_reproductive`,
using each title and the `matched` column. Note any paper whose human
experiments clearly used other tissue.

- **Result:**

## Sign-off

- **Checked by:**
- **Date:**
- **Summary:** for example, "10 of 10 pass" or "9 pass, 1 problem (E-4: …)".
