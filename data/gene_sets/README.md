# Candidate-gene sets

This folder holds the Priority 2 evidence tables and the gene sets built from
them. The rules for what counts as evidence are in
[`docs/gene_set_methods.md`](../../docs/gene_set_methods.md); people apply
those rules. The helper [`gene_sets/gene_sets.py`](../../gene_sets/gene_sets.py)
only maps identifiers, builds membership files, checks schemas, and computes
the overlap. It never decides whether a gene counts.

## Files

| File | Made by | One row means |
|---|---|---|
| `risk_evidence.tsv` | `extract_gwas.py` (G4–G5), then `map`, then people (G7) | One evidence assertion from one source |
| `dioxin_like_chemicals.tsv` | People (E1), frozen | One dioxin-like chemical with its exact CTD ID |
| `environment_evidence.tsv` | `extract_ctd.py` (E3–E4), `map`, then `review_ctd_context.py apply` (E6) | One measured response to one chemical in one context |
| `environment_mapping_review.tsv` | People (E5 step 6), applied by `review` | One reviewed label: mapped by hand or excluded, with a reason |
| `risk_mapping_review.tsv` | People (G7), applied by `review` | One reviewed risk row (keyed by `evidence_id`): mapped, split across overlapping genes, or excluded, with a reason |
| `risk_audit_sample.tsv` | `audit_risk.py` (G9), then a person | One sampled risk gene with its check against the sources |
| `environment_audit_sample.tsv` | `audit_environment.py` (E8), then a person | One sampled evidence row with its check against the cited papers |
| `ctd_publication_tissue.tsv` | `review_ctd_context.py tag-publications` (E6) | One cited PubMed paper with its tissue tag and the terms that decided it |
| `endometriosis_risk.tsv` | `dedupe --kind risk` | One gene |
| `environment_response.tsv` | `dedupe --kind environment` | One gene × chemical × tissue × context |
| `overlap.tsv` | `overlap` | One gene in both sets |
| `candidate_genes.tsv` | `overlap` | One gene in either set, with category flags |
| `environment_genes.tsv` | `overlap` | One environmental gene per row (565), so the set can be a `stats` gene family |
| `*_mapping.tsv` | `map` | One original label and its candidate IDs, for G7 review |

`risk_evidence.tsv` holds the GWAS Catalog and Koller Supplementary Table 4
extraction from November 17, 2026 after the G7 review (1,119 rows; 228 mapped
to 103 genes, the rest excluded with reasons). `endometriosis_risk.tsv` is the
audited G8–G10 membership. `overlap.tsv` (7 genes) and `candidate_genes.tsv`
(661 genes) are the E9 outputs; `candidate_genes.tsv` is the one-row-per-gene
file to pass to `splice-girl stats --genes`. `environment_evidence.tsv`
holds the CTD extraction from November 17, 2026 (7,958 rows) after the E6
context review: 1,028 co-treatment rows are `excluded` with a reason, and
every row has a paper-level tissue tag. The 12 labels `map` could not
resolve were reviewed by hand (1 mapped, 11 excluded). `environment_response.tsv`
is the E7 membership under the final publication rule (decision 2026-11-17). See
`docs/gene_set_methods.md`. Do not add example
rows to either file.

## Shared rules

- Tab-separated, UTF-8, one header row. Columns must appear in the order shown
  by the template.
- No blank cells. Write `NA` for a value that is missing or does not apply.
- `gene_id` is the versioned GENCODE v50 ID, for example `ENSG00000123456.7`.
  It is the only join key. `gene_symbol` is for reading and is never joined on.
- `original_id` is copied from the source exactly as written and never edited.
- Lists inside one cell use `|`, for example `PMID:1|PMID:2`.

## Evidence columns

Both evidence tables share `evidence_id` (unique), `source`, `source_version`,
`source_record` (table and row, accession, or similar), `original_id`,
`gene_id`, `gene_symbol`, `mapping_status`, `mapping_method`, and `notes`.

**Risk only:** `variant_id`, `p_value` (0–1), `evidence_type`,
`evidence_level`, `population`.

**Environment only:** `publications` (required, for the publication rule),
`exposure` (CTD chemical name), `exposure_id` (exact CTD ID such as
`MESH:D013749`), `dose_duration`, `species`, `tissue`, `assay`,
`response_type`, `evidence_type`, `context`.

### Allowed values

| Column | Values |
|---|---|
| `mapping_status` | `pending` (not mapped yet), `mapped`, `ambiguous`, `unresolved`, `excluded` |
| `evidence_type` (risk) | `gwas_association` |
| `evidence_level` (risk) | `variant`, `locus`, `gene` |
| `evidence_type` (environment) | `measured_response` |
| `response_type` (environment) | `expression`, `protein`, `activity` |

`gene_id` and `gene_symbol` must be `NA` unless `mapping_status` is `mapped`.

## Workflow

Run from the repository root. `references/` is the ignored local folder named
in [`docs/data_methods.md`](../../docs/data_methods.md).

```sh
# 0a. Risk: draft rows from the GWAS Catalog release and Koller et al. 2026
#     Supplementary Table 4 recorded in docs/gene_set_methods.md. --fetch
#     downloads Catalog responses into the cache folder. This overwrites
#     risk_evidence.tsv, so run it only before any review edits.
python3 gene_sets/extract_gwas.py \
  --catalog-dir /path/to/gwas_cache --fetch \
  --koller-xlsx /path/to/41588_2026_2582_MOESM4_ESM.xlsx \
  --output data/gene_sets/risk_evidence.tsv

# 0b. Environment only: draft rows from the CTD release recorded in
#    docs/gene_set_methods.md (download it outside the repository and check
#    its SHA-256 first). This overwrites environment_evidence.tsv, so run it
#    only before any review edits.
python3 gene_sets/extract_ctd.py \
  --ctd /path/to/CTD_chem_gene_ixns.tsv.gz \
  --chemicals data/gene_sets/dioxin_like_chemicals.tsv \
  --output data/gene_sets/environment_evidence.tsv

# 1. People enter rows with mapping_status = pending. Then map identifiers.
#    For environment, follow map with the E6 review step:
#      python3 gene_sets/review_ctd_context.py tag-publications \
#        --evidence data/gene_sets/environment_evidence.tsv \
#        --pubmed-dir /path/to/pubmed_xml --fetch \
#        --output data/gene_sets/ctd_publication_tissue.tsv
#      python3 gene_sets/review_ctd_context.py apply \
#        --evidence data/gene_sets/environment_evidence.tsv \
#        --publications data/gene_sets/ctd_publication_tissue.tsv \
#        --output data/gene_sets/environment_evidence.tsv
python3 gene_sets/gene_sets.py map --kind risk \
  --evidence data/gene_sets/risk_evidence.tsv \
  --gtf references/gencode_v50/gencode.v50.primary_assembly.annotation.gtf.gz \
  --hgnc references/hgnc/hgnc_complete_set.txt \
  --output data/gene_sets/risk_evidence.tsv \
  --mapping-table data/gene_sets/risk_mapping.tsv

# 2. People review every ambiguous and unresolved row (G7 / E5 step 6) and
#    record each decision in a review table, then apply it. The first column
#    is original_id (one decision per label) or evidence_id (one per row);
#    then decision (mapped|excluded), gene_id, and reason. A mapped gene_id
#    may list several IDs with "|" to split a row across genes the source
#    names equally. A gene GENCODE leaves unnamed gets gene_symbol NA.
#    Only rows still ambiguous or unresolved change, so reruns are safe:
python3 gene_sets/gene_sets.py review --kind environment \
  --evidence data/gene_sets/environment_evidence.tsv \
  --review data/gene_sets/environment_mapping_review.tsv \
  --gtf references/gencode_v50/gencode.v50.primary_assembly.annotation.gtf.gz \
  --output data/gene_sets/environment_evidence.tsv

# 3. Check, then build membership:
python3 gene_sets/gene_sets.py check --kind risk_evidence \
  --table data/gene_sets/risk_evidence.tsv \
  --gtf references/gencode_v50/gencode.v50.primary_assembly.annotation.gtf.gz
python3 gene_sets/gene_sets.py dedupe --kind risk \
  --evidence data/gene_sets/risk_evidence.tsv \
  --output data/gene_sets/endometriosis_risk.tsv

# Environment: same steps with --kind environment. The final publication
# rule (decision 2026-11-17) is applied here: at least 2 publications per
# gene-chemical pair, at least 1 of them not a large screen (a publication
# cited by 100+ mapped rows):
python3 gene_sets/gene_sets.py dedupe --kind environment --min-publications 2 \
  --min-targeted 1 --screen-rows 100 \
  --evidence data/gene_sets/environment_evidence.tsv \
  --output data/gene_sets/environment_response.tsv

# Risk audit (G9): re-checks both risk tables against the cached Catalog
# responses, the supplement, the review table, GENCODE v50 and v19, and writes
# a fixed-seed sample of genes for a person to read. Exits 1 if a check fails.
python3 gene_sets/audit_risk.py \
  --evidence data/gene_sets/risk_evidence.tsv \
  --membership data/gene_sets/endometriosis_risk.tsv \
  --review data/gene_sets/risk_mapping_review.tsv \
  --catalog-dir /path/to/gwas_cache \
  --koller-xlsx /path/to/41588_2026_2582_MOESM4_ESM.xlsx \
  --gtf references/gencode_v50/gencode.v50.primary_assembly.annotation.gtf.gz \
  --gtf-v19 /path/to/gencode.v19.annotation.gtf.gz \
  --sample-output data/gene_sets/risk_audit_sample.tsv

# Environment audit (E8): re-checks both tables against the raw CTD file and
# references, rebuilds the membership, and writes a fixed-seed sample for a
# person to check against the cited papers. Exits 1 if any check fails.
python3 gene_sets/audit_environment.py \
  --evidence data/gene_sets/environment_evidence.tsv \
  --membership data/gene_sets/environment_response.tsv \
  --chemicals data/gene_sets/dioxin_like_chemicals.tsv \
  --ctd /path/to/CTD_chem_gene_ixns.tsv.gz \
  --gtf references/gencode_v50/gencode.v50.primary_assembly.annotation.gtf.gz \
  --hgnc references/hgnc/hgnc_complete_set.txt \
  --min-publications 2 --min-targeted 1 --screen-rows 100 \
  --sample-output data/gene_sets/environment_audit_sample.tsv

# 4. Overlap and category counts (E9–E10):
python3 gene_sets/gene_sets.py overlap \
  --risk data/gene_sets/endometriosis_risk.tsv \
  --environment data/gene_sets/environment_response.tsv \
  --overlap data/gene_sets/overlap.tsv \
  --candidates data/gene_sets/candidate_genes.tsv \
  --environment-genes data/gene_sets/environment_genes.tsv
```

`map`, `dedupe`, and `overlap` print counts as JSON. Copy them into the
extraction notes and audit records.

### How `map` decides

It tries each step in order and stops at the first match:

1. The label is an Ensembl gene ID: look it up in the GTF.
2. The label matches exactly one GTF `gene_name`. If GENCODE gives that name
   to several gene IDs, HGNC's Ensembl link for the approved symbol is used
   only when it picks out exactly one of them (`gencode_symbol_hgnc_link`).
3. The label is an HGNC approved, then previous, then alias symbol for
   exactly one gene. For example `BIM` → `BCL2L11`.

A label naming several genes (`GREB1, KDR`), a name shared by several genes,
or an alias used by several genes is marked `ambiguous`, and the candidates
are listed in `notes`, capped at 20 for multi-copy names such as `Y_RNA`. The
mapping table's `n_candidates` column gives the full count. Anything else is
`unresolved`. The script never picks one gene for you. A PAR_Y copy is ignored
when the primary copy exists. Genes that sit only on alternate haplotypes
(for example `HLA-DRB4` and `KIR2DL5A`) are not in the primary-assembly GTF,
so they come out `unresolved` as `hgnc_approved_not_in_release`.
