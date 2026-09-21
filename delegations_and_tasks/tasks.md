# Splice Girl — Task Guide

> **Shared task list:** this file contains general task.

This file expands the task IDs in `ROADMAP.md` into small tasks that anyone
can pick up one at a time. It is a shared checklist for the team.
Split the tasks in whatever way is convenient; this file only describes the
work and its expected outputs.

## How to use this file

1. Pick one task number.
2. Follow the small numbered steps under that task.
3. Save the file, table, code, or result the task asks for.
4. Compare the result with the task's **Done when** sentence.
5. Mark the task complete in the shared checklist.

If something is missing, write down what is missing and move to another task.
You do not need to understand the entire project before starting one small task.

The letters are simple categories: **G** genetics, **E** environmental evidence,
**D** data, **M** math, **R** Rust, **C** sample comparison, **S** statistics,
**I** interpretation, **U** usability, **V** presentation, **F** final documents,
and **X** optional extensions. For example, `M12` means math task 12 and `R5`
means Rust task 5.

## How to use this list

Split task IDs according to available time and interest. A task can be completed
alone or worked on together. A simple shared table is enough:

| Task ID | Status | Notes or link to output |
|---|---|---|
| P1.1 | Not started | |
| G1 | Not started | |
| M1 | Not started | |
| R1 | Not started | |

Add rows as needed. When a task is finished, link its output and change the
status to **Done**. When a task is blocked, write the missing information in
the Notes column and choose another task until the blocker is resolved.

Every task has a **Source** line. Open those links first to see the project
decision or explanation behind the task. If the task requires a new paper,
database, or dataset, record its exact citation in the output instead of
leaving the source in a chat message.

## What the words mean

These definitions are enough to begin the technical tasks.

- A **gene** is a named region of biological information. The project will
  analyze selected genes rather than every gene immediately.
- A **gene ID** is a stable identifier for a gene. A gene symbol is the short,
  human-readable name. Store both when possible.
- A **junction** is an observed connection between two exons or splice sites.
  In the input table, `from` and `to` name its endpoints.
- A **sample** is one biological specimen or one separately measured library.
  Never treat all patients as one sample just because they have the same label.
- **Case** means an endometriosis sample in the selected dataset. **Control**
  means an eligible comparison sample. Document the labels in the data
  manifest.
- A **TSV** is a text table whose columns are separated by tab characters.
- A **directed graph** has an edge with a start and an end. For example,
  `E1 -> E2` is different from `E2 -> E1`.
- The **incidence matrix** `B1` records which node each directed edge leaves
  and enters. The source receives `-1`; the destination receives `+1`.
- A **gradient component** is the portion of an edge signal explainable by
  assigning values to nodes.
- A **cycle-space component** is the mathematically circulating portion left
  after the gradient projection. It is not automatically proof of circular RNA.
- A **fixture** is a tiny known input used to test code.
- A **tolerance** is the small numerical error accepted because computers store
  decimal values approximately.
- A **baseline** is a simpler measurement used to check whether the Hodge-based
  summary adds useful information.
- A **p-value** measures how unusual an observed difference is under a stated
  randomization procedure. An **FDR/q-value** adjusts for testing many genes.

## Dependencies at a glance

```text
Phase 1: project setup, decisions, and file formats
       |
       +--> G tasks: risk genes -----------+
       |                                    |
       +--> E tasks: environmental genes ---+--> candidate gene tables
       |
       +--> D tasks: metadata and junction tables
       |
       +--> M tasks: synthetic math fixtures
                         |
                         +--> R1-R9: Rust synthetic MVP
                                          |
candidate genes + real data + MVP --------+--> R10-R14: integration and outputs
                                                         |
                                                         +--> S tasks: statistics
                                                         +--> V tasks: figures and demo
```

Genetics, environmental evidence, dataset suitability, synthetic math, parser
design, and presentation planning can begin at the same time. Real biological
interpretation must wait for data QC, validated per-sample output, and the
statistics review.

## Phase 1 — Organize the project

These tasks create the shared structure before anyone depends on a private
spreadsheet or an undocumented convention.

### P1.1 — Make a simple task list
**Source:** [ROADMAP.md](../ROADMAP.md), Phase 1 — project setup.

**Output:** a shared task checklist.

1. Copy the task IDs from this file into a shared table.
2. Add a status column with **Not started**, **In progress**, **Blocked**, and
   **Done**.
3. Add a Notes or link column for the file, table, code, or result produced.
4. Start with a few small tasks instead of trying to do everything at once.
5. Leave tasks that have not started as **Not started**.

**Done when:** every task has a visible row and the team can see what still
needs to be done.

### P1.2 — Record dates and checkpoints
**Source:** [ROADMAP.md](../ROADMAP.md), Phase 1 — project setup.

**Output:** project dates in the task board.

1. Write down the actual course start, midpoint, rehearsal, and final dates.
2. Put the midpoint checklist from `ROADMAP.md` into the board.
3. Put the final checklist into the board.
4. Add a target date for the first completed tasks.
5. Add a weekly review time.

**Done when:** the team can answer what should exist by the midpoint and what
must be ready before the final presentation.

### P1.3 — Recover the original Overleaf source
**Source:** [ROADMAP.md](../ROADMAP.md), Phase 1 — project setup.

**Output:** `docs/source_map.md`.

1. Search the repository for `.tex` files.
2. Search the project notes for the location of the main Overleaf source if it is absent.
3. Record the main `.tex` filename and its folder.
4. Record the document class, margins, packages, and section input pattern.
5. Check whether `roadmap.tex` should be included with `\input{roadmap}` or a
   `sections/` path.
6. Do not reconstruct the proposal silently if the source cannot be found.
7. Record the missing-source limitation in `docs/source_map.md`.

**Done when:** a new reader knows which file controls the Overleaf document.

### P1.4 — Freeze the first data contract
**Source:** [ROADMAP.md](../ROADMAP.md), Phase 1 — project setup.

**Output:** `docs/data_contract.md`.

1. Create a list of required input columns.
2. Start with `sample_id`, `gene_id`, `from`, `to`, and `count`.
3. Define what each column means in one sentence.
4. State which fields may be empty. Prefer none for the first prototype.
5. State that `count` must be finite and nonnegative.
6. State whether duplicate directed edges are rejected or combined.
7. State how a zero-support group is reported.
8. State the reference or coordinate convention once it is known.
9. Give the contract a version number.
10. Check the contract against two example tables before real-data use.

**Done when:** the same valid example table can be created without asking a
hidden question.

### P1.5 — Create the decisions log
**Source:** [ROADMAP.md](../ROADMAP.md), Phase 1 — project setup.

**Output:** `docs/decisions.md`.

1. Make a table with columns: date, decision, participants, reason, and
   affected files.
2. Add the current CLI name: `splice-girl`.
3. Add decisions about node representation, identifiers, normalization, and
   missing values as they are approved.
4. Never overwrite an old decision; add a new dated entry if it changes.
5. Link each decision to the task ID it unblocks.

**Done when:** a new reader can understand why a convention exists without
searching chat history.

### P1.6 — Create the task board
**Source:** [ROADMAP.md](../ROADMAP.md), Phase 1 — project setup.

**Output:** a board or table with all task IDs.

1. Enter IDs P1.1–P1.6, G1–G10, E1–E11, D1–D12, M1–M16, R1–R17,
   S1–S12, and V1–V11.
2. Add status, dependency, due date, and output columns.
3. Copy the task's completion condition into a short note.
4. Move a task to **Blocked** only with a named blocker.
5. Review the board at each team meeting.

**Done when:** every item in this document has a visible place on the board.

### P1.7 — Set up checking and file naming rules
**Source:** [ROADMAP.md](../ROADMAP.md), Phase 1 — project setup.

**Output:** a short section in `docs/decisions.md`.

1. Choose lowercase filenames with underscores for tables.
2. Choose where raw, processed, results, figures, and documentation belong.
3. Check gene lists, data manifests, math fixtures, and result tables before
   using them in later tasks.
4. Require a run manifest for real analysis.
5. Decide whether large raw files stay outside Git.

**Done when:** new files have predictable locations and no one needs to guess
whether a table is raw or processed.

## Phase 2A — Endometriosis genetics tasks

These tasks create the genetic-susceptibility gene set already described in the
proposal. Code should not be used to decide which genes count as biological
evidence.

### G1 — Write inclusion and exclusion rules
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2A — genetics.

**Output:** the criteria section of `docs/gene_set_methods.md`.

1. Write what “genetic susceptibility” means for this project.
2. Decide which evidence types qualify.
3. Decide what expression-only evidence does not qualify as inherited risk.
4. Decide how variant, locus, and gene-level evidence will be recorded.
5. Write how ambiguous evidence will be marked.
6. Write what will be excluded and why.
7. Ask someone unfamiliar with the rules to read them without explanation.
8. Revise unclear words.
9. Save the dated version of the rules.

**Done when:** a reader can classify a new evidence row using only the
written rules.

### G2 — List the evidence sources
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2A — genetics.

**Output:** a source list in `docs/gene_set_methods.md`.

1. List every proposed literature or database source before extracting genes.
2. Record the source name and exact access location.
3. Record the version, release, or publication date.
4. Record the date the source was accessed.
5. Keep a separate extraction block for each source.
6. Do not describe a source as used until its evidence rows are saved.

**Done when:** every future risk-gene row can point to a named source.

### G3 — Make the evidence table
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2A — genetics.

**Output:** an empty table template such as
`data/gene_sets/risk_evidence.tsv`.

1. Create the header row.
2. Include stable gene ID, gene symbol, original identifier, evidence type,
   source, source version, mapping method, and notes.
3. Add an example row only if it is clearly labeled as an example.
4. Remove example rows before real extraction.
5. Document whether one row means one evidence assertion.
6. Check that the source can be traced from each row.

**Done when:** the template can preserve multiple sources supporting one gene.

### G4 — Extract evidence from the first source
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2A — genetics.

**Output:** evidence rows with provenance.

1. Open the approved source.
2. Copy the source identifier into the evidence table.
3. Record the original gene or locus label exactly as shown.
4. Record the evidence type using the approved vocabulary.
5. Add the tissue, study, or population context when available.
6. Do not normalize the identifier by guessing.
7. Mark a missing mapping for review.
8. Save the table and record how many rows were added.

**Done when:** every extracted row can be traced back to its source.

### G5 — Extract remaining sources
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2A — genetics.

**Output:** the completed raw risk evidence table.

1. Give each source a separate extraction block or source label.
2. Use the same column definitions as G4.
3. Do not silently change inclusion rules for a new source.
4. Record a reason when a source produces no usable genes.
5. Count rows by source.
6. Compare the source counts with the extraction notes.

**Done when:** all agreed sources have either evidence rows or a documented
reason for having none.

### G6 — Normalize identifiers
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2A — genetics.

**Output:** mapped stable IDs with a mapping record.

1. Choose the stable gene-ID namespace and reference release recorded in the
   data contract.
2. Create a mapping table from each original identifier to the stable ID.
3. Keep the original identifier; never delete it after mapping.
4. Keep the readable gene symbol in its own column.
5. Mark unmapped identifiers as unresolved.
6. Mark one-to-many mappings for human review.
7. Record the mapping method and release.
8. Inspect a sample of mappings a second time.

**Done when:** every included gene has a stable ID or a documented reason it is
still unresolved.

### G7 — Resolve ambiguous mappings
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2A — genetics.

**Output:** an ambiguity decision table.

1. Filter the mapping table for missing, one-to-many, and conflicting rows.
2. For each row, record the possible mappings.
3. Check the source context before selecting one.
4. If no safe selection exists, exclude the row with a reason.
5. Do not replace an ambiguous ID with a convenient symbol.
6. Count resolved and excluded rows.

**Done when:** no ambiguous row disappears without a recorded decision.

### G8 — Deduplicate membership
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2A — genetics.

**Output:** `data/gene_sets/endometriosis_risk.tsv` plus retained evidence.

1. Group rows by stable gene ID.
2. Create one membership row per gene.
3. Keep all evidence rows in the evidence table.
4. Combine source names without losing source-specific detail.
5. Keep a readable symbol, but use the stable ID for joins.
6. Count genes before and after deduplication.
7. Check for duplicate stable IDs in the final membership file.

**Done when:** membership is unique by stable ID and provenance remains available.

### G9 — Perform an independent audit
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2A — genetics.

**Output:** a review note.

1. Check the criteria, evidence table, mapping table, and final set together.
2. Inspect random included genes.
3. Inspect every unresolved mapping.
4. Ask whether any expression-only evidence was accidentally included.
5. Record every correction.
6. Re-export the table after corrections.
7. Record the check date and any corrections.

**Done when:** the criteria, provenance, mapping, and uniqueness checks are
recorded.

### G10 — Export and schema-check the risk set
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2A — genetics.

**Output:** `data/gene_sets/endometriosis_risk.tsv`.

1. Check that the header has the required fields.
2. Check that every row has a stable ID and evidence type.
3. Check that stable IDs are unique in the membership file.
4. Check that symbols are separate from IDs.
5. Check that source names are not blank.
6. Check that the file can be joined to the data tables and Rust outputs.
7. Save a count and checksum if the team uses checksums.

**Done when:** the file passes the schema check and its provenance audit.

## Phase 2B — Environmental-response tasks

These tasks prioritize genes. Unless the dataset contains patient-level
exposure information, they do not prove that an exposure caused an observed
splicing pattern.

### E1 — Choose the exposure scope
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2B — environmental evidence.

**Output:** the exposure definition in `docs/gene_set_methods.md`.

1. List the exposure or exposure class being considered.
2. State whether the scope is one chemical, a class, or a response category.
3. State what counts as relevant evidence.
4. State what is outside the scope.
5. Record who approved the scope and the date.
6. Do not expand the scope halfway through extraction without recording a new
   decision.

**Done when:** the inclusion rules make it clear whether a source belongs in
the project.

### E2 — Create the context template
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2B — environmental evidence.

**Output:** `data/gene_sets/environment_evidence.tsv` template.

1. Copy stable ID, symbol, original ID, source, version, and mapping fields
   from the risk template.
2. Add exposure, dose or duration when available, species, tissue, assay,
   response type, and context.
3. Define what to write when context is missing.
4. Keep measured response separate from indirect association.
5. Review the required fields against two example rows.

**Done when:** an evidence row says what exposure and biological context it
represents.

### E3 — Record sources and extract the first one
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2B — environmental evidence.

**Output:** first environmental evidence block.

1. Record source name, release or publication date, and access date.
2. Copy original gene labels exactly.
3. Record exposure and context for every extracted row.
4. Record the response measurement or evidence type.
5. Mark missing mapping information for E5–E6.
6. Count extracted rows.

**Done when:** the first source can be independently audited.

### E4 — Extract remaining environmental sources
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2B — environmental evidence.

**Output:** completed raw environmental evidence table.

1. Use the same template for each source.
2. Keep sources separate in the source column.
3. Do not merge genes before preserving context.
4. Mark indirect or weak evidence instead of silently treating it as equal.
5. Record sources that produced no usable rows and why.

**Done when:** all approved sources have rows or documented empty results.

### E5 — Normalize stable IDs
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2B — environmental evidence.

**Output:** mapped environmental evidence.

1. Use the same namespace and reference release as G6.
2. Map original identifiers to stable IDs.
3. Keep original IDs and symbols.
4. Record mapping method.
5. Flag unresolved and one-to-many mappings.
6. Review ambiguous mappings against the evidence source.

**Done when:** the environment table can join to the risk table by stable ID.

### E6 — Record tissue and context mismatches
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2B — environmental evidence.

**Output:** context review note.

1. Filter rows with different species, tissue, assay, or response context.
2. Decide whether each mismatch is allowed by E1.
3. Keep allowed mismatches labeled.
4. Exclude disallowed rows with a reason.
5. Explain that context affects interpretation.

**Done when:** a reader can see where environmental evidence is indirect.

### E7 — Deduplicate without losing context
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2B — environmental evidence.

**Output:** deduplicated environment membership plus raw evidence.

1. Group membership by stable gene ID.
2. Create one membership row per gene/context combination if context matters.
3. Keep all supporting evidence rows.
4. Keep exposure labels and tissue fields.
5. Count genes and contexts before and after merging.

**Done when:** duplicate membership is removed but exposure context remains.

### E8 — Audit and export the environmental set
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2B — environmental evidence.

**Output:** `data/gene_sets/environment_response.tsv`.

1. Check required headers.
2. Check stable IDs and source fields.
3. Check exposure and context fields.
4. Check duplicate handling.
5. Check a sample of the rows against the evidence sources.
6. Record set counts.

**Done when:** the output passes the same provenance standard as G10.

### E9 — Compute overlap and candidate union
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2B — environmental evidence.

**Output:** `overlap.tsv` and `candidate_genes.tsv`.

1. Load the risk membership by stable ID.
2. Load environmental membership by stable ID.
3. Compute the intersection.
4. Compute the union.
5. Add `is_risk`, `is_environment`, and `is_overlap` flags.
6. Keep risk-only and environment-only genes distinguishable.
7. Save the source fields for the overlap.

**Done when:** every overlap row can be reproduced from the two input sets.

### E10 — Check set sizes and categories
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2B — environmental evidence.

**Output:** a category count note.

1. Count risk-only genes.
2. Count environment-only genes.
3. Count overlap genes.
4. Count the union.
5. Check that the categories add up correctly.
6. Investigate unexpected zero or enormous counts.
7. Do not discard the separate sets if the overlap is small.

**Done when:** category counts and joins are internally consistent.

### E11 — Write the interpretation boundary
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 2B — environmental evidence.

**Output:** a methods paragraph.

1. State that environmental evidence prioritizes genes initially.
2. State whether patient-level exposure data exist.
3. If they do not, state that exposure causation cannot be claimed.
4. Add this boundary to the methods and presentation notes.
5. Check the wording against the evidence boundary.

**Done when:** no figure or result caption implies exposure causation without
patient-level exposure measurements.

## Phase 3 — RNA-seq and data tasks

These tasks prepare the inputs. The mathematical definition must stay fixed even
when the data are inconvenient.

### D1 — Assess dataset suitability
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 3 — RNA-seq data.

**Output:** a dataset assessment table.

1. List possible public datasets.
2. Record tissue type for cases and controls.
3. Record sample count and subject count.
4. Record whether case/control labels are available.
5. Record library type and sequencing information.
6. Record whether junction counts, BAMs, or FASTQs are available.
7. Record metadata completeness and batch information.
8. Flag any tissue or design mismatch.
9. Check suitability against the project requirements.

**Done when:** the team can explain why the selected dataset can answer the
case/control comparison.

### D2 — Select and document the dataset
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 3 — RNA-seq data.

**Output:** an accession decision in `docs/data_methods.md`.

1. Compare candidates using D1.
2. Choose the dataset that satisfies the approved tissue/design requirements.
3. Record exact accessions and source pages.
4. Record access restrictions and expected file sizes.
5. Record the genome build and annotation if known.
6. Record why the chosen dataset was selected.
7. Do not download or analyze data before this decision is reviewed.

**Done when:** the same dataset can be located from the recorded accession.

### D3 — Obtain data reproducibly
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 3 — RNA-seq data.

**Output:** raw-data manifest and checksums.

1. Download only the approved files.
2. Keep large raw files outside Git unless the team explicitly approves them.
3. Record the source URL or accession for every file.
4. Record the download date.
5. Record a checksum if practical.
6. Record software versions used to unpack or convert files.
7. Do not rename files without recording the original name.

**Done when:** the recorded manifest identifies exactly which raw files were used.

### D4 — Build the sample manifest
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 3 — RNA-seq data.

**Output:** `data/metadata/samples.tsv`.

1. Give every sample a unique `sample_id`.
2. Record subject ID separately from sample ID.
3. Record condition as case or control.
4. Record tissue, accession, library type, batch, and file path.
5. Record genome build and annotation version.
6. Mark unknown metadata explicitly.
7. Check that no sample ID appears twice.
8. Check that every junction file maps to one sample.

**Done when:** the manifest tells the Rust tool which sample each row belongs to.

### D5 — Freeze reference conventions
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 3 — RNA-seq data.

**Output:** reference section in `docs/data_methods.md`.

1. Record genome build.
2. Record annotation release.
3. Record chromosome and coordinate conventions.
4. Decide whether graph nodes represent exons or splice sites.
5. Record strand handling.
6. Use the same conventions in node tables and junction tables.

**Done when:** the same endpoint always receives the same node ID.

### D6 — Obtain or derive junction counts
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 3 — RNA-seq data.

**Output:** raw junction-count files or a reproducible derivation record.

1. Check whether the dataset already provides junction counts.
2. If counts are provided, document their columns and software origin.
3. If counts must be derived, record the selected established processing tool.
4. Record its version and command settings.
5. Keep sample identity attached during processing.
6. Save a small test output before processing everything.

**Done when:** every count can be traced to a sample and a documented process.

### D7 — Check back-splice support
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 3 — RNA-seq data.

**Output:** a detection note.

1. Find out whether the input process can detect back-splice junctions.
2. Record any filtering that removes them.
3. Record whether strand and alignment information are available.
4. Explain that failure to detect a back-splice is not proof that none exists.
5. Record the limitation in the math and presentation notes.

**Done when:** the team knows what “back-splice absent” means in this dataset.

### D8 — Perform junction QC
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 3 — RNA-seq data.

**Output:** `qc/junction_qc.tsv` and `qc/data_qc.md`.

1. Count rows before filtering.
2. Check for missing sample, gene, endpoint, and count values.
3. Check for negative or non-finite counts.
4. Check duplicate rows and duplicate directed edges.
5. Check low-support junctions.
6. Check samples with unusually low total support.
7. Check endpoint and gene mapping failures.
8. Record every exclusion rule.
9. Save counts before and after each filter.

**Done when:** the filtered table can be reproduced from the QC note.

### D9 — Apply count and coverage rules
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 3 — RNA-seq data.

**Output:** filtered junction data and parameter record.

1. Use the threshold recorded in the approved statistics and math notes.
2. Apply the same rule to cases and controls.
3. Record the exact threshold.
4. Record how many rows and samples were removed.
5. Keep an exclusion table.
6. Distinguish missing measurement from measured zero support.

**Done when:** filtering cannot be changed silently after seeing the results.

### D10 — Create node and junction tables
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 3 — RNA-seq data.

**Output:** `data/processed/nodes.tsv` and `data/processed/junctions.tsv`.

1. Map each endpoint to a stable node ID.
2. Record endpoint coordinates and representation.
3. Copy sample ID, gene ID, source endpoint, destination endpoint, and count.
4. Check all gene IDs against the candidate set when filtering is active.
5. Check that every endpoint appears in the node table.
6. Check that counts remain numeric and nonnegative.

**Done when:** the Rust parser can join every junction to its sample, gene,
and endpoints.

### D11 — Check graph comparability
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 3 — RNA-seq data.

**Output:** comparability report.

1. Compare node and edge support across samples.
2. List genes present only in one group.
3. Decide how shared support or varying support will be handled.
4. Record the rule before case/control statistics.
5. Check the rule against two independent example comparisons.

**Done when:** the comparison rule does not depend on which group had the larger
observed result.

### D12 — Export an integration subset
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 3 — RNA-seq data.

**Output:** a small test subset and the full standardized candidate input.

1. Select a small number of valid samples and genes for integration testing.
2. Copy rows without changing their schema.
3. Include the sample manifest needed to interpret them.
4. Run the Rust parser against the subset.
5. Fix schema problems before scaling to the full candidate set.
6. Save the full input only after the subset passes.

**Done when:** the small real-data input runs without manually editing rows.

## Phase 4 — Mathematics tasks

These tasks define the expected answers. The first implementation uses
graph-level decomposition. The optional `B2`/full
simplicial extension is not a prerequisite for the MVP.

### M1 — Fix node and edge ordering
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** ordering rules in `docs/math_spec.md`.

1. Decide how node names are sorted or assigned.
2. Decide how edges are sorted.
3. State that the same ordering is used in vectors and matrices.
4. Write a two-edge example with its exact order.
5. Check that the implementation can follow the rule.

**Done when:** shuffled input produces the same matrix and metrics.

### M2 — Fix edge orientation
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** an orientation example.

1. Draw `E1 -> E2`.
2. Put `-1` in the `E1` row of that edge's B1 column.
3. Put `+1` in the `E2` row.
4. Put zero in every other row.
5. Explain that reversing an edge reverses its column and signal coordinate.
6. Test the sign with a hand-written matrix.

**Done when:** the team agrees on source and destination signs.

### M3 — Fix normalization and exceptional inputs
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** normalization section and error rules.

1. Define total support as the sum of counts for one sample and gene.
2. Define each normalized edge value as count divided by that total.
3. Decide what happens when total support is zero.
4. Decide what happens for an empty group.
5. Decide what happens to duplicate edges.
6. Decide what happens to self-loops.
7. Decide how disconnected components are handled.

**Done when:** the parser and the math fixtures agree on every edge case.

### M4 — Write the gradient least-squares problem
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** equation and dimensions in `docs/math_spec.md`.

1. Define `F` as an edge-length vector.
2. Define `phi` as a node-length vector.
3. Write `F_grad = B1^T phi`.
4. Write the minimization of `||F - B1^T phi||^2`.
5. Derive the normal equation or choose a trusted solver method.
6. Explain that node potentials are not unique up to a component constant.

**Done when:** a beginner can identify the dimensions of every matrix.

### M5 — Define the cycle projection
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** cycle-space definition and checks.

1. Define the cycle component as `F - F_grad` for the orthogonal projection.
2. State that it should satisfy `B1 F_cycle = 0` within tolerance.
3. State that a residual must not be called cyclic unless this check passes.
4. Explain the difference between directed biological cycles and graph cycle
   space.
5. Give one small example.

**Done when:** the software output uses the term consistently.

### M6 — Define energy and cycle fraction
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** formula and zero-signal rule.

1. Define signal energy as the squared Euclidean norm of `F`.
2. Define cycle energy as the squared norm of `F_cycle`.
3. Define the cycle fraction as cycle energy divided by signal energy.
4. Decide what status is returned for zero signal.
5. State that the value is a mathematical summary, not a causal claim.

**Done when:** the denominator and interpretation are documented.

### M7 — Choose tolerance and solver rules
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** numerical settings.

1. Choose floating-point precision.
2. Choose absolute and relative tolerances.
3. Define what counts as solver failure.
4. Define what happens if a matrix is singular.
5. Define convergence or pivot checks if an iterative or elimination solver is used.
6. Save settings with real run manifests later.

**Done when:** two implementations would reject the same failed fixture.

### M8 — Verify decomposition identities
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** validation checklist.

1. Reconstruct `F_grad + F_cycle`.
2. Compare it to `F` using the scale-aware tolerance.
3. Multiply `B1 F_cycle`.
4. Check that the result is near zero.
5. Compute the dot product of gradient and cycle components.
6. Check orthogonality.
7. Check energy additivity.
8. Check that the cycle fraction is between zero and one up to roundoff.

**Done when:** every decomposition reports pass/fail checks.

### M9 — Check incidence dimensions and cycle rank
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** hand-worked matrix examples.

1. Count nodes `n` and edges `m`.
2. Check that B1 has `n` rows and `m` columns.
3. Count connected components `k`.
4. Compute expected cycle rank `m - n + k`.
5. Compare the expected rank with a null-space calculation on fixtures.

**Done when:** the graph and matrix dimensions agree on every fixture.

### M10 — Validate a line and branching tree
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** linear and branching fixtures with expected results.

1. Draw `E1 -> E2 -> E3 -> E4`.
2. Add a branching tree without closing a graph loop.
3. Assign simple counts.
4. Normalize the counts.
5. Calculate the expected cycle fraction.
6. Confirm it is approximately zero.

**Done when:** trees do not produce a false cycle component.

### M11 — Validate exon skipping
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** exon-skipping fixture and independent expectation.

1. Add `E1 -> E2`, `E2 -> E3`, and `E1 -> E3`.
2. Record that the directed graph has no directed return path.
3. Calculate B1 and the projection independently.
4. Do not assume the cycle-space component must be zero.
5. Compare the implementation to the independent result.

**Done when:** the team understands why an acyclic directed graph may still
have graph cycle-space structure.

### M12 — Validate back-splicing
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** directed cycle fixture and expected result.

1. Create `E1 -> E2`, `E2 -> E3`, `E3 -> E1`.
2. Use equal weights first.
3. Calculate B1 by hand.
4. Verify that the equal circulation is in `ker(B1)`.
5. Verify that its gradient component is zero within tolerance.
6. Verify that its cycle fraction is one within tolerance.

**Done when:** the Rust fixture reproduces the hand result.

### M13 — Test noise
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** fixed-seed noise fixtures and a noise note.

1. Choose a random seed and save it.
2. Add small perturbations to a known signal.
3. Keep the original graph structure documented.
4. Measure how the cycle fraction changes.
5. Try a larger perturbation separately.
6. Explain when noise becomes indistinguishable from a signal change.

**Done when:** the noisy fixture can be reproduced from its saved seed.

### M14 — Test changing edge weights
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** weight-sweep table.

1. Keep the graph support fixed.
2. Change the back-splice weight through `0.01`, `0.05`, `0.10`, `0.20`,
   and `0.40` as proposed.
3. Record all other weights.
4. Record the normalization denominator.
5. Calculate the cycle fraction at each value.
6. Compare the values with independently computed projections.
7. Do not assume monotonicity without checking the exact normalization.

**Done when:** the table shows how the weighted signal affects the metric.

### M15 — Add edge-case fixtures
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** edge-case test inputs.

1. Create a zero-support input.
2. Create a disconnected graph.
3. Create a duplicate edge.
4. Create a self-loop if the contract rejects it.
5. Create a graph with a different row order.
6. Create a graph with reversed orientation and reversed signal coordinate.
7. Record the expected status for each.

**Done when:** errors are explicit and valid reordered inputs agree.

### M16 — Review reference expectations
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 4 — mathematics.

**Output:** signed-off math fixture folder.

1. Put every fixture and expected value in the shared fixture folder.
2. Explain which values were calculated by hand or an independent method.
3. Mark approximate values with their tolerance.
4. Resolve disagreements before copying expectations into Rust tests.
5. Record the check date and any corrections.

**Done when:** tests are based on reviewed expectations rather than the code's
own output.

## Phase 5 — Rust implementation tasks

These tasks build a dependable small tool. Large files, parallel processing,
advanced statistics, and optional Hodge extensions
come later.

### R1 — Freeze the Rust interfaces
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** Rust interface section in `docs/data_contract.md`.

1. Record the settled package and command name: `splice-girl`.
2. Record the current direct command: `splice-girl <junctions.tsv>`.
3. Record the current start command: `splice-girl --start`.
4. Record the planned future commands separately as not implemented.
5. Record input columns and output fields.
6. Agree on how errors are printed and how exit codes work.

**Done when:** the data, math, and code boundaries are written in one place.

### R2 — Set up the Rust project
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** `Cargo.toml`, source folders, and build instructions.

1. Confirm Rust is installed with `rustc --version` and `cargo --version`.
2. Create or verify the Cargo package.
3. Record the Rust toolchain if the team pins one.
4. Create `src/` and `tests/` folders.
5. Add a `.gitignore` for `target/`.
6. Run `cargo test` before adding project logic.
7. Record how a new reader builds the project.

**Done when:** a clean checkout can compile the empty or basic project.

### R3 — Choose numerical machinery
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** solver decision in `docs/math_spec.md`.

1. List possible linear-algebra approaches.
2. Test the preferred approach on a singular graph Laplacian.
3. Test a disconnected graph.
4. Test the linear and cycle fixtures.
5. Check whether dependencies are acceptable for the course.
6. Record precision, pivoting, and failure behavior.
7. Do not choose speed over a failed mathematical identity.

**Done when:** the solver behavior and its checks are recorded.

### R4 — Build the CLI shell
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** argument handling in `src/main.rs`.

1. Accept a junction TSV path.
2. Accept `--start`.
3. Print a useful usage message for missing arguments.
4. Return a nonzero exit code for invalid input.
5. Keep the startup banner separate from analysis logic.
6. Make the interactive menu exit cleanly on `q` or end-of-input.

**Done when:** the executable starts, prints usage on a bad command, and
supports the one current analysis feature.

### R5 — Implement parser validation
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** `src/parser.rs` and parser tests.

1. Read the header line.
2. Check exact required column names.
3. Split data rows on tabs.
4. Check the number of fields.
5. Check nonempty identifiers.
6. Parse counts as finite numbers.
7. Reject negative counts.
8. Apply self-loop and duplicate rules from the contract.
9. Include line numbers in error messages.
10. Test a valid row and several invalid rows.

**Done when:** malformed data fails clearly and valid data is represented in
typed Rust values.

### R6 — Group rows and construct graphs
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** grouping and graph modules.

1. Group rows by `sample_id`.
2. Within each sample, group by `gene_id`.
3. Confirm samples are not silently pooled.
4. Create a deterministic node list.
5. Create a deterministic edge list.
6. Preserve raw counts.
7. Calculate normalized signal values.
8. Reject zero-total groups as specified.
9. Test that shuffled rows produce the same graph.

**Done when:** one sample/gene group produces a stable graph and signal vector.

### R7 — Implement B1 and projection
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** matrix and Hodge modules.

1. Allocate an `n` by `m` matrix.
2. Put `-1` at each edge's source.
3. Put `+1` at each edge's destination.
4. Calculate the fitted gradient signal.
5. Calculate the cycle residual.
6. Handle disconnected components according to M3/M4.
7. Apply tolerance checks from M7/M8.
8. Return a clear failure instead of silently emitting bad values.

**Done when:** the implementation matches the reviewed M10 and M12 results.

### R8 — Connect one-sample analysis to output
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** per-sample metric output.

1. Run parser, grouping, graph construction, and decomposition in order.
2. Print sample and gene identifiers.
3. Print node and edge counts.
4. Print signal, gradient, and cycle energies.
5. Print cycle fraction.
6. Print per-edge components.
7. Print validation status.
8. Return failure if checks fail.

**Done when:** a beginner can run one fixture and see the values described in
the math specification.

### R9 — Pass the synthetic end-to-end gate
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** Rust tests and a short run note.

1. Run linear fixture.
2. Run branching fixture when M10 is ready.
3. Run exon-skipping fixture.
4. Run back-splice fixture.
5. Run noise fixture.
6. Run weight-sweep fixture.
7. Confirm every check passes.
8. Save command examples and expected key values.

**Done when:** no real biological result is claimed before all required synthetic
math checks pass.

### R10 — Implement one-gene case/control reading
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** comparison input handling.

1. Read the sample manifest.
2. Verify case/control labels.
3. Load per-sample files without pooling them.
4. Select one gene.
5. Report samples included and excluded.
6. Calculate one metric per eligible sample.
7. Test with a tiny known fixture.

**Done when:** one gene can be compared while preserving sample identity.

### R11 — Connect statistics outputs
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** statistics module interface.

1. Decide which per-sample columns the statistics checks need.
2. Export one row per sample/gene.
3. Include missing and excluded statuses.
4. Include run ID and parameter information.
5. Run the statistics checks on a tiny output.
6. Do not implement a p-value without an approved statistical specification.

**Done when:** the statistics checks can run without parsing terminal text.

### R12 — Implement ranking and junction explanations
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** deterministic ranking and component export.

1. Use the ranking rule recorded in the statistics and evidence notes.
2. Keep effect size separate from significance.
3. Join candidate-gene category flags.
4. Sort ties deterministically.
5. Select a documented junction contribution summary.
6. Preserve the sign of gradient and cycle components.
7. Export ranked genes and explanatory junctions.

**Done when:** a ranked row can be traced back to a sample, gene, and junction.

### R13 — Add result exports and manifests
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** machine-readable result files and `run_manifest.json`.

1. Define output headers.
2. Write sample metrics to TSV.
3. Write junction components to TSV.
4. Write comparisons to TSV.
5. Record input paths and checksums.
6. Record tool version and command.
7. Record solver settings, thresholds, and random seed.
8. Record exclusions and failures.

**Done when:** the run manifest identifies exactly how a result was produced.

### R14 — Handle errors safely
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** error tests and user documentation.

1. Test missing files.
2. Test malformed headers.
3. Test invalid counts.
4. Test missing sample metadata.
5. Test zero-support groups.
6. Test solver failure.
7. Test insufficient comparison samples.
8. Make each error name the problem and likely fix.
9. Avoid silently overwriting output files.

**Done when:** a beginner can recover from common input mistakes.

### R15 — Measure performance
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** small benchmark note.

1. Choose a small and a candidate-sized input.
2. Measure runtime with the same machine and command.
3. Measure memory if practical.
4. Identify the slowest step.
5. Only add sparse or parallel processing if profiling justifies it.
6. Confirm optimizations preserve the synthetic answers.

**Done when:** the team knows whether performance is a real blocker.

### R16 — Document installation and examples
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** installation section in README and `examples/`.

1. State the required Rust version.
2. Document `cargo test`.
3. Document `cargo run -- --start`.
4. Document `cargo install --path .` for local installation.
5. Explain the TSV schema.
6. Include a tiny example file.
7. Explain expected cycle fractions for the synthetic fixtures.
8. Document current limitations.

**Done when:** a new reader can install and run the tool from a clean checkout.

### R17 — Perform a clean-checkout review
**Source:** [Cargo project](../Cargo.toml), [Rust source](../src/), and [ROADMAP.md](../ROADMAP.md), Phase 5 — implementation.

**Output:** independent reproducibility note.

1. Ask someone who did not write the code to clone or receive the repository.
2. Have them install Rust if needed.
3. Have them run the documented build command.
4. Have them run tests.
5. Have them run `splice-girl --start`.
6. Have them run the synthetic analysis.
7. Record every confusing instruction.
8. Fix documentation or error messages.

**Done when:** the example can be run from the instructions without developer
help.

## Phase 6 — Real-data integration tasks

### D/R6.1 — Start with a small real subset
**Source:** [ROADMAP.md](../ROADMAP.md), Phase 6 — real-data integration.

**Output:** integration run folder.

1. Select a few eligible samples.
2. Select a small candidate-gene subset.
3. Copy the approved metadata.
4. Run the parser.
5. Compare parsed counts with the source table.
6. Inspect one graph manually.
7. Inspect one decomposition's checks.
8. Record exclusions and failures.

**Done when:** one real gene flows through the tool without manual file edits.

### D/R6.2 — Expand only after review
**Source:** [ROADMAP.md](../ROADMAP.md), Phase 6 — real-data integration.

1. Review numerical checks.
2. Review sample joins.
3. Review candidate-gene joins.
4. Fix issues found in the small subset.
5. Expand to all approved candidate genes.
6. Save the run manifest.

**Done when:** the full candidate input runs reproducibly and failure counts are
known.

## Phase 7 — Sample comparison tasks

### C1 — Freeze eligible samples
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 7 — sample comparison.

**Output:** eligible-sample list.

1. Read the sample manifest.
2. Apply tissue eligibility rules.
3. Apply QC rules.
4. Record case and control counts.
5. Record excluded samples and reasons.
6. Freeze the list before looking at final rankings.

**Done when:** the same sample list is used by every comparison.

### C2 — Calculate per-sample values
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 7 — sample comparison.

1. Choose one gene.
2. Collect its cycle fraction for every eligible sample.
3. Check that each sample appears once.
4. Keep missing values marked as missing.
5. Do not replace missing values with zero.
6. Repeat for every candidate gene.

**Done when:** `results/sample_metrics.tsv` has one auditable row per sample/gene.

### C3 — Calculate descriptive group effects
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 7 — sample comparison.

1. Separate case and control values.
2. Count values in each group.
3. Calculate group means or the pre-approved summary.
4. Calculate case minus control difference.
5. Record variation and missingness.
6. Do not call a large difference significant yet.

**Done when:** every effect has a direction and sample count.

### C4 — Produce provisional ranks
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 7 — sample comparison.

1. Apply the documented effect-ranking rule.
2. Break ties deterministically.
3. Keep p-value and q-value columns empty or marked unavailable until S tasks.
4. Save the ranking rule with the table.

**Done when:** rankings are descriptive and cannot be mistaken for inference.

## Phase 8 — Statistics and baseline tasks

### S1 — Identify the biological unit
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 8 — statistics and baselines.

1. Determine whether samples are independent, paired, or repeated.
2. Identify the subject ID column.
3. Record any paired structure.
4. Ask whether label permutations are allowed.
5. Write the restriction in `docs/statistics_spec.md`.

**Done when:** the team knows what can and cannot be shuffled.

### S2 — Freeze analysis rules
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 8 — statistics and baselines.

1. Freeze eligible samples.
2. Freeze eligible genes.
3. Freeze support and coverage rules.
4. Freeze missing-value handling.
5. Freeze graph-support comparison policy.
6. Date and sign the specification.

**Done when:** rules cannot be changed after seeing the top results without a
documented sensitivity analysis.

### S3 — Validate per-sample input
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 8 — statistics and baselines.

1. Check one row per sample/gene.
2. Check sample IDs against the manifest.
3. Check group labels.
4. Check metric ranges.
5. Count missing and excluded rows.
6. Stop if patients were accidentally pooled.

**Done when:** the statistics checks pass on the input table.

### S4 — Define the primary effect
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 8 — statistics and baselines.

1. Write the case summary statistic.
2. Write the control summary statistic.
3. Define the difference as case minus control if that remains approved.
4. Record how variation is shown.
5. Record how missing data are handled.

**Done when:** the same effect can be calculated from the TSV.

### S5 — Specify permutation testing
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 8 — statistics and baselines.

1. Define the observed effect.
2. Define which labels may be shuffled.
3. Choose the number of permutations.
4. Choose and save a random seed.
5. Define one- or two-sided testing.
6. Define the p-value calculation.
7. Respect paired or blocked designs.

**Done when:** the permutation procedure can be written as numbered instructions.

### S6 — Validate the permutation implementation
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 8 — statistics and baselines.

1. Create a no-difference fixture.
2. Create a known-difference fixture.
3. Test missing values.
4. Test a small group where all permutations can be enumerated.
5. Compare exact and sampled results when possible.
6. Check that the saved seed reproduces the result.

**Done when:** the code does not produce a zero p-value merely because too few
permutations were requested.

### S7 — Run tests only when design permits
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 8 — statistics and baselines.

1. Check replication and exchangeability.
2. If inference is valid, run the approved tests.
3. If inference is not valid, report descriptive effects only.
4. Keep an explicit “inference unavailable” status.
5. Do not manufacture biological replication by treating junctions as patients.

**Done when:** every reported p-value has a defensible sample design.

### S8 — Apply FDR correction
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 8 — statistics and baselines.

1. Define the family of genes being tested.
2. Collect raw p-values for that family.
3. Apply the approved correction, such as Benjamini–Hochberg.
4. Keep raw p-values and q-values.
5. Mark untested genes explicitly.
6. Record the correction family and method.

**Done when:** a reader can tell which genes were included in the correction.

### S9 — Calculate baselines
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 8 — statistics and baselines.

1. Define junction-usage differences on the same eligible samples.
2. Calculate graph cycle rank.
3. Include back-splice support only if detection is valid.
4. Include expression comparison only if compatible expression data exist.
5. Apply the same filtering and sample rules.

**Done when:** the baseline comparison uses the same data restrictions as the
Hodge summary.

### S10 — Compare methods
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 8 — statistics and baselines.

1. Compare effect sizes.
2. Compare rankings.
3. Compare sensitivity to support thresholds.
4. Compare cycle fraction with cycle rank.
5. Record whether Hodge decomposition adds information.
6. Accept a null or no-advantage result as a valid result.

**Done when:** the report does not claim novelty merely because the method is
more complicated.

### S11 — Run sensitivity checks
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 8 — statistics and baselines.

1. List the approved alternate thresholds.
2. Re-run without changing the biological question.
3. Record changes in eligible samples and genes.
4. Compare effect directions and top ranks.
5. Keep all attempted settings in the run record.

**Done when:** readers can see whether conclusions depend on one arbitrary filter.

### S12 — Export the statistics report
**Source:** [math specification](../docs/math_spec.md) and [ROADMAP.md](../ROADMAP.md), Phase 8 — statistics and baselines.

1. Export effects, counts, p-values, q-values, and statuses.
2. Export baseline results.
3. Export permutation settings and seed.
4. Write a plain-language methods paragraph.
5. Write a plain-language result paragraph.
6. State limitations and unavailable inference.

**Done when:** two readers can use the report without
guessing what a column means.

## Phase 9 — Interpretation tasks

### I1 — Join results to provenance
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 9 — interpretation.

1. Load the final gene results.
2. Join risk/environment/overlap flags by stable ID.
3. Check that no gene silently lost its evidence.
4. Link important junctions to graph and sample output.
5. Save an interpretation table.

**Done when:** every discussed gene has both computational and biological
provenance.

### I2 — Review candidate categories
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 9 — interpretation.

1. Summarize risk-only genes.
2. Summarize environment-only genes.
3. Summarize overlap genes.
4. Compare counts and effects descriptively.
5. Avoid claiming the overlap is automatically the most important group.

**Done when:** category comparisons reflect the actual candidate set.

### I3 — Write limitations
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Phase 9 — interpretation.

1. State that graph cycles do not prove circular RNA.
2. State alignment, coverage, and ambiguous-read limitations.
3. State that topology differences do not prove disease causation.
4. State that environmental-response genes do not prove exposure causation.
5. State whether Hodge summaries outperform baselines.
6. State sample-size and tissue limitations.

**Done when:** every strong claim in the report has its matching boundary.

## Phase 10 — Tool cleanup and usability tasks

### U1 — Make installation repeatable
**Source:** [README.md](../README.md), [Cargo project](../Cargo.toml), and [ROADMAP.md](../ROADMAP.md), Phase 10 — usability.

1. Test a clean checkout.
2. Record Rust prerequisites.
3. Test `cargo build`.
4. Test `cargo install --path .`.
5. Test the installed `splice-girl` command.
6. Document common PATH problems.

**Done when:** a new user can install without editing source code.

### U2 — Make the input contract readable
**Source:** [README.md](../README.md), [Cargo project](../Cargo.toml), and [ROADMAP.md](../ROADMAP.md), Phase 10 — usability.

1. Show the header.
2. Explain every column.
3. Include a valid example.
4. Include an invalid example.
5. Explain sample versus gene grouping.
6. Explain missing versus zero.

**Done when:** the documentation can be handed to a beginner.

### U3 — Make errors actionable
**Source:** [README.md](../README.md), [Cargo project](../Cargo.toml), and [ROADMAP.md](../ROADMAP.md), Phase 10 — usability.

1. Include the filename.
2. Include the line number when relevant.
3. Name the bad column.
4. Explain the expected value.
5. Give a likely correction.
6. Return a nonzero status for failure.

**Done when:** an error message tells a beginner what to inspect next.

### U4 — Add reproducibility metadata
**Source:** [README.md](../README.md), [Cargo project](../Cargo.toml), and [ROADMAP.md](../ROADMAP.md), Phase 10 — usability.

1. Assign a run ID.
2. Record command and input paths.
3. Record version and commit.
4. Record thresholds and tolerances.
5. Record seed and exclusions.
6. Write a manifest next to results.

**Done when:** a result can be traced to one software state and one input state.

### U5 — Independent usability review
**Source:** [README.md](../README.md), [Cargo project](../Cargo.toml), and [ROADMAP.md](../ROADMAP.md), Phase 10 — usability.

1. Give the instructions to someone unfamiliar with the code.
2. Ask them to install and run the synthetic example.
3. Watch where they get stuck.
4. Record their exact question.
5. Improve the README or error message.
6. Repeat once.

**Done when:** the workflow can be completed from the written instructions.

## Phase 11 — Visualization and presentation tasks

### V1 — Draw the complete pipeline
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [README.md](../README.md), and [ROADMAP.md](../ROADMAP.md), Phase 11 — presentation.

1. Start with genetic and environmental evidence.
2. Draw candidate-gene categories.
3. Draw RNA-seq junction data.
4. Draw standardized tables.
5. Draw the Rust graph engine.
6. Draw gradient/cycle decomposition.
7. Draw sample comparison and ranking.
8. Draw reusable tool output.
9. Label which steps are completed and which are future.

**Done when:** a beginner can explain the project by following the arrows.

### V2 — Draw a splice graph
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [README.md](../README.md), and [ROADMAP.md](../ROADMAP.md), Phase 11 — presentation.

1. Label nodes as exons or splice sites.
2. Draw directed arrows.
3. Label example support values.
4. Show an alternative junction.
5. Show a back-splice only as an example if its meaning is explained.
6. Match orientation to B1.

**Done when:** the figure agrees with the actual synthetic fixture.

### V3 — Show decomposition components
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [README.md](../README.md), and [ROADMAP.md](../ROADMAP.md), Phase 11 — presentation.

1. Show the original signal.
2. Show the gradient component.
3. Show the cycle component.
4. Keep edge order the same in all three views.
5. Explain that components may be signed.
6. Label the figure as mathematical output.

**Done when:** no viewer has to infer which edge belongs to which value.

### V4 — Show one-gene case/control output
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [README.md](../README.md), and [ROADMAP.md](../ROADMAP.md), Phase 11 — presentation.

1. Choose a gene only after data review.
2. Show individual sample points or rows.
3. Show case and control labels.
4. Show group summaries.
5. Show exclusions.
6. Include the cycle fraction and its units or interpretation.

**Done when:** sample variation is visible instead of only two averages.

### V5 — Prepare top-result figures
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [README.md](../README.md), and [ROADMAP.md](../ROADMAP.md), Phase 11 — presentation.

1. Run the statistics checks and record their result.
2. Select results using a documented rule.
3. Include effect size and uncertainty where available.
4. Include p/q status only when valid.
5. Link each gene to its provenance category.
6. Mark illustrative examples clearly.

**Done when:** every plotted value comes from a versioned result table.

### V6 — Write biological explanations
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [README.md](../README.md), and [ROADMAP.md](../ROADMAP.md), Phase 11 — presentation.

1. Describe what a junction means.
2. Describe what a cycle fraction means mathematically.
3. Explain what it does not prove.
4. Explain the environmental prioritization boundary.
5. Use plain language before equations.
6. Check the wording with the plain-language definitions above.

**Done when:** a basic biology/computer-science student can follow the result.

### V7 — Plan presentation speaking parts
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [README.md](../README.md), and [ROADMAP.md](../ROADMAP.md), Phase 11 — presentation.

1. List the problem introduction.
2. List gene prioritization.
3. List data preparation.
4. List the graph and math explanation.
5. List the CLI demonstration.
6. List statistics interpretation.
7. List limitations.
8. Write the handoff sentence between speakers.

**Done when:** every section has a written speaking plan and backup notes.

### V8 — Select a real demo gene
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [README.md](../README.md), and [ROADMAP.md](../ROADMAP.md), Phase 11 — presentation.

1. Wait until real-data QC and statistics review.
2. List eligible genes with valid output.
3. Choose using a recorded rule such as clear supported output.
4. Accept a null-result gene if it is the clearest demo.
5. Record the gene ID and reason.

**Done when:** the demo gene is reproducible and not selected only because it
looks exciting after repeated searching.

### V9 — Package the tiny demo
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [README.md](../README.md), and [ROADMAP.md](../ROADMAP.md), Phase 11 — presentation.

1. Copy only the small synthetic and one-gene inputs needed live.
2. Keep the schema valid.
3. Write exact commands.
4. Record expected output.
5. Measure runtime.
6. Remove machine-specific paths.

**Done when:** the demo finishes in seconds on the presentation machine.

### V10 — Rehearse and make a backup
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [README.md](../README.md), and [ROADMAP.md](../ROADMAP.md), Phase 11 — presentation.

1. Rebuild the tool on the presentation machine.
2. Run the synthetic command.
3. Run the real-gene command.
4. Time both commands.
5. Save terminal output and figures.
6. Make an offline recording or screenshot backup.
7. Practice explaining the backup as precomputed output.

**Done when:** the demonstration can continue without network access.

### V11 — Audit figures and captions
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [README.md](../README.md), and [ROADMAP.md](../ROADMAP.md), Phase 11 — presentation.

1. Link every figure to a run ID.
2. Check labels against source tables.
3. Check that examples are not called experimental results.
4. Check that cycle language does not imply circular RNA.
5. Check environmental language does not imply causation.
6. Have someone who did not make the figure proofread it.

**Done when:** captions are accurate without a verbal correction.

## Phase 12 — Final report and presentation tasks

These are shared tasks after the core outputs exist.

### F1 — Assemble methods
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [roadmap.tex](../roadmap.tex), and [ROADMAP.md](../ROADMAP.md), Phase 12 — final documents.

1. Describe candidate-gene evidence and provenance.
2. Describe dataset and sample eligibility.
3. Describe junction standardization.
4. Describe graph orientation and B1.
5. Describe normalization and projections.
6. Describe synthetic validation.
7. Describe comparisons and statistics.
8. Include versions and run identifiers.

### F2 — Assemble results
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [roadmap.tex](../roadmap.tex), and [ROADMAP.md](../ROADMAP.md), Phase 12 — final documents.

1. Report QC counts.
2. Report eligible samples and genes.
3. Report per-sample and group outputs.
4. Report statistics and baselines.
5. Report ranking and junction explanations.
6. Distinguish null and positive results.

### F3 — Assemble software release
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [roadmap.tex](../roadmap.tex), and [ROADMAP.md](../ROADMAP.md), Phase 12 — final documents.

1. Include source code.
2. Include installation instructions.
3. Include input schema.
4. Include synthetic fixtures.
5. Include commands and expected outputs.
6. Include limitations and version information.

### F4 — Compile and proof the documents
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [roadmap.tex](../roadmap.tex), and [ROADMAP.md](../ROADMAP.md), Phase 12 — final documents.

1. Compile the available LaTeX source.
2. Check for missing references.
3. Check table widths.
4. Check code blocks.
5. Check overfull and underfull boxes.
6. Open the generated PDF and inspect every page.
7. Confirm that figures and numbers match the final run.

### F5 — Final acceptance meeting
**Source:** [project proposal](../assets/hodgetheorynew.pdf), [roadmap.tex](../roadmap.tex), and [ROADMAP.md](../ROADMAP.md), Phase 12 — final documents.

1. Run the complete checklist.
2. Check each completed deliverable against its task description.
3. Record unresolved limitations.
4. Decide which optional extensions are deferred.
5. Freeze the version used in the presentation.

**Done when:** the report, presentation, tool, data record, and limitations all
refer to the same validated run.

## Optional Phase 13 tasks

These tasks are **STRETCH** work. Do not start them if they put the core tool,
real-data comparison, statistics, or presentation at risk.

### X1 — Patient-level exposure comparison
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Optional Phase 13.

1. Confirm that exposure measurements exist for the same patients.
2. Check that the exposure data can be joined to RNA-seq samples.
3. Define low/high or other groups before looking at results.
4. Reuse the approved graph and statistics method.
5. Record missing exposure values.
6. Keep the extension separate from the prioritization-only analysis.

### X2 — Full simplicial Hodge extension
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Optional Phase 13.

1. Define a biologically meaningful two-dimensional face.
2. Explain what each face means biologically.
3. Define `B2` and its orientation.
4. Validate boundary identities.
5. Add the curl and harmonic components only after review.
6. Compare with the graph-level MVP.

### X3 — Genome-wide scaling
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Optional Phase 13.

1. Estimate input size.
2. Profile the candidate-gene implementation.
3. Design memory limits and sparse structures.
4. Test on progressively larger subsets.
5. Preserve the same numerical checks.
6. Record which genes fail and why.

### X4 — Other-disease application
**Source:** [project proposal](../assets/hodgetheorynew.pdf) and [ROADMAP.md](../ROADMAP.md), Optional Phase 13.

1. Freeze the endometriosis release first.
2. Choose a compatible disease only if time remains.
3. Reuse the same data contract and graph method.
4. Document every disease-specific change.
5. Keep the extension separate from the primary report.

## First tasks to start this week

If you are a beginner and do not know where to begin, choose one of these
small tasks:

1. Complete P1.1 and P1.2.
2. Complete P1.4 using the existing synthetic TSV files.
3. Complete M1, M2, and M12 by hand.
4. Complete R4 and R5, then run `cargo test`.
5. After these are moving, start G1/G3, E1/E2, S1/S2, and V1 in any order.

These tasks can begin in parallel. Do not wait for real RNA-seq data to start
the synthetic math, parser, graph, documentation, or presentation work.

## Definition of a complete core project

The core project is complete when all of the following are true:

- Risk, environmental, overlap, and candidate gene tables have provenance.
- Case/control data have reviewed metadata and reproducible junction processing.
- Synthetic math fixtures pass reconstruction, kernel, orthogonality, and energy
  checks.
- Splice Girl reads the standardized input and preserves sample identity.
- Per-sample decomposition and case/control comparison outputs exist.
- Statistics and baselines are valid for the available sample design, or the
  report explicitly says inference is unavailable.
- The Rust tool exports machine-readable results and a run manifest.
- A beginner can install and run the synthetic example.
- One small real-gene demonstration and an offline backup are ready.
- The final report states limitations, null results, and unresolved decisions.

## What not to do yet

- Do not download a random dataset before D1 and D2 are reviewed.
- Do not choose genes based on an unrecorded preference.
- Do not call every residual a cycle without the M5 checks.
- Do not pool patients into one sample.
- Do not call environmental-response evidence proof of exposure causation.
- Do not call a graph cycle proof of circular RNA.
- Do not add `B2` before defining meaningful faces.
- Do not build the final statistical story from a single example gene.
- Do not replace the approved project with a different research question.
