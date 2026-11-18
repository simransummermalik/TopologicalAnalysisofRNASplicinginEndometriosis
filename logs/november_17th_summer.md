# November 17, 2026 — Summer

**Status: prepared locally; no commit or push performed by the assistant.**

I connected more of the tool so we can go from local read files to a saved
Splice Girl report. The preprocessing code uses the actual Cleaver tool for
cleaning and BAM work, with HISAT2 for alignment, and it counts annotated
junctions while recording what gets filtered out. There is also a batch
command that handles samples one at a time and keeps their counts separate.
I added simpler graph/count comparisons beside the Hodge results, optional
gene sorting with consistent ties, and a way to rerun the analysis with
different support thresholds without turning missing results into zeros.
The original gene list and its evidence columns stay with the output. I added
one command to check the project and tested the actual preprocessing-to-report
path on tiny synthetic paired reads, along with 91 Rust and 41 Python tests.
I also fixed an occasional temporary-folder collision in the tests and updated
the guides and task statuses. The full human-data analysis, reviewed evidence
tables and final scientific settings still need to be completed; these tests
do not count as biological results.

## Real-data pilot and evidence joins

I built the full human HISAT2 index from our verified GENCODE files and ran
the first one million already-cleaned read pairs from our local control
sample through alignment and junction extraction. The alignment rate was
97.23%. A rule written before looking at the Hodge values selected CDH11,
and I checked its counts directly against the BAM and its math against a
separate NumPy calculation. Those checks passed, and I saved a basic report,
graph, count table and threshold checks in `data/pilots/`. This gene's retained
graph has no loops, so its cycle fraction is correctly zero; this is a working
real-data test, not a disease finding. The QC also showed a large number of
ambiguous annotation matches, which we need to review before processing the
whole group. I connected the existing gene-evidence labels directly to ranked
and statistical outputs, made long coordinate labels readable in the figures,
and added checks against incomplete indexes and mismatched paired reads. The
current combined suite passes 178 tests, including 95 Rust tests. No commit or
push was performed by the assistant.
