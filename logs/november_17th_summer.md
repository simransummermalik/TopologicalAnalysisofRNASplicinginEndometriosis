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
