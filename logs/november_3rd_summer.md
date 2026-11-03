# November 3, 2026 — Summer

**Status:** Prepared locally; not committed or pushed

I added a report command so people can actually look through Splice Girl's
results without reading terminal output or opening every table separately.
It creates a plain academic HTML page with the gene comparisons, a point for
each individual sample, exact cycle fractions, and explanations for excluded
samples. Each eligible sample also has a link to its observed, gradient, and
cycle graph panels. The report works offline and comes with the same TSVs and
input-checksum record as the existing comparison command. I added it to the
startup menu as option 6, checked that its numbers match the original exports,
tested missing groups and unusual sample names, and checked the rendered page
in Chrome. The Rust suite now has 62 tests. The real biological analysis and
statistics still need the group's checked junction data and final rules.
