# October 20, 2026 — Summer

**Status:** Pushed

**Commits:** `d22d4eb`, `a0fb1c2`, `4a9bb00`, `f374576`, `27813a1`, and
`637a7da`

I built the first plain HTML and CSS website for Splice Girl and added a graph
view to the Rust tool so a completed analysis can be opened as an SVG report.
The graph is created from the selected sample and gene, so its nodes, arrows,
read support, normalized signal, and cycle component come from the input data.
I then simplified the styling so it looked more like a basic academic report
and updated the README and demo instructions to explain how to generate it.

I also organized the first real RNA-seq quality-control results. An initial
MultiQC download did not contain the complete reports, so I moved it into an
`incorrect` folder instead of presenting it as valid. I then added the usable
control and endometriosis MultiQC HTML reports and their summary tables after
the FastQC collections were flattened. These reports document raw-read quality
for the four controls and four endometriosis samples; they are quality checks,
not splice-junction results.
