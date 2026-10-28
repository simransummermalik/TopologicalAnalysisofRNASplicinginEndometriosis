# October 24, 2026 — Summer

**Status:** Pushed

**Commits:** `6a6fb87`, `7d22a2e`, and `cea7764`

I expanded the Splice Girl graph output so the same sample and gene can be
viewed as three matching panels: the original normalized junction signal, the
gradient part, and the cycle part. The panels use the same node positions and
edge order, which makes the mathematical split easier to compare. I updated
the demo instructions, input documentation, task list, priorities, README,
and example SVG so the feature is explained where someone would look for it.

I also made a large Mermaid project map and linked it at the very beginning of
the README. The first version was a normal pipeline chart, but I changed it to
look more like a SQL schema because that format shows the project pieces,
fields, and connections more clearly. The final map connects the study and
samples to preprocessing, Rust input, splice graphs, Hodge results, sample
comparisons, statistics, interpretation, figures, and reproducible outputs.

