# November 2, 2026 — Summer

**Status:** Prepared locally; not committed or pushed

I finished the software housekeeping that was still missing around the Rust
tool. I added a source map so the team can see which proposal files are local
and which parent Overleaf file still needs to be supplied, plus a dated
decisions log for the reference, cohort, graph conventions, and open research
choices. I also added a release checklist and a short changelog so the current
prototype is clearly separated from the unfinished real-data study.

The repository now has GitHub Actions checks for Rust formatting, tests,
Clippy, and the Python preprocessing tests. I added a repeatable benchmark that
creates temporary synthetic data at 200 genes and eight samples, runs the
all-gene comparison ten times, and reports timing without adding fake biology
to the repository. The benchmark and developer checks passed locally, along
with the Rust suite and all six preprocessing tests. Nothing was committed or
pushed during this work.
