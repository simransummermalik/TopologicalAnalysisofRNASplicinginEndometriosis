# November 16, 2026 — Summer

**Status:** Prepared locally; not committed or pushed

I connected the permutation and multiple-testing correction code to a usable
`stats` command and added it as option 7 in the startup menu. It takes a gene
list, test direction, independent-sample declaration, and minimum sample count
explicitly, then saves effects, raw and adjusted p-values, permutation counts,
and reasons for untested genes. The offline report shows those results next
to the sample plots and Hodge graphs, and the run manifest records the gene
list checksum and settings. I tested known answers, missing genes, small
groups, repeated subjects, large permutation spaces, and the menu workflow.
The project now has 89 tests. R11 is complete as a software task; real-data
QC, study-design approval, baselines, and the final biological results still
need to be finished by the group.
