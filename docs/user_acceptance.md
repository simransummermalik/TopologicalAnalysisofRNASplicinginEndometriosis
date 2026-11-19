# First-time researcher acceptance exercise

**Status: independent human review pending.** Automated developer checks are
recorded separately in `docs/software_validation.md`. No tester feedback or
completion time has been invented.

Ask someone who has not built the tool to try these tasks using the installed
CLI and its written instructions. Record where they need help before giving
it. No personal or human-read data are required. Testing this example is not
validation of a biological finding.

## Tasks for the reader

1. Install the experimental source release using its crate README. Record
   your operating system, `rustc --version` and `splice-girl --version`.
2. Use command help to create the supplied example in a new folder.
3. Find the saved configuration and check it without starting analysis.
4. Run that configuration and open the resulting report.
5. Explain why TREE, LOW_COVERAGE and ABSENT are excluded from the primary
   statistical test, and identify the gene that remains eligible.
6. Find the primary cycle fraction for control S7. Explain why a zero-count
   shortcut can coexist with a nonzero cycle component on a shared graph.
7. Try the same output folder again. Does the error explain what to change?
8. Change `output` to a new folder, set `statistics` to null, check and run
   again. Can you locate the descriptive results and exclusions?
9. Find the license for the installed version. Explain how its license notices
   or attribution differ from the separate request to cite the software in a
   paper.

## Record actual observations

Use the repository's `user_testing/FEEDBACK_TEMPLATE.md` or a plain text file.
For each task, record whether it worked, any confusing wording, the exact
command/error, and how much assistance was needed. Include expected versus
observed behavior and a minimal reproduction for any bug. Do not add private
sample paths or identifying biological data to public feedback.

The author should correct reproducible issues, add targeted regressions for
behavioral bugs, and have the reader retry the affected steps. A real person
must complete this before R17/U5 can be marked done. Success on an automated
fresh install does not satisfy this independent review.

## Evaluator reference

The supplied primary setup has one eligible gene, TRIANGLE, across eight
invented samples. TREE has no cycle, LOW_COVERAGE lacks retained support in
S7, and ABSENT has no shared edges. The S7/TRIANGLE cycle fraction is 2/3 on
the common graph. Seven sensitivity scenarios retain every requested gene in
the audit; three skip statistics because no genes qualify. Existing output
folders are never replaced. Use these expectations to assess understanding,
not to claim the independent test has already happened.
