# September 24, 2026 — Summer

**Status:** Pushed

I wrote down the first data rules for Splice Girl and added them to
`docs/data_contract.md`. The contract explains the five columns the program
expects and what counts as a bad row, such as a negative count, a duplicate
edge, an empty value, or a self-loop. I also added parser tests for these bad
inputs so the program gives an error instead of quietly using wrong data. I
started the preprocessing notes from Dr. White's advice and added a short
plain-language explanation in the README for why the project uses Hodge
theory. These changes connected the biology notes, the math idea, and the
software rules in one place.

Commits: `7640a88`, `4db4516`, `dce1173`, `a686340`, and `934fcc2`.
