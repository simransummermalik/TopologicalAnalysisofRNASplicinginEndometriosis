# October 5, 2026 — Summer

**Status:** Pushed

I checked Sydney's sample sheet and turned it into a cleaner sample manifest.
The new reader checks that every sample has one read-one file and one read-two
file, makes sure the paired rows agree, and gives the samples simpler condition
labels. It wrote a 48-sample E-MTAB-15117 manifest, removed the broken old
reader file, and updated the research README with the sample counts and the
tissue-location warning. This makes the source data easier to trace before any
raw reads are processed.

Commit: `01abbe2`.
