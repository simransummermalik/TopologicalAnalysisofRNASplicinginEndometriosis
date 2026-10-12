# September 30, 2026 — Summer

**Status:** Pushed

I built the first small preprocessing workflow for raw RNA sequencing reads. It
can check whether the needed programs are installed, print a plan before doing
work, run one sample, clean reads with the Cleaver-X command, optionally remove
PhiX reads with Bowtie2, align reads with STAR, index the BAM file, and convert
STAR splice junctions into the five-column file that Splice Girl reads. I also
added the junction-to-gene helper, small test fixtures, preprocessing notes,
and tests that check the commands without needing large biological files.

Commits: `afaf59b` and `6f57ea9`.
