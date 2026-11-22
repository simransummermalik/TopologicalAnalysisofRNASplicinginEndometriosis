# November 21, 2026 — Summer

**Status:** HPC sample completed; this log is prepared locally. The run outputs remain in HPC scratch and have not been added to the repository.

I checked the completion of the SRR36710385 preprocessing run (Slurm job 27893687), using the same human-reference alignment and junction-counting settings as the earlier SRR36710386 and SRR36710387 runs. Slurm reported `COMPLETED`, exit code `0:0`, and an elapsed runtime of 2 hours 53 minutes 36 seconds. The run manifest ended with `complete`; `junctions.tsv` was 19 MB, `alignment_summary.txt` was 637 bytes, and `junction_qc.json` was about 1.4 KB. The remaining five samples—two cases and three controls—still need to be processed before a four-versus-four comparison can be made. The sample outputs are still on HPC scratch; no biological conclusion is claimed.
