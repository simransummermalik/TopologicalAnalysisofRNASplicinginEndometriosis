# Running the full cohort outside this laptop

The eight-sample run needs the external drive folder or server address and
working folder. The user confirmed that no external drive or server is currently
available. No full-cohort download or alignment was launched; the instructions
below are for when a suitable location becomes available.

The saved QC tables contain **545,210,220 raw read pairs**. Extrapolating the
pilot's compressed reads, SAM and BAM sizes gives about **768 GB (716 GiB)**
for the current workflow if it retains every sample's raw reads, cleaned reads,
SAM and BAM, including a 25% allowance. This is a planning estimate, not a
minimum for every possible workflow. References and extra tool temporary files
need additional space. Sequential processing with a separately agreed retention
policy could use less; the current pipeline does not automatically delete data.
The [recorded preflight](../data/pilots/SRR36710387_prefix1m_nov17/cohort_preflight.json)
contains the calculation and exact pilot sizes.

## On a mounted drive

1. Create a folder for the project on the actual drive. Keep large reads,
   references, indexes and run outputs there, outside Git.
2. Save a TSV sample sheet with columns `sample_id`, `condition`, `read1`,
   `read2`. Use the eight IDs and conditions from
   [samples.tsv](../data/metadata/samples.tsv). Each mate must be a distinct
   real file; paths may be absolute or relative to the sheet.
3. Check availability and capacity before starting:

```sh
python3 preprocessing/cohort_readiness.py \
  --samples /Volumes/YOUR_DRIVE/splice-girl/samples.tsv \
  --storage /Volumes/YOUR_DRIVE/splice-girl \
  --output /tmp/splice-girl-cohort-preflight.json
```

`YOUR_DRIVE` is a placeholder, not an existing drive name. This command only
reads files and writes the requested small JSON report. It does not download
reads or start alignment. Without `--samples`, input presence is unassessed;
false presence fields do not mean the already cleaned local control disappeared.
The pilot sizing inputs must remain accessible, or be supplied explicitly with
`--pilot-manifest`, `--pilot-sam`, and `--pilot-bam`.

## On a server

Run the preflight on that server using server-local paths and its filesystem
capacity. A path on this Mac cannot measure a remote server's free space.
Transfer source and small manifests, install the documented tools, and verify
the reference/checksums before reusing an index. Keep a stable working session
and logs so a browser or SSH disconnect does not lose the job.

## Before starting the actual batch

Resolve the scientific settings in the [dataset review](dataset_review.md)
and [pilot audit](../data/pilots/SRR36710387_prefix1m_nov17/assignment_audit/README.md).
Then use the existing `batch_preprocessing.py plan` command with explicit
strandedness, anchors, support and reference paths; inspect that plan before
using `run`. The batch adapter validates every sample first and processes
sequentially. Its guide is in [real_pilot.md](real_pilot.md) and
[preprocessing/README.md](../preprocessing/README.md).

Raw samples need Cleaver cleaning; inputs already cleaned by the recorded
procedure use `--skip-cleaning`. Because that setting applies to the whole
batch, do not mix raw and cleaned inputs in one sheet. Preserve the cleaning
manifest and sample IDs when combining separately processed batches.

After processing, verify each BAM and converter QC, compare node identities
and eligibility across samples, and save the reviewed junction table. Only
then run the frozen four-versus-four statistics and baseline/sensitivity
commands. A prefix, a missing gene, or an excluded sample must not be filled
with invented zeros to make a statistical test run.
