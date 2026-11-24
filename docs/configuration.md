# Run from one saved configuration

The installed CLI accepts a strict JSON configuration. Unknown/duplicate
settings and wrong types fail instead of being silently ignored. All input
and output paths are relative to the configuration file's folder.

## Start with the supplied example

```sh
splice-girl example --output demo-inputs
splice-girl check --config demo-inputs/analysis.json
splice-girl run --config demo-inputs/analysis.json
```

The first command writes invented inputs, `analysis.json` and a walkthrough.
The check creates no output files. The run repeats preflight and saves an
offline `demo-inputs/configured-results/report.html` plus the original config,
resolved paths, checked hashes, eligibility decisions and sensitivity outputs.
The report labels the run as synthetic. Open it directly in a browser.

To repeat a run, change `output` in the configuration to a new folder. Existing
folders are protected. No raw-read processing or downloads are started.

## Configuration fields

| Field | Meaning |
|---|---|
| `schema_version` | Must be `1`. |
| `input_kind` | `synthetic` for invented inputs, or `pipeline` for verified completed pipeline exports. |
| `junctions`, `manifest`, `genes` | Paths to complete cohort counts, sample labels/metadata and the requested gene family. `genes` may be null for descriptive analysis only. |
| `output` | A new destination folder. |
| `complete_counts` | Must be true after reviewing exports; the pipeline mode additionally verifies counts against upstream records. |
| `shared_graph` | Explicit `minimum_edge_count`, `minimum_samples_per_edge`, `minimum_gene_support` integers. No hidden biological defaults. |
| `sensitivity` | Grids `edge_counts`, `sample_counts`, `gene_support`; set null to run only the primary setting. |
| `statistics` | Explicit `family`, `alternative`, `design`, `min_per_group`; null means descriptive analysis only. |
| `provenance` | Null for synthetic inputs; required metadata below for pipeline inputs. |

See `examples/shared_graph/analysis.json` for a complete example. The settings
there exercise a proposal; choose primary settings before seeing biological
case/control effects. The example's thresholds match the study values
approved on 2026-11-22 and confirmed final on 2026-11-24.

## Verify biological pipeline exports

Use the same top-level structure but set `input_kind` to `pipeline`, point
`junctions` to the completed batch's combined count table, and add:

```json
{
  "batch_manifest": "batch/batch_manifest.json",
  "reference_manifest": "reference/reference_manifest.json",
  "gene_map_manifest": "reference/junction_exon_map.tsv.manifest.json",
  "node_rule": "shared-splice-site",
  "strandedness": "RF",
  "minimum_anchor": 8,
  "exclude_duplicates": false
}
```

That object replaces `provenance: null`; it is not a standalone full config.
The eight-base anchor is an example setting, not an approved general default.
`genes` can point to your reviewed family; the rich analysis sample manifest
can add subject IDs, but its IDs/conditions must match the batch.

Preflight checks:

1. All intended samples have complete batch/sample records and unique labels.
2. Sample manifest hashes and individual count-file hashes match the batch
   records, and every per-sample row exactly matches the combined count table.
3. All eight recorded index-part hashes and the splice hints match the chosen
   completed reference manifest.
4. The map's GTF matches that reference's GTF, the map hash matches every
   sample's conversion input, and the node rule matches the configuration.
5. Every sample used the configured strand, anchor and duplicate policy,
   retained read-counter support 1, and completed its recorded stages. A
   sample whose run manifest records a strand check of `hold_for_review` is
   refused.
   Recorded counter/alignment arguments must agree with those method fields,
   and the recorded index prefix must identify the same eight index inputs.
6. Gene-family/schema/coverage/cycle rules and any independent-subject/sample-
   minimum requirements pass before outputs are created.

Reference preparation and standalone map building now write
`<map-path>.manifest.json`. Older maps without this record need rebuilding
into a new path from the original GTF and node rule. Existing historical
pilot results are not rewritten or silently relabeled.

Keep the batch's per-sample `run_manifest.json` and `junctions.tsv` files under
their sample folders. Local files are used when a batch folder is moved;
recorded paths remain provenance. Preflight reads these small records/tables,
not FASTQ, BAM or large genome/index files. It compares recorded reference
hashes, so it does not independently prove original read quality or authenticity.
Synthetic mode is a caller declaration and never claims biological verification.

All checked file hashes are saved in `preflight.json` and rechecked before and
after execution. A change fails the configured run; its root manifest remains
incomplete/failed. Treat nested completed outputs as partial evidence whenever
the root `run_manifest.json` is not `complete`.

`configuration.json` preserves the original file exactly.
`resolved_configuration.json` is runnable with the same schema and absolute
input/provenance paths. Copy it outside the saved result folder and change
the copy's `output` before rerunning, preserving the original result hashes.
The original relative-path configuration remains useful when moving the
complete input folder to another machine. Both require the corresponding data.

## Report and reproducibility

`report.html` links the primary descriptive/statistical reports, sensitivity
table, eligibility audits and saved configuration. Missing/no-cycle/low-coverage
genes stay visible. A run with no eligible genes produces an exclusion report
and skips inference. This does not approve the scientific cutoffs or establish
that Hodge measurements outperform simpler junction measurements.
