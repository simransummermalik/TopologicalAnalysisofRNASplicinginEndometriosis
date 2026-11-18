# Technical pilot figures

Plain scientific exports for slides or Overleaf. Each figure is available as
PDF, SVG and PNG. All use saved tables; none show case/control differences.

| Figure | Caption and source |
|---|---|
| [Junction assignment](junction_assignment.pdf) | Assignment-stage outcomes for the first one million cleaned pairs of control SRR36710387. Bars count read-junction observations after upstream filters, not subjects or independent molecules. Most excluded observations identify one gene but several full exon pairs. Source: `../assignment_audit/summary.json`. |
| [Mate orientation](mate_orientation.pdf) | Orientation of primary unique spliced alignments passing the audit filters, counted once per alignment. R1 is approximately 99.72% antisense and R2 99.72% sense. This is consistent with RF, pending protocol confirmation. Source: the same summary's strand counts. |
| [CDH11 support sensitivity](cdh11_support_sensitivity.pdf) | Retained support and cycle fraction at four exploratory thresholds. The cycle-fraction axis spans 0–1 so numerical roundoff is not magnified. All four graphs have cycle rank zero. Source: `../verified/support_sensitivity.tsv`. |

The labels, counts, units and scope were checked against those source tables;
all three rendered PNGs were visually inspected. The
[figure manifest](figure_manifest.json) records input/script/output hashes.
For the main graph, use [verified/graph.svg](../verified/graph.svg), whose node
aliases have an explicit coordinate key. A final case/control figure remains
pending because there is only one biological sample in this pilot.

To regenerate, install Matplotlib in the plotting environment and run from
the repository root with a new output directory:

```sh
python3 scripts/plot_pilot_audit.py \
  --summary data/pilots/SRR36710387_prefix1m_nov17/assignment_audit/summary.json \
  --sensitivity data/pilots/SRR36710387_prefix1m_nov17/verified/support_sensitivity.tsv \
  --output-dir /tmp/splice-girl-pilot-figures
```

Matplotlib is optional for figure generation. It is not a Rust dependency.
