# Presentation notes

This file contains a shared speaking plan, plain-language explanations and
demo commands. The four sections below are available to any of the four team
members; they are not assigned roles. Target time is about five minutes,
excluding questions. Rehearsal and an independent audience check are pending.

## 1. Question and gene evidence — about 50 seconds

> "RNA splicing joins pieces of a gene in different ways. We are building
> Splice Girl to compare those connection patterns in endometriosis and
> control samples. Our candidate list combines published genetic-risk
> evidence with genes reported to respond to dioxin-like chemicals. We have
> 661 candidate genes, including seven in both lists. These labels help us
> choose what to examine; they do not show that these patients were exposed."

Handoff: "Before we can compare those patterns, we need reliable junction
counts from the sequencing reads."

## 2. Reads and the completed pilot — about 55 seconds

> "We have identified four case and four control samples from the same tissue
> and menstrual phase. The tools first inspect and clean the reads, align
> them to the human reference, and count supported exon connections. We have
> completed a bounded technical run using one million cleaned read pairs
> from one control. It aligned successfully and produced junction tables.
> We also checked why many observations were excluded: usually the gene was
> clear, but several annotated exon pairs shared the same splice boundary.
> We are reviewing that before processing the full cohort."

Handoff: "Once the junction table is checked, the Rust program turns it into
something we can measure and display."

## 3. Graph, Hodge theory and live commands — about two minutes

> "Exons become graph nodes, and supported splice junctions become arrows.
> We divide each junction's count by the total for that sample and gene so
> the graph describes relative usage. Hodge decomposition separates the
> signal into a part explained by values on the nodes and a remaining
> cycle-space part. The cycle fraction measures how much squared signal is
> in that remaining part. It is a mathematical summary, not a percentage of
> circular RNA or a direct measure of disease."

Run these from the repository root. Build beforehand; every output folder
below is new, so the tool never needs to overwrite an old result.

```sh
cargo build --release
demo_dir="$(mktemp -d /tmp/splice-girl-demo.XXXXXX)"
target/release/splice-girl validate examples/demo/all_genes.tsv
target/release/splice-girl report examples/demo/all_genes.tsv examples/demo/manifest.tsv --output "$demo_dir/synthetic"
open "$demo_dir/synthetic/report.html"
```

Say: "This first example is synthetic. It is deliberately small so we can
check the answers. The program reads the table, groups rows by sample and
gene, computes the decomposition, checks the math identities, and writes
tables plus a report. Each graph is generated from that run's input."

The standalone validation above checks junction structure and math. The demo
manifest intentionally includes `control_03` without junction data so the
report can demonstrate an exclusion. Adding that manifest to `validate` is an
intentional mismatch check and exits with failure; that is not the command
used for the successful opening demonstration.

Then run the already checked real-gene input:

```sh
target/release/splice-girl analyze data/pilots/SRR36710387_prefix1m_nov17/verified/junctions.tsv --output "$demo_dir/real"
open data/pilots/SRR36710387_prefix1m_nov17/verified/report.html
```

Say: "This is CDH11 from the real control prefix. We selected it with a
written coverage and graph-size rule. An independent BAM recount matches its
counts, and an independent calculation agrees with Rust. Its retained graph
has no independent loop, so its cycle fraction is zero. That is the expected
answer for this graph; it does not mean the patient has no alternative
splicing. There is no disease comparison in this one-sample result."

The saved HTML is the original independently verified pilot report; the new
`$demo_dir/real` folder is the live recomputation. The selected input already
has the pilot support-three filter applied: four edges, six nodes, 103 support.

Handoff: "The next step is to apply the reviewed rules consistently to all
eight samples and compare the groups."

## 4. Statistics, limits and next steps — about 55 seconds

> "The software can compare each gene's average cycle fraction between the
> groups, export simpler baselines, and repeat the analysis across declared
> coverage settings. Our planned primary test summarizes candidate-gene
> differences and checks all 70 ways to relabel four cases and four controls,
> under the independent-group assumptions. The full biological result is
> still pending. We need the remaining reads, reviewed preprocessing settings
> and comparable graphs. This is a small study, and a graph cycle signal
> alone cannot establish circular RNA, disease causation or exposure effects."

Closing line: "We now have a tested tool and a checked real-input example;
the remaining work is to establish what the complete cohort actually supports."

## Offline backup and likely questions

- Open `examples/demo/README.md` for expected synthetic values if a live
  command fails. The installed release binary and local inputs need no internet.
- Keep the saved real pilot HTML, graph and three audit PNGs open before the
  talk. A recording and a timed rehearsal still need to be made by the team.
- **Why Rust?** It gives us a standalone compiled command with explicit input
  errors, reusable analysis code and automated checks. It does not improve
  biological validity by itself.
- **Why Hodge?** It summarizes the relationships among junctions, beyond one
  total read count. Baselines will test whether that extra summary is useful.
- **What is a cycle here?** An undirected graph loop or alternative path can
  create cycle space even when every biological arrow moves forward.
- **Are the results significant?** No disease significance is claimed from
  the one-control pilot. The shipped synthetic demo checks calculations.
- **Is everything finished?** The core tool and bounded technical validation
  work. Full-cohort processing, biological comparison and independent review remain.

## Interpretation boundary: environmental evidence (E11)

Source: the interpretation boundary in
[`gene_set_methods.md`](gene_set_methods.md#interpretation-boundary-e11).

**Line to say when the environmental gene set is introduced:**

> "We used published dioxin studies to choose which genes to look at. Our
> eight patients have no exposure measurements, so we can't say whether any
> of them were exposed, or that exposure caused anything we see. And most of
> that dioxin evidence comes from liver and breast-cancer cells, not
> endometrium."

**Points the speakers should have ready:**

- The environmental set is a way to *prioritize* genes. It is not a finding
  about these patients.
- GSE315857 records only diagnosis, tissue, and cycle phase. There are no
  dioxin levels, questionnaires, or exposure histories.
- An "overlap gene" has both genetic-risk and dioxin-response evidence. That
  is not a demonstrated gene–environment interaction.
- Only 21 genes have dioxin evidence from endometrial or uterine studies. For
  a specific gene, the evidence table shows which kind of evidence it has.

**Phrases to avoid on slides, in captions, and when answering questions:**

- "caused by exposure"
- "exposure-driven"
- "dioxin-affected patients"
- "gene–environment interaction"
- "dioxin changes splicing in endometriosis"

If asked "Does dioxin cause these splicing differences?": "We can't test
that with this data. There are no exposure measurements for these women. The
dioxin evidence only told us which genes were worth examining."
