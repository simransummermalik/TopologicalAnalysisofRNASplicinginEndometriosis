# Splice Girl beginner walkthrough

This folder records an actual assistant-run attempt to use Splice Girl by
following the beginner instructions, from a fresh source copy. It also gives
you the same steps to follow yourself. The test uses the small **synthetic**
examples supplied with the project. No patient reads are downloaded or aligned.

**Date:** November 17, 2026.  
**Review type:** assistant-run walkthrough, with commands and evidence saved.
It is not a claim that someone new to the project independently reviewed it.
That separate team task remains open.

## What is in this folder?

| File or folder | What it is for |
|---|---|
| [RESULTS.md](RESULTS.md) | What passed, what was confusing, what was fixed, and the limits of this test. |
| [run_walkthrough.py](run_walkthrough.py) | Repeats the checks in a new temporary source copy and saves its findings. |
| [runs/2026-11-17](runs/2026-11-17/) | First attempt. The CLI checks passed; optional browser startup failed in the restricted execution environment. Retained so the first attempt is visible. |
| [runs/2026-11-17-checked](runs/2026-11-17-checked/) | Completed CLI walkthrough after documentation/checker corrections; one log per command and copies of the results. |

The detailed guide below is for running the tool yourself. All paths in the
commands are relative to the **project root**, the folder containing
`Cargo.toml` and the main `README.md`. Do not run them from `beginner_test/`.

## Before you start

You need the downloaded project and Rust 1.85 or newer. Cargo comes with Rust.
Python is only needed for the optional walkthrough runner at the end; it is
not needed for the normal Rust demo. Rust was already installed on the test
machine, so downloading/installing Rust itself was not tested.

Some words used below:

- **Terminal:** the app where you type commands.
- **Command:** one instruction you type, followed by Return/Enter.
- **Current folder:** where the terminal looks for a relative file path.
- **TSV:** a text table whose columns are separated by tabs.
- **CLI:** a program you use by typing commands.
- **Exit status:** `0` means a command succeeded; a nonzero number means it
  reported a problem. Some examples below deliberately test a problem.

You do not type the terminal's prompt, such as `(base) summer@Mac %`.
Only copy the commands inside the code blocks. In a command split across
lines, the backslash `\` means the command continues on the next line.

## 1. Make sure the terminal is in the project

Open the project folder in your terminal, then run:

```sh
pwd
ls Cargo.toml README.md
```

`pwd` prints the current folder's full path. The second command should list
both files. If it says they do not exist, you are in the wrong folder.
Use `cd` followed by the actual project path on your computer; put a path in
quotes if it has spaces. For example, `cd "/actual/path/to/project"` uses a
placeholder path that you must replace.

**Behind the scenes:** no analysis happens. This checks that commands like
`tests/fixtures/cycle.tsv` will point to the intended files.

## 2. Check that Rust and Cargo are available

```sh
rustc --version
cargo --version
```

Both commands should print a version. `rustc` must be at least 1.85.
If the terminal says `command not found`, use the Rust installation link in
the [main guide](../README.md#1-check-rust), then reopen the terminal.

**Behind the scenes:** `rustc` is the compiler; Cargo builds the project and
manages Rust commands. These commands only report the installed versions.

## 3. Install the Splice Girl command

From the project root:

```sh
cargo install --path .
splice-girl --version
splice-girl --help
```

The dot means "use the Rust project in this current folder." Cargo builds
the program, then installs its executable. The version should include
`splice-girl 0.1.0`; help should list commands such as `analyze`, `validate`,
`compare`, `report`, `stats` and `visualize`.

Cargo normally installs commands in `$HOME/.cargo/bin`. If installation
succeeds but the shell cannot find `splice-girl`, that folder may be missing
from `PATH`, the list of folders your shell searches for programs. You can
still use `cargo run --` from this project without changing `PATH`.

**What the test did:** it ran the same install command with Cargo's installation
root set to a temporary folder. It then verified that the tested executable
came from that folder. It did not replace your account's existing installation.

## 4. Know the two ways to run it

These both run the same analysis:

```sh
splice-girl analyze tests/fixtures/cycle.tsv
cargo run -- analyze tests/fixtures/cycle.tsv
```

`splice-girl` runs the installed executable. `cargo run --` builds/runs the
source in this checkout; everything after `--` is passed to Splice Girl.
The source form is useful after editing code because an older installed
executable does not automatically update when source files change.

The remaining steps use the source form for consistency with the README.

## 5. Open the menu and analyze a simple cycle

```sh
cargo run -- --start
```

Follow each prompt in order:

1. At **Choose an option**, type `1` and press Return.
2. At **Junction TSV path**, type `tests/fixtures/cycle.tsv` and press Return.
3. At **Output folder**, press Return without typing anything. This means
   print the answer in the terminal without saving a result folder.
4. Read the answer. Expect `cycle_fraction: 1.0000000000` and `checks: PASS`.
5. At **Press Enter to return to the menu**, press Return.
6. Type `q` and press Return to leave the menu.

**Data:** this is a made-up three-node directed cycle with equal junction
counts. It is built to have a known answer.

**Behind the scenes:** the code reads the rows, groups them by sample/gene,
builds the graph, normalizes counts, runs the Hodge decomposition and checks
its mathematical identities. The cycle fraction is one for this example.
This means all of this graph's signal lies in its cycle space; it is not a
measurement of the fraction of circular RNA in a patient.

## 6. Run an example directly

```sh
cargo run -- tests/fixtures/exon_skip.tsv
```

Expect `cycle_fraction: 0.1111111111` and `checks: PASS`.

**Data:** a synthetic exon-skipping shape: E1→E2, E2→E3 and E1→E3, with
equal support. Its mathematical answer is `1/9`, about 11.1%.

**Behind the scenes:** the program runs the same analysis without menu
prompts. The alternative path creates an undirected graph loop even though
the arrows do not lead back to E1. A cycle-space component does not require
a directed biological return path.

## 7. Validate the demo input

```sh
cargo run -- validate examples/demo/junctions.tsv
```

Expect `validation_status: PASS`.

**Behind the scenes:** validation uses the actual parser, graph builder and
mathematical checks. It is more than a check that the file exists. It does
not establish that a synthetic file is real biological evidence.

## 8. See an intentional metadata error

```sh
cargo run -- validate examples/demo/junctions.tsv examples/demo/manifest.tsv
```

This time, **expect failure**. The output should name `orphan_01` and
`control_03`, show `validation_status: FAIL`, and return exit status `1`.

`manifest.tsv` is the small sample-label table. It deliberately lists
`control_03` without junction data. `junctions.tsv` deliberately has an
`orphan_01` sample that is missing from that label table. The error means
the program noticed both mismatches. Do not edit the files to hide them.

After a command in a macOS/Linux shell, `echo $?` shows its exit status.
For this deliberate example it should print `1`; for step 7 it should print
`0`. This value refers only to the immediately preceding command.

## 9. Compare the two synthetic genes

```sh
cargo run -- compare examples/demo/all_genes.tsv examples/demo/manifest.tsv
```

Expect these case-minus-control cycle-fraction differences:

| Gene | Expected difference |
|---|---:|
| `GENE_HIGH_CYCLE` | `+0.5555555556` |
| `GENE_LOW_CYCLE` | `-0.5555555556` |

**Data:** `all_genes.tsv` has 22 junction rows, two made-up genes and four
measured samples: two cases and two controls. The fifth manifest sample,
`control_03`, has no measurements and is excluded explicitly.

**Behind the scenes:** the program analyzes each sample/gene graph and
compares group summaries. `compare` produces descriptive differences, not
p-values. These synthetic labels and values are not disease results.

## 10. Make a graph you can open

```sh
cargo run -- visualize case_02 GENE_DEMO examples/demo/junctions.tsv splice-girl-demo.svg
```

The file should appear in the project root. On macOS:

```sh
open splice-girl-demo.svg
```

On another system, open the file in a browser using its file manager or
the browser's Open File command. No web server is needed.

The figure should show the observed, gradient and cycle signals using the
same graph layout, with a table of junction values. `case_02` and `GENE_DEMO`
choose the sample/gene in the input; the figure is calculated from that data.

If you repeat the generation command with the same destination, expect an
`already exists` error. Use a new name such as `splice-girl-demo-2.svg`.
The walkthrough checked that existing figures are protected.

## 11. Save one gene's results

```sh
cargo run -- gene GENE_DEMO examples/demo/junctions.tsv examples/demo/manifest.tsv \
  --output demo_results
```

Look for a new `demo_results/` folder containing:

- `sample_metrics.tsv`: each analyzed sample's graph and Hodge measurements;
- `junction_components.tsv`: each junction's normalized, gradient and cycle values;
- `comparisons.tsv`: included/excluded samples and the reasons;
- `run_manifest.json`: the command, input hashes and settings.

The input has five mathematically analyzable sample graphs, including
`orphan_01`. Only four have both junctions and usable case/control metadata
for the group comparison. The exclusions remain visible.

Repeat this command with the same `demo_results` destination and expect an
`already exists` error. That protection prevents an old result being silently
replaced. Pick a different folder for a new run.

## 12. Save the comparison for both genes

```sh
cargo run -- compare examples/demo/all_genes.tsv examples/demo/manifest.tsv \
  --output demo_all_genes
```

This produces the comparison bundle plus `gene_summaries.tsv`, where you can
check the two effects from step 9. A TSV opens as text or as a spreadsheet;
use tabs as the separator if your spreadsheet asks.

## 13. Run the Rust tests

```sh
cargo test
```

Expect **100 passing tests** and no failed tests for this revision. Cargo
prints several groups; add their counts rather than expecting one line
saying 100. The zero-test groups are normal executable/doc-test groups.

The tests cover known mathematical answers, invalid inputs, saved outputs,
reports and statistics. This command does not require the large human data,
Cleaver, HISAT2, a Galaxy account or Python packages.

## 14. Open the full offline report

```sh
cargo run -- report examples/demo/all_genes.tsv examples/demo/manifest.tsv \
  --output target/demo-report
```

On macOS:

```sh
open target/demo-report/report.html
```

Otherwise open `report.html` in a browser. It should show both genes, their
sample measurements, group differences and exclusions. Click **View
decomposition** to open a sample's graph. There should be eight graph files:
four measured samples multiplied by two genes.

**Behind the scenes:** `report` performs the descriptive comparison and
writes tables, HTML and generated graph SVGs. It does not upload data or
start a website. Keep the entire report folder together when sharing it,
because its HTML links to files inside that folder.

## 15. Optional: try the exact statistical example

This is a longer command, so it is a separate step after the basic demo:

```sh
cargo build --release
cargo run --release -- stats examples/demo/all_genes.tsv \
  examples/demo/manifest.tsv --genes examples/demo/stats_genes.tsv \
  --family synthetic-demo --alternative two-sided --design independent \
  --min-per-group 2 --rank-by absolute-effect --output target/demo_complete
```

Open `target/demo_complete/report.html`. Both genes should have raw and
BH-adjusted p-values of `0.3333333333`, approximately one third.

| Part of the command | Plain-language meaning |
|---|---|
| `--release` | Use the optimized build. |
| `stats` | Run the explicitly configured statistical workflow. |
| `--genes` | Give the fixed gene list to test. |
| `--family synthetic-demo` | Name this testing family in the saved results. |
| `--alternative two-sided` | Look for a difference in either direction. |
| `--design independent` | Declare that the samples are independent, not paired. |
| `--min-per-group 2` | Require two eligible samples in each group for this toy example. |
| `--rank-by absolute-effect` | Sort by effect magnitude; sorting is not a significance claim. |
| `--output` | Choose a new destination for the results and report. |

There are only six ways to relabel this toy example's two cases and two
controls. The walkthrough checked that the saved table reports six
assignments per gene and the expected p-values. The actual study's frozen
policy uses four per group; the toy example does not change that policy.

## 16. Check that mistakes get useful responses

These were also tested. They should report errors and exit with status `1`:

```sh
splice-girl analyze does-not-exist.tsv
splice-girl analyze tests/fixtures/zero_support.tsv
splice-girl analyze tests/fixtures/duplicate_edge.tsv
```

The first should name the missing file. The second should explain that
positive total support is required. The third should identify a duplicate
edge. None should silently invent a result.

The walkthrough also tested a correct-looking relative path from the wrong
folder. The program cannot find it there; returning to the project root fixes
the location problem. In the menu, a bad path should show an error and let
you try again or type `q`, rather than closing the whole program.

## Repeat the recorded walkthrough automatically

From the project root, choose an output folder that does not already exist:

```sh
python3 beginner_test/run_walkthrough.py \
  --output-dir beginner_test/runs/my-new-run
```

This runner needs Git, Python and Rust already installed. It copies tracked
working-tree files to a new temporary folder, without old `target/` builds,
and runs the commands there. A source manifest identifies that copy. It
does not reset, stage, commit or push your repository. It sets an isolated
Cargo install destination and uses offline mode; it does not test a fresh
Rust installation or a download from GitHub.

Afterward, open `run_manifest.json` in the output folder. Each entry names
the command, supplied menu input, expected and actual exit code, elapsed time
and full log. `artifacts/` keeps the demo graph, reports, tables and example
inputs. Large temporary compiled files stay outside the repository.

The optional `--chrome /actual/path/to/chrome` argument requests a headless
browser screenshot too. It is not required to test the CLI. A restricted
execution environment may block browser startup even when the CLI and saved
HTML are working; the first attempt documents that case.

## What this does and does not establish

This checks that the documented commands work in a fresh local source copy,
that expected synthetic answers agree, and that common mistakes are reported.
The saved records distinguish expected errors from unexpected failures.
The reviewer already knows the project and the machine already has Rust.
It does not replace a first-time person's feedback, test every operating
system, validate all eight biological samples or establish a disease finding.
