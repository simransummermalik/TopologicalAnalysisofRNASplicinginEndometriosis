# Beginner walkthrough results — November 17, 2026

**Result: the completed CLI walkthrough passed all 25 checks.** The separate
browser follow-up also passed, and its screenshot was visually inspected.
The walkthrough included `cargo test`: all **100 Rust tests passed**.
These are overlapping checks, not 125 separate automated unit tests.

Start with the [step-by-step guide](README.md). The saved
[comparison report](runs/2026-11-17-checked/artifacts/report/report.html),
[statistics report](runs/2026-11-17-checked/artifacts/statistics/report.html),
[graph](runs/2026-11-17-checked/artifacts/splice-girl-demo.svg), and
[browser screenshot](runs/2026-11-17-checked/artifacts/report-preview.png)
are available without rerunning the code.

## How the test was conducted

The reviewer was the coding assistant, acting through the documented beginner
steps. It already knows the project. No first-time human participant was
involved, so the independent beginner sign-off remains pending.

- System: macOS 26.6.2, ARM64.
- Rust: 1.97.1; Cargo: 1.97.1. Rust was already installed.
- Version tested: Splice Girl 0.1.0.
- Source base: commit `a58a2e1ce26d31d9ff96159c9403ef9461dd47d0` plus the
  recorded working-tree documentation changes.
- Setup: 288 tracked files copied into a new temporary source folder, without
  old build outputs. The [run manifest](runs/2026-11-17-checked/run_manifest.json)
  records the file hashes, source state, environment and command results.
- Installation: `cargo install --path .`, directed by an environment variable
  to an isolated temporary installation. The resolved executable was checked.
  The user's existing installed program was not replaced.
- Network: Cargo was set to offline mode. GitHub cloning, downloading Rust,
  and installing prerequisites on a new operating system were not tested.
- Inputs: supplied synthetic graph/junction tables. No FASTQ downloads,
  Cleaver processing, human-reference alignment or cohort inference ran.

The sum of the completed CLI command durations was about **22.1 seconds** on
this machine, including compilation. Source copying, reviewing instructions,
the earlier attempt and visual inspection are outside that number. This is
not an estimate of a beginner's reading time or another machine's runtime.

## Every command check

Each numbered log includes the exact command, working folder, menu input,
expected exit status, actual exit status, stdout and stderr. A deliberate
error counts as a successful check only when the expected error occurred.

| # | Check | Expected and observed result | Evidence |
|---|---|---|---|
| 1 | Rust version | Installed version printed; exit 0. | [Log](runs/2026-11-17-checked/logs/01_rust_version.txt) |
| 2 | Cargo version | Installed version printed; exit 0. | [Log](runs/2026-11-17-checked/logs/02_cargo_version.txt) |
| 3 | Current project folder | `Cargo.toml` and `README.md` both present. | [Log](runs/2026-11-17-checked/logs/03_project_folder.txt) |
| 4 | Fresh installation | Isolated installation succeeded; exit 0. | [Log](runs/2026-11-17-checked/logs/04_install.txt) |
| 5 | Installed command version | `splice-girl 0.1.0`. | [Log](runs/2026-11-17-checked/logs/05_installed_version.txt) |
| 6 | Help | Menu, graph and statistics commands described. | [Log](runs/2026-11-17-checked/logs/06_installed_help.txt) |
| 7 | Menu analysis | Banner appears; choice 1 analyzes the cycle; cycle fraction 1 and math checks pass; Return then `q` exits. | [Log](runs/2026-11-17-checked/logs/07_menu_cycle.txt) |
| 8 | Direct exon-skipping analysis | Cycle fraction 0.1111111111, math checks pass. | [Log](runs/2026-11-17-checked/logs/08_direct_exon_skip.txt) |
| 9 | Validate junction table | Validation passes; exit 0. | [Log](runs/2026-11-17-checked/logs/09_validate_junctions.txt) |
| 10 | Intentionally mismatched manifest | Names `orphan_01` and `control_03`; validation fails; expected exit 1. | [Log](runs/2026-11-17-checked/logs/10_intentional_manifest_mismatch.txt) |
| 11 | Compare both genes | Opposite effects of approximately +0.5556 and −0.5556. | [Log](runs/2026-11-17-checked/logs/11_compare_all_genes.txt) |
| 12 | Generate SVG | Graph file written and valid SVG/XML. | [Log](runs/2026-11-17-checked/logs/12_draw_graph.txt) |
| 13 | Repeat SVG destination | Existing figure protected; expected exit 1. | [Log](runs/2026-11-17-checked/logs/13_protect_graph.txt) |
| 14 | Save one gene | Tables and run manifest written; five sample graphs pass math checks. | [Log](runs/2026-11-17-checked/logs/14_save_gene.txt) |
| 15 | Repeat result destination | Existing result folder protected; expected exit 1. | [Log](runs/2026-11-17-checked/logs/15_protect_results.txt) |
| 16 | Save both-gene comparison | Result bundle includes gene summaries. | [Log](runs/2026-11-17-checked/logs/16_save_comparison.txt) |
| 17 | Rust tests | 59 library, 9 error, 11 statistics CLI, 9 workflow and 12 synthetic math tests pass: 100 total. | [Log](runs/2026-11-17-checked/logs/17_rust_tests.txt) |
| 18 | Generate report | HTML, tables and eight graph SVGs written. | [Log](runs/2026-11-17-checked/logs/18_browser_report.txt) |
| 19 | Optimized build | Release build succeeds. | [Log](runs/2026-11-17-checked/logs/19_release_build.txt) |
| 20 | Exact statistics | Both genes have 2 cases/2 controls, six assignments, raw and adjusted p-values of 1/3. | [Log](runs/2026-11-17-checked/logs/20_statistics.txt) |
| 21 | Missing file | Error names the requested missing path; expected exit 1. | [Log](runs/2026-11-17-checked/logs/21_missing_file.txt) |
| 22 | Zero support | Positive total support required; expected exit 1. | [Log](runs/2026-11-17-checked/logs/22_zero_support.txt) |
| 23 | Duplicate edge | Duplicate reported; expected exit 1. | [Log](runs/2026-11-17-checked/logs/23_duplicate_edge.txt) |
| 24 | Right relative path, wrong folder | File cannot be found there; expected exit 1. | [Log](runs/2026-11-17-checked/logs/24_wrong_folder.txt) |
| 25 | Recover from a bad menu path | Error appears; menu remains usable; `q` exits with status 0. | [Log](runs/2026-11-17-checked/logs/25_menu_error_recovery.txt) |

The checker also parsed the saved tables rather than relying only on the
terminal text. It confirmed both group counts and numeric effects; both
raw/adjusted p-values; six assignments per gene; valid SVG files; and working
relative links/anchors. After copying the artifacts into this folder, it
checked the links again: **17 local links and eight SVGs** for the descriptive
report, **26 local links and eight SVGs** for the statistics report.

## What was confusing, and what changed

| Finding | Change made |
|---|---|
| The guide assumed the terminal was already in the project. | Added `pwd`, a check for `Cargo.toml`/`README.md`, and a plain explanation of relative paths. |
| The menu instructions did not explain the output-folder question or how to return/quit. | Added the exact key sequence: 1, path, blank output, Return, `q`. |
| The installed command and `cargo run --` could look like different tools. | Explained installed versus source execution and what the `--` separator does. |
| The intentional failed manifest validation could look like a broken demo. | Explicitly labeled expected exit 1 and explained both mismatched samples. |
| The graph command did not give an explicit macOS opening command. | Added `open splice-girl-demo.svg`, browser instructions and repeat-run naming advice. |
| `compare` wording sounded as if the statistics feature had not been implemented. | Clarified that `compare` is descriptive and linked the separately configured `stats` workflow. |
| The complete run guide still said the human index and bounded real pilot were missing. | Updated it to distinguish completed technical work from pending full-cohort work. |

These were instruction/checker improvements. No application math or analysis
logic was changed during this walkthrough.

## First attempt and browser follow-up

The [first attempt](runs/2026-11-17/run_manifest.json) passed its 25 CLI checks
but recorded failure when headless Chrome aborted at the optional screenshot
step inside the execution sandbox. That original failed record is retained.
The checker also originally supplied the report folder rather than its
`report.html` file to the browser; that harness path was corrected. Neither
issue was a failed Hodge calculation or a broken `report` command.

The completed CLI run was recorded separately, without the optional browser
step. Chrome then rendered the copied `report.html` using a new temporary
profile outside the sandbox restriction. It returned exit 0, and its rendered
DOM contained both expected gene names and a report heading. The
[browser record](runs/2026-11-17-checked/browser_check.json),
[browser log](runs/2026-11-17-checked/logs/26_browser_followup.txt) and
[screenshot](runs/2026-11-17-checked/artifacts/report-preview.png) preserve
that follow-up. The screenshot shows the report, expected values and readable
academic formatting; it was visually inspected. This is automated browser
rendering, not an observation of a human clicking every control.

## What remains for another person

A teammate who did not write the tool should try the guide without live
coaching and record where they hesitate, misunderstand a term or cannot find
a file. That is the independent R17/U5 usability review. Windows/Linux behavior,
fresh Rust installation and slow/unavailable network installation also remain
outside this macOS walkthrough.

No final real-data analysis, disease association or exposure effect follows
from passing these synthetic checks. No commit or push was performed.
