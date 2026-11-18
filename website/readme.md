# Splice Girl website

Software documentation for Splice Girl, including a beginner quickstart,
concepts, output reference and glossary. It uses plain HTML, one CSS file
and a small JavaScript file. No build step, package installation, framework,
downloaded font, account or backend is needed.

The layout uses a compact dark header, a left navigation sidebar, system fonts
and a neutral reading column with restrained plum links. The homepage starts
with source installation, an example command and the input/output format.
The complete original `assets/splice_girl_banner.txt` artwork is preserved
inside the homepage's collapsible Terminal banner section. The quickstart
still shows one step at a time; with JavaScript off, all steps remain readable.

## Preview locally

From the repository root:

```sh
python3 -m http.server 8000 --bind 127.0.0.1 --directory website
```

Then open:

```text
http://127.0.0.1:8000/
```

Leave that terminal open while viewing the site. Press Control+C to stop the
server. If port 8000 is in use, use 8001 in both the command and address.
The server serves only `website/` and binds to this computer's loopback address.

You can also open `website/index.html` directly in a browser. Local-file
clipboard permissions vary; if Copy cannot access the clipboard, the page
selects the command and explains how to copy it manually. Copying does not
execute anything. Visitors run demo commands in their own terminal, from
the repository root, not from `website/`.

## Files

- `index.html`: installation, example usage, input/output formats, research status and terminal banner.
- `learn.html`: biology, graphs, normalization, Hodge theory, Rust and limits.
- `guide.html`: seven demo steps, commands, menu prompts and expected results.
- `results.html`: reading reports, output files, statistics context and common errors.
- `dictionary.html`: short definitions for biology, math, files and terminal terms.
- `styles.css`: header/sidebar layout, shared colors, system typography, mobile and print styles.
- `favicon.svg`: a small SG text icon for browser tabs.
- `script.js`: copy buttons, step navigation and browser-history handling.

Edit teaching text directly in its HTML file. Navigation is ordinary HTML
repeated on each page; update it in all five pages if you add or rename a page.
Colors live at the top of `styles.css`. No generated bundles need rebuilding.

## Content boundaries

The demo uses synthetic inputs and known answers from `examples/demo/` and
`tests/fixtures/`. Its group labels are not real patient results. The actual
one-control CDH11 prefix is a bounded technical check; the full eight-sample
biological comparison remains in progress. Cycle fraction is not a circular-RNA
percentage, and environmental evidence does not establish patient exposures.

The basic website guide does not start preprocessing, download reads, change
system settings, or run any commands in a visitor's browser. Links to the
repository and Rust setup are ordinary external links.

The website is local only. Nothing in this folder deploys or publishes it.

## Checked November 18, 2026

- All five pages load successfully, with working local links and section anchors.
- Desktop (1280 px) and phone (390 px) layouts have no horizontal page overflow.
- The original embedded ASCII matches `assets/splice_girl_banner.txt` exactly.
- Guide navigation works with keyboard Enter, Previous/Next, direct step links,
  browser Back, URL fragments and direct `file://` opening.
- With JavaScript disabled, all seven steps remain readable and inactive
  copy/navigation controls stay hidden. Native troubleshooting disclosures work.
- Copy success and denied-permission fallback were checked with controlled
  browser clipboard implementations; the operating system clipboard was not changed.
- The ten setup/analysis/menu/export command invocations copied from the guide
  passed in an isolated local source copy. Expected cycle fractions, both
  comparison effects and eight exported report SVGs were checked. The macOS
  `open` commands were checked for existing destinations, without launching
  the user's normal browser profile. Optional global installation was not rerun.

These are developer checks, not independent first-time reader feedback.
