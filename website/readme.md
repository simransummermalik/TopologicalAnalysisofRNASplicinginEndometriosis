# Splice Girl website

A small educational website explaining Splice Girl to beginners. It uses
plain HTML, one CSS file and a small JavaScript file. No build step, package
installation, framework, downloaded font, account or backend is needed.

The site keeps the presentation's pink/plum colors and the complete original
ASCII artwork from `assets/splice_girl_banner.txt`. Topics are split into
separate pages. The demo shows one step at a time instead of putting all the
instructions into a long page; with JavaScript off, all steps remain readable.

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

- `index.html`: what the project is, its current status, team and original ASCII.
- `learn.html`: biology, graphs, normalization, Hodge theory, Rust and limits.
- `guide.html`: seven demo steps, commands, menu prompts and expected results.
- `results.html`: reading reports, output files, statistics context and common errors.
- `dictionary.html`: short definitions for biology, math, files and terminal terms.
- `styles.css`: shared colors, typography, page layout, mobile and print styles.
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
