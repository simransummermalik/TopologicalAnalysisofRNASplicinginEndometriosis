#!/usr/bin/env python3
"""Render the "inside the code" slide: four annotated excerpts from Splice Girl.

The excerpts are copied from src/graph.rs, src/hodge.rs and
src/permutation.rs. Some lines are dedented or joined, and `…` marks
omitted code; no code is otherwise changed. The image fills the area under
the deck's own heading and uses the deck's colours. Each card circles its
important lines with a short plain-language note.

Writes a static PNG and, optionally, a slow GIF in which the cards and
circles appear one at a time. Requires matplotlib, NumPy and Pillow:

    python3 scripts/make_code_slide.py \
        --png assets/presentation/code_overview.png \
        --gif assets/presentation/code_overview.gif
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle, FancyBboxPatch  # noqa: E402
from PIL import Image  # noqa: E402

BG, PLUM, ROSE, INKNOTE = "#eddfe3", "#6b4159", "#90456b", "#3f2534"
CODE_BG, MARK = "#121314", "#ff6fae"
SYNTAX = {
    "keyword": "#ff7b72", "function": "#d2a8ff", "type": "#4ec9b0", "number": "#79c0ff",
    "string": "#a5d6ff", "comment": "#8b949e", "plain": "#d6dce3",
}
KEYWORDS = {"pub", "fn", "let", "mut", "for", "in", "if", "return", "as", "impl"}
PRIMITIVES = {"f64", "usize", "bool", "str"}
SANS, MONO = "Helvetica Neue", "Menlo"

WIDTH, HEIGHT = 16.0, 6.4
MARGIN, GAP_X, GAP_Y = 0.8, 0.35, 0.16
CARD_W = (WIDTH - 2 * MARGIN - GAP_X) / 2
CARD_H = 3.02
LINE_H, CODE_PT = 0.192, 9.2
ROWS = 8
CARD_GAP_S = 6.0           # seconds between cards in the GIF
DURATION = 30.0

CARDS = [
    {
        "title": "Build the graph",
        "where": "graph.rs  ·  from_junctions",
        "lines": [
            "pub fn from_junctions(rows: &[Junction]) -> Result<Self, String> {",
            "    let total: f64 = rows.iter().map(|row| row.count).sum();",
            "    if !total.is_finite() || total <= 0.0 { return Err(…); }",
            "    let mut node_names = BTreeSet::new();",
            "    let mut edge_keys = BTreeSet::new();",
            "    …",
            "    let signal = edges.iter().map(|edge| edge.count / total).collect();",
        ],
        "marks": [
            (1, 2, "Add up the reads; reject a gene with none"),
            (3, 4, "Sorted sets: the same data always gives the same graph"),
            (6, 6, "Each junction becomes its share of the reads"),
        ],
    },
    {
        "title": "Split the flow",
        "where": "hodge.rs  ·  decompose",
        "lines": [
            "pub fn decompose(graph: &SpliceGraph) -> Result<Decomposition, String> {",
            "    let incidence = incidence_matrix(graph);",
            "    let potential = fit_potential(graph, &incidence)?;",
            "    let gradient = transpose_vector(&incidence, &potential);",
            "    let cycle: Vec<f64> = graph.signal.iter().zip(&gradient)",
            "        .map(|(signal, fitted)| signal - fitted)",
            "        .collect();",
        ],
        "marks": [
            (1, 1, "Write the graph as a table of exons and junctions"),
            (2, 3, "Best straight-through fit: the gradient part"),
            (4, 6, "Whatever is left over: the loop (cycle) part"),
        ],
    },
    {
        "title": "Check the math",
        "where": "hodge.rs  ·  decompose",
        "lines": [
            "    let reconstruction_error = norm(&reconstruction_delta);",
            "    let kernel_error = norm(&matrix_vector(&incidence, &cycle));",
            "    let orthogonality_error = dot(&gradient, &cycle).abs();",
            "    let passed = reconstruction_error <= TOLERANCE * scale",
            "        && kernel_error <= TOLERANCE * scale",
            "        && orthogonality_error <= TOLERANCE * scale * scale;",
            "    …",
            "        cycle_fraction: cycle_energy / signal_energy,",
        ],
        "marks": [
            (0, 2, "Parts add back, the loop balances, the parts don't overlap"),
            (3, 5, "Every gene gets a pass or fail"),
            (7, 7, "Cycle fraction: the number we compare"),
        ],
    },
    {
        "title": "Compare cases and controls",
        "where": "permutation.rs  ·  permutation_test",
        "lines": [
            "pub fn permutation_test(samples: &[IncludedSample], sidedness: Sidedness)",
            "    …",
            "    for case_subset in combinations(n, k) {",
            "        let effect = group_difference(&values, &case_subset, n);",
            "        if is_at_least_as_extreme(effect, observed_effect, sidedness, tolerance) {",
            "            extreme_count += 1;",
            "    …",
            "        p_value: extreme_count as f64 / permutations_considered as f64,",
        ],
        "marks": [
            (2, 5, "Try every way to relabel cases and controls (70 for 4 vs 4)"),
            (7, 7, "p-value: how often a relabeling looks at least this different"),
        ],
    },
]

TOKEN = re.compile(r"(//.*$)|(\"(?:[^\"\\]|\\.)*\")|(\b\d+(?:\.\d+)?\b)|([A-Za-z_][A-Za-z0-9_]*!?)|(\s+)|(.)")


def tokens(line):
    """Yield (column, text, colour) for a Rust source line."""
    for match in TOKEN.finditer(line):
        comment, string, number, word, space, other = match.groups()
        column = match.start()
        if space:
            continue
        if comment or other == "…":
            yield column, match.group(), SYNTAX["comment"]
        elif string:
            yield column, string, SYNTAX["string"]
        elif number:
            yield column, number, SYNTAX["number"]
        elif word:
            rest = line[match.end():].lstrip()
            if word in KEYWORDS:
                colour = SYNTAX["keyword"]
            elif rest.startswith("("):
                colour = SYNTAX["function"]
            elif word[0].isupper() or word in PRIMITIVES:
                colour = SYNTAX["type"]
            else:
                colour = SYNTAX["plain"]
            yield column, word, colour
        else:
            yield column, other, SYNTAX["plain"]


def smooth(x):
    x = min(max(x, 0.0), 1.0)
    return x * x * (3 - 2 * x)


def ramp(t, start, end):
    return smooth((t - start) / (end - start))


def rounded_rect(x0, y0, x1, y1, r, n=10):
    corners = [((x1 - r, y1 - r), 90, 0), ((x1 - r, y0 + r), 0, -90),
               ((x0 + r, y0 + r), -90, -180), ((x0 + r, y1 - r), 180, 90)]
    points = [(x0 + r, y1)]
    for (cx, cy), a0, a1 in corners:
        angles = np.radians(np.linspace(a0, a1, n))
        points += list(zip(cx + r * np.cos(angles), cy + r * np.sin(angles)))
    return np.array(points)


def partial(points, fraction):
    if fraction >= 1:
        return points
    steps = np.linalg.norm(np.diff(points, axis=0), axis=1)
    along = np.concatenate([[0.0], np.cumsum(steps)])
    stop = fraction * along[-1]
    keep = points[along <= stop]
    index = len(keep)
    if index < len(points):
        u = (stop - along[index - 1]) / max(steps[index - 1], 1e-9)
        keep = np.vstack([keep, points[index - 1] + u * (points[index] - points[index - 1])])
    return keep


class Slide:
    def __init__(self, ax, char_width):
        self.ax = ax
        self.cw = char_width

    def card(self, index, card, t, start):
        ax = self.ax
        alpha = ramp(t, start, start + 0.8)
        if alpha <= 0.01:
            return
        column, row = index % 2, index // 2
        x0 = MARGIN + column * (CARD_W + GAP_X)
        top = HEIGHT - 0.15 - row * (CARD_H + GAP_Y)

        ax.add_patch(Circle((x0 + 0.14, top - 0.17), 0.135, facecolor=PLUM, edgecolor="none", alpha=alpha))
        ax.text(x0 + 0.14, top - 0.172, str(index + 1), fontfamily=SANS, fontweight="bold", fontsize=11,
                color=BG, ha="center", va="center", alpha=alpha)
        ax.text(x0 + 0.38, top - 0.17, card["title"], fontfamily=SANS, fontweight="bold", fontsize=15,
                color=PLUM, va="center", alpha=alpha)
        ax.text(x0 + CARD_W, top - 0.17, card["where"], fontfamily=MONO, fontsize=9, color=ROSE,
                ha="right", va="center", alpha=alpha)

        panel_top = top - 0.4
        panel_h = 0.14 + ROWS * LINE_H + 0.1
        ax.add_patch(FancyBboxPatch((x0, panel_top - panel_h), CARD_W, panel_h,
                                    boxstyle="round,pad=0,rounding_size=0.07",
                                    facecolor=CODE_BG, edgecolor="none", alpha=alpha))
        text_x = x0 + 0.36

        def line_y(k):
            return panel_top - 0.14 - LINE_H * (k + 0.5)

        for k, line in enumerate(card["lines"]):
            for col, text, colour in tokens(line):
                ax.text(text_x + col * self.cw, line_y(k), text, fontfamily=MONO, fontsize=CODE_PT,
                        color=colour, va="center", ha="left", alpha=alpha)

        notes_top = panel_top - panel_h - 0.21
        for m, (first, last, note) in enumerate(card["marks"]):
            begin = start + 1.6 + 1.4 * m
            draw = ramp(t, begin, begin + 1.0)
            if draw <= 0.01:
                continue
            chosen = card["lines"][first:last + 1]
            indent = min(len(s) - len(s.lstrip()) for s in chosen)
            width = max(len(s.rstrip()) for s in chosen)
            rx0 = text_x + indent * self.cw - 0.07
            rx1 = text_x + width * self.cw + 0.07
            ry1 = line_y(first) + LINE_H / 2 - 0.012
            ry0 = line_y(last) - LINE_H / 2 + 0.012
            path = partial(rounded_rect(rx0, ry0, rx1, ry1, 0.05), draw)
            ax.plot(path[:, 0], path[:, 1], color=MARK, lw=1.7, solid_capstyle="round", zorder=5)
            letter = "abc"[m]
            shown = ramp(t, begin + 0.5, begin + 1.1)
            badge_y = (ry0 + ry1) / 2
            ax.add_patch(Circle((x0 + 0.18, badge_y), 0.085, facecolor=MARK, edgecolor="none",
                                alpha=shown, zorder=6))
            ax.text(x0 + 0.18, badge_y - 0.004, letter, fontfamily=SANS, fontweight="bold", fontsize=8.5,
                    color=CODE_BG, ha="center", va="center", alpha=shown, zorder=7)
            note_y = notes_top - 0.235 * m
            ax.add_patch(Circle((x0 + 0.09, note_y), 0.085, facecolor=MARK, edgecolor="none", alpha=shown))
            ax.text(x0 + 0.09, note_y - 0.004, letter, fontfamily=SANS, fontweight="bold", fontsize=8.5,
                    color=CODE_BG, ha="center", va="center", alpha=shown)
            ax.text(x0 + 0.27, note_y, note, fontfamily=SANS, fontsize=12, color=INKNOTE,
                    va="center", alpha=shown)

    def draw(self, t):
        ax = self.ax
        ax.clear()
        ax.set_xlim(0, WIDTH)
        ax.set_ylim(0, HEIGHT)
        ax.axis("off")
        for index, card in enumerate(CARDS):
            self.card(index, card, t, 0.4 + CARD_GAP_S * index)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--png", type=Path, required=True, help="new static .png path")
    parser.add_argument("--gif", type=Path, help="optional new animated .gif path")
    parser.add_argument("--gif-width", type=int, default=1600)
    parser.add_argument("--fps", type=float, default=10.0)
    args = parser.parse_args()
    for path in (args.png, args.gif):
        if path and path.exists():
            raise SystemExit(f"{path} already exists; choose a new path")

    figure = plt.figure(figsize=(WIDTH, HEIGHT), dpi=200, facecolor=BG)
    ax = figure.add_axes([0, 0, 1, 1])
    probe = ax.text(0, 0, "M" * 50, fontfamily=MONO, fontsize=CODE_PT)
    figure.canvas.draw()
    char_width = probe.get_window_extent().width / 50 / figure.dpi
    slide = Slide(ax, char_width)

    def render(t):
        slide.draw(t)
        figure.canvas.draw()
        return Image.fromarray(np.asarray(figure.canvas.buffer_rgba())[..., :3])

    args.png.parent.mkdir(parents=True, exist_ok=True)
    final = render(DURATION)
    final.save(args.png)
    print(f"Wrote {args.png} {final.size}")

    if args.gif:
        # Render GIF frames directly at GIF size; layout is in inches, so it is unchanged.
        figure.set_dpi(args.gif_width / WIDTH)
        palette = render(DURATION).quantize(colors=128, method=Image.Quantize.MEDIANCUT)
        total = int(round(DURATION * args.fps))

        def frames():
            for frame in range(total):
                yield render(frame / args.fps).quantize(palette=palette, dither=Image.Dither.NONE)

        stream = frames()
        next(stream).save(args.gif, save_all=True, append_images=stream,
                          duration=round(1000 / args.fps), loop=0, optimize=False)
        print(f"Wrote {args.gif} ({total} frames, {DURATION:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
