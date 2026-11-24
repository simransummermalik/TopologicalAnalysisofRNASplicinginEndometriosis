#!/usr/bin/env python3
"""Render a full-slide animated code walkthrough of Splice Girl.

The code is read verbatim from src/graph.rs, src/hodge.rs and
src/permutation.rs at the line ranges below, with real line numbers; a dashed
row marks skipped lines. The view scrolls slowly from file to file, and at
each stop pink circles draw around the important lines one at a time, each
with a short plain-language note.

Requires matplotlib, NumPy and Pillow:

    python3 scripts/make_code_walkthrough.py --output assets/presentation/code_walkthrough.gif
"""

from __future__ import annotations

import argparse
import re
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from PIL import Image  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CODE_BG, TAB_BG, BAND = "#121314", "#1b1c1f", "#1d1f23"
MARK, PINK, LIGHT, GUTTER, RULE = "#ff6fae", "#e0becf", "#f6eff1", "#6e7681", "#30363d"
SYNTAX = {
    "keyword": "#ff7b72", "function": "#d2a8ff", "type": "#4ec9b0", "number": "#79c0ff",
    "string": "#a5d6ff", "comment": "#8b949e", "plain": "#d6dce3",
}
KEYWORDS = {"pub", "fn", "let", "mut", "for", "in", "if", "return", "as", "impl"}
PRIMITIVES = {"f64", "usize", "bool", "str"}
SANS, MONO = "Helvetica Neue", "Menlo"

# Layout in design pixels (1920 x 1080); fonts in points at 120 dpi.
W, H = 1920, 1080
TAB_H, TOP_PAD, LINE_H = 64, 22, 42
CODE_X, NUMBER_X, CODE_RIGHT = 128, 92, 1456
NOTE_X, NOTE_W_CHARS = 1504, 33
CODE_PT, NOTE_PT = 16, 13.5
VISIBLE = (H - TAB_H - TOP_PAD) // LINE_H

FILES = [
    ("src/graph.rs", [(31, 34), (38, 44), (56, 57), (77, 79)]),
    ("src/hodge.rs", [(47, 59), (78, 84), (92, 92)]),
    ("src/permutation.rs", [(59, 62), (79, 80), (84, 89), (96, 96)]),
]
# Each stop: file, then (first line, last line, note) for every circle.
STOPS = [
    ("src/graph.rs", [
        (32, 34, "Add up the reads. A gene with no reads is rejected."),
        (40, 41, "Sorted sets, so the same data always builds the same graph."),
        (79, 79, "Each junction becomes its share of the gene's reads."),
    ]),
    ("src/hodge.rs", [
        (48, 48, "Write the graph as a table: which junction leaves and enters each exon."),
        (49, 52, "One value per exon that best explains the flow: the straight-through (gradient) part."),
        (54, 59, "Whatever is left over is the loop (cycle) part."),
    ]),
    ("src/hodge.rs", [
        (78, 80, "Three checks: the parts add back, the loop balances, the parts don't overlap."),
        (82, 84, "Every gene gets a pass or a fail."),
        (92, 92, "Cycle fraction: the share of the flow that goes around loops."),
    ]),
    ("src/permutation.rs", [
        (79, 79, "The real difference between cases and controls."),
        (84, 89, "Try every way to relabel the samples: 70 ways for 4 cases vs 4 controls."),
        (96, 96, "p-value: how often a relabeling looks at least this different."),
    ]),
]
LEAD, MARK_GAP, DRAW, HOLD, FADE, SCROLL, END_HOLD = 1.0, 2.6, 1.0, 2.8, 0.5, 1.8, 4.0

TOKEN = re.compile(r"(//.*$)|(\"(?:[^\"\\]|\\.)*\")|(\b\d+(?:\.\d+)?\b)|([A-Za-z_][A-Za-z0-9_]*!?)|(\s+)|(.)")


def tokens(line):
    for match in TOKEN.finditer(line):
        comment, string, number, word, space, other = match.groups()
        if space:
            continue
        if comment:
            colour = SYNTAX["comment"]
        elif string:
            colour = SYNTAX["string"]
        elif number:
            colour = SYNTAX["number"]
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
        else:
            colour = SYNTAX["plain"]
        yield match.start(), match.group(), colour


def build_rows():
    rows, where = [], {}
    for index, (name, ranges) in enumerate(FILES):
        source = (ROOT / name).read_text().splitlines()
        if index:
            rows.append({"kind": "blank"})
        rows.append({"kind": "header", "file": name})
        for part, (first, last) in enumerate(ranges):
            if part:
                rows.append({"kind": "gap"})
            for number in range(first, last + 1):
                where[(name, number)] = len(rows)
                rows.append({"kind": "code", "file": name, "number": number,
                             "text": source[number - 1].rstrip()})
    return rows, where


def smooth(x):
    x = min(max(x, 0.0), 1.0)
    return x * x * (3 - 2 * x)


def ramp(t, start, end):
    return smooth((t - start) / (end - start))


def rounded_rect(x0, y0, x1, y1, r, n=10):
    corners = [((x1 - r, y0 + r), -90, 0), ((x1 - r, y1 - r), 0, 90),
               ((x0 + r, y1 - r), 90, 180), ((x0 + r, y0 + r), 180, 270)]
    points = [(x0 + r, y0)]
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


class Walkthrough:
    def __init__(self, ax, char_width):
        self.ax, self.cw = ax, char_width
        self.rows, self.where = build_rows()
        max_offset = max(0, len(self.rows) - VISIBLE + 1)
        self.timeline, t, seen = [], 1.0, set()
        for name, marks in STOPS:
            first_row = self.where[(name, marks[0][0])]
            target = first_row - 3
            if name not in seen:
                header = next(i for i, r in enumerate(self.rows) if r.get("file") == name and r["kind"] == "header")
                target = min(target, header - 1)
                seen.add(name)
            offset = min(max(target, 0), max_offset)
            starts = [t + LEAD + MARK_GAP * j for j in range(len(marks))]
            depart = starts[-1] + DRAW + HOLD
            self.timeline.append({"file": name, "marks": marks, "offset": offset,
                                  "arrive": t, "starts": starts, "depart": depart})
            t = depart + FADE + SCROLL
        self.duration = self.timeline[-1]["depart"] - HOLD + END_HOLD

    def offset(self, t):
        stops = self.timeline
        for here, after in zip(stops, stops[1:]):
            if t < here["depart"] + FADE:
                return here["offset"]
            if t < after["arrive"]:
                u = smooth((t - here["depart"] - FADE) / SCROLL)
                return here["offset"] + (after["offset"] - here["offset"]) * u
        return stops[-1]["offset"]

    def row_y(self, row, offset):
        return TAB_H + TOP_PAD + (row - offset + 0.5) * LINE_H

    def draw(self, t):
        ax = self.ax
        ax.clear()
        ax.set_xlim(0, W)
        ax.set_ylim(H, 0)
        ax.axis("off")
        offset = self.offset(t)
        reveal = ramp(t, 0.0, 1.0)
        current = max((s for s in self.timeline if t >= s["arrive"] - SCROLL), key=lambda s: s["arrive"],
                      default=self.timeline[0])

        # Circles, line bands and notes for every stop that has started.
        circled = set()
        for stop in self.timeline:
            notes = []
            for (first, last, note), start in zip(stop["marks"], stop["starts"]):
                draw = ramp(t, start, start + DRAW)
                if draw <= 0:
                    continue
                r0, r1 = self.where[(stop["file"], first)], self.where[(stop["file"], last)]
                texts = [self.rows[r]["text"] for r in range(r0, r1 + 1)]
                indent = min(len(s) - len(s.lstrip()) for s in texts if s.strip())
                width = max(len(s) for s in texts)
                x0 = CODE_X + indent * self.cw - 12
                x1 = CODE_X + width * self.cw + 12
                y0 = self.row_y(r0, offset) - LINE_H / 2 + 3
                y1 = self.row_y(r1, offset) + LINE_H / 2 - 3
                ax.add_patch(Rectangle((CODE_X - 8, y0), CODE_RIGHT - CODE_X + 16, y1 - y0,
                                       facecolor=BAND, edgecolor="none", alpha=draw, zorder=0))
                path = partial(rounded_rect(x0, y0, x1, y1, 9), draw)
                ax.plot(path[:, 0], path[:, 1], color=MARK, lw=2.2, solid_capstyle="round", zorder=4)
                circled.update(range(r0, r1 + 1))
                if stop is current:
                    alpha = ramp(t, start + 0.5, start + 1.2) * (1 - ramp(t, stop["depart"], stop["depart"] + FADE))
                    if stop is self.timeline[-1]:
                        alpha = ramp(t, start + 0.5, start + 1.2)
                    notes.append((alpha, x1, (y0 + y1) / 2, textwrap.wrap(note, NOTE_W_CHARS)))
            floor = TAB_H + 10
            for alpha, x1, centre, lines in notes:
                height = 30 * len(lines)
                top = max(centre - height / 2, floor)
                floor = top + height + 22
                if alpha <= 0.01:
                    continue
                ax.plot([x1, NOTE_X - 12], [centre, top + height / 2], color=MARK, lw=1.2,
                        alpha=0.8 * alpha, zorder=4)
                ax.add_patch(Rectangle((NOTE_X - 6, top), 4, height, facecolor=MARK,
                                       edgecolor="none", alpha=alpha, zorder=4))
                for k, text in enumerate(lines):
                    ax.text(NOTE_X + 10, top + 15 + 30 * k, text, fontfamily=SANS, fontsize=NOTE_PT,
                            color=LIGHT, va="center", alpha=alpha, zorder=5)

        # The code itself.
        first_row = int(np.floor(offset)) - 1
        for row in range(max(0, first_row), min(len(self.rows), first_row + VISIBLE + 3)):
            item, y = self.rows[row], self.row_y(row, offset)
            if item["kind"] == "header":
                ax.text(CODE_X, y + 2, item["file"], fontfamily=MONO, fontweight="bold",
                        fontsize=CODE_PT - 2, color=PINK, va="center", alpha=reveal)
                ax.plot([CODE_X + (len(item["file"]) + 1.5) * self.cw * 0.88, CODE_RIGHT], [y + 2, y + 2],
                        color=RULE, lw=1.2, alpha=reveal)
            elif item["kind"] == "gap":
                ax.plot([CODE_X, CODE_RIGHT], [y, y], color=RULE, lw=1.2, ls=(0, (4, 5)), alpha=reveal)
            elif item["kind"] == "code":
                ax.text(NUMBER_X, y, str(item["number"]), fontfamily=MONO, fontsize=CODE_PT - 3,
                        color=LIGHT if row in circled else GUTTER, ha="right", va="center", alpha=reveal)
                for column, text, colour in tokens(item["text"]):
                    ax.text(CODE_X + column * self.cw, y, text, fontfamily=MONO, fontsize=CODE_PT,
                            color=colour, va="center", ha="left", alpha=reveal, zorder=3)

        # Editor tab bar, drawn last so scrolled code slides underneath it.
        ax.add_patch(Rectangle((0, 0), W, TAB_H, facecolor=TAB_BG, edgecolor="none", zorder=8))
        x = 24
        for name, _ in FILES:
            label = name.split("/")[-1]
            width = len(label) * 12.5 + 48
            active = name == current["file"]
            if active:
                ax.add_patch(Rectangle((x, 8), width, TAB_H - 8, facecolor=CODE_BG, edgecolor="none", zorder=9))
                ax.add_patch(Rectangle((x, 8), width, 3, facecolor=MARK, edgecolor="none", zorder=10))
            ax.text(x + width / 2, TAB_H / 2 + 5, label, fontfamily=MONO, fontsize=13,
                    color=LIGHT if active else GUTTER, ha="center", va="center", zorder=10)
            x += width + 6
        ax.text(W - 28, TAB_H / 2 + 4, "Splice Girl  ·  inside the code", fontfamily=SANS, fontsize=14,
                color=PINK, ha="right", va="center", zorder=10)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, required=True, help="new .gif path")
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--fps", type=float, default=12.5,
                        help="GIF delays are whole 1/100 s, so use 10, 12.5, 20 or 25")
    parser.add_argument("--still", type=float, action="append", default=[],
                        help="render only this time in seconds to a PNG beside --output")
    args = parser.parse_args()

    figure = plt.figure(figsize=(16, 9), dpi=args.width / 16, facecolor=CODE_BG)
    ax = figure.add_axes([0, 0, 1, 1])
    probe = ax.text(0, 0, "M" * 50, fontfamily=MONO, fontsize=CODE_PT)
    figure.canvas.draw()
    char_width = probe.get_window_extent().width / 50 * W / figure.bbox.width
    walk = Walkthrough(ax, char_width)

    def render(t):
        walk.draw(t)
        figure.canvas.draw()
        return Image.fromarray(np.asarray(figure.canvas.buffer_rgba())[..., :3])

    if args.still:
        for t in args.still:
            path = args.output.with_name(f"{args.output.stem}_{t:05.2f}s.png")
            render(t).save(path)
            print(path)
        print(f"duration {walk.duration:.1f} s; stops at " +
              ", ".join(f"{s['arrive']:.1f}" for s in walk.timeline))
        return 0

    if args.output.exists():
        raise SystemExit(f"{args.output} already exists; choose a new path")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    samples = [render(s["depart"] - 0.1) for s in walk.timeline]
    strip = Image.new("RGB", (samples[0].width, samples[0].height * len(samples)))
    for i, image in enumerate(samples):
        strip.paste(image, (0, i * image.height))
    palette = strip.quantize(colors=160, method=Image.Quantize.MEDIANCUT)
    total = int(round(walk.duration * args.fps))

    def frames():
        for frame in range(total):
            yield render(frame / args.fps).quantize(palette=palette, dither=Image.Dither.NONE)

    stream = frames()
    next(stream).save(args.output, save_all=True, append_images=stream,
                      duration=round(1000 / args.fps), loop=0, optimize=False)
    print(f"Wrote {args.output} ({total} frames, {walk.duration:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
