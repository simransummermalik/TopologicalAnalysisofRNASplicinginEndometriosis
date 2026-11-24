#!/usr/bin/env python3
"""Render the animated "What we built" slide: why Splice Girl's analysis is fast.

Three slow scenes, one line of text each. Every number comes from the
repository's own records:

- one case sample, SRR36710386 (logs/november_19th_summer.md): alignment plus
  junction extraction 1 h 41 min 47 s (Slurm job 27844959); Rust Hodge
  analysis 5 s (Slurm job 27866517);
- 25 passing CLI walkthrough checks (docs/software_validation.md, November 17);
- 125 passing Rust tests (`cargo test` on commit 7594e06, November 24).

The small graphs in the first scene are an illustration, not real genes.
The colours match the presentation deck. Requires matplotlib, NumPy, Pillow:

    python3 scripts/make_speed_animation.py --output assets/presentation/what_we_built.gif
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from PIL import Image  # noqa: E402

BG, PINK, LIGHT, DIM = "#49283f", "#e0becf", "#f6eff1", "#7a5a6f"
WIDTH, HEIGHT = 16.0, 9.0
DURATION = 40.0
LEFT, RIGHT = 1.4, 14.6

ALIGN_SECONDS, HODGE_SECONDS = 1 * 3600 + 41 * 60 + 47, 5
TESTS, WALKTHROUGH = 125, 25

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["cmr10"],
        "mathtext.fontset": "cm",
        "axes.unicode_minus": False,
        "axes.formatter.use_mathtext": True,
    }
)
warnings.filterwarnings("ignore", message=".*cmr10.*")


def smooth(x):
    x = min(max(x, 0.0), 1.0)
    return x * x * (3 - 2 * x)


def ease_out(x):
    x = min(max(x, 0.0), 1.0)
    return 1 - (1 - x) ** 3


def ramp(t, start, end):
    return smooth((t - start) / (end - start))


def window(t, start, end, fade_in=1.0, fade_out=1.0):
    return ramp(t, start, start + fade_in) * (1 - ramp(t, end - fade_out, end))


def rgb(color):
    color = color.lstrip("#")
    return np.array([int(color[i:i + 2], 16) / 255 for i in (0, 2, 4)])


def mix(a, b, u):
    a = rgb(a) if isinstance(a, str) else np.asarray(a)
    b = rgb(b) if isinstance(b, str) else np.asarray(b)
    return tuple(a + (b - a) * min(max(u, 0.0), 1.0))


def make_mini_graphs(columns=10, rows=4, seed=7):
    """Small made-up splice graphs, used only as an illustration."""
    rng = np.random.default_rng(seed)
    graphs = []
    cell_w, cell_h = (RIGHT - LEFT) / columns, 0.98
    for row in range(rows):
        for column in range(columns):
            n = int(rng.integers(3, 7))
            x0 = LEFT + column * cell_w + 0.14 * cell_w
            edges = [(i, i + 1) for i in range(n - 1)]
            for _ in range(int(rng.integers(0, 3))):
                i = int(rng.integers(0, n - 2))
                j = int(rng.integers(i + 2, n))
                if (i, j) not in edges:
                    edges.append((i, j))
            graphs.append({
                "xs": x0 + np.linspace(0, 0.72 * cell_w, n),
                "y": 4.75 - row * cell_h,
                "edges": edges,
                "appear": 1.8 + 2.6 * rng.random(),
                "solve": 5.6 + 0.42 * column + 0.12 * row,
            })
    return graphs


class Scene:
    def __init__(self):
        self.graphs = make_mini_graphs()

    @staticmethod
    def headline(ax, text, alpha):
        if alpha > 0.01:
            ax.text(LEFT, 6.35, text, fontsize=30, color=LIGHT, alpha=alpha, va="center")

    def graphs_scene(self, ax, t):
        fade = window(t, 0.8, 14.0)
        if fade <= 0.01:
            return
        self.headline(ax, "Each gene becomes its own small graph", fade)
        for g in self.graphs:
            appear = ramp(t, g["appear"], g["appear"] + 0.8) * fade
            if appear <= 0.01:
                continue
            flash = max(0.0, 1 - abs(t - g["solve"] - 0.3) / 0.5)
            color = mix(mix(DIM, PINK, ramp(t, g["solve"], g["solve"] + 0.6)), LIGHT, 0.8 * flash)
            xs, y = g["xs"], g["y"]
            for i, j in g["edges"]:
                if j == i + 1:
                    ax.plot([xs[i], xs[j]], [y, y], color=color, lw=1.8, alpha=appear,
                            solid_capstyle="butt")
                else:
                    s = np.linspace(0, 1, 30)
                    height = 0.16 + 0.09 * (j - i)
                    ax.plot(xs[i] + (xs[j] - xs[i]) * s, y + height * 4 * s * (1 - s),
                            color=color, lw=1.5, alpha=appear)
            ax.scatter(xs, [y] * len(xs), marker="s", s=60, color=color, alpha=appear,
                       edgecolors="none", zorder=3)
        ax.text(RIGHT, 1.3, "(illustration)", fontsize=12, color=DIM, alpha=fade,
                ha="right", va="center")

    def race(self, ax, t):
        fade = window(t, 14.0, 27.0)
        if fade <= 0.01:
            return
        self.headline(ax, "The slow part is not the math", fade)
        full = RIGHT - LEFT
        bars = [
            (4.3, 15.6, 20.6, ALIGN_SECONDS, "Alignment", DIM),
            (2.3, 21.8, 22.1, HODGE_SECONDS, "Splice Girl", PINK),
        ]
        for y, start, end, seconds, label, color in bars:
            shown = ramp(t, start - 0.8, start) * fade
            if shown <= 0.01:
                continue
            ax.text(LEFT, y + 0.7, label, fontsize=24, color=LIGHT, alpha=shown, va="center")
            grow = smooth((t - start) / (end - start))
            width = max(0.07, full * seconds / ALIGN_SECONDS) * grow
            if width > 0:
                ax.add_patch(Rectangle((LEFT, y - 0.28), width, 0.56, facecolor=color,
                                       edgecolor="none", alpha=shown))
            if seconds >= 60:
                minutes = int(round(seconds * grow / 60))
                text = f"{minutes} min" if minutes < 60 else f"{minutes // 60} h {minutes % 60:02d} min"
                ax.text(LEFT + width - 0.18, y, text, fontsize=24, color=LIGHT,
                        alpha=shown * min(1, 3 * grow), ha="right", va="center")
            else:
                ax.text(LEFT + 0.4, y, "5 seconds", fontsize=40, color=PINK,
                        alpha=shown * ramp(t, end, end + 0.8), va="center")

    def stats(self, ax, t):
        fade = ramp(t, 27.0, 28.0)
        if fade <= 0.01:
            return
        self.headline(ax, "And it checks its own work", fade)
        width = (RIGHT - LEFT - 2 * 0.6) / 3
        columns = [
            (28.2, lambda u: f"{int(round(TESTS * u))}", "tests passing"),
            (29.8, lambda u: f"{int(round(WALKTHROUGH * u))}/{WALKTHROUGH}", "walkthrough checks"),
            (31.4, lambda u: f"{int(round(HODGE_SECONDS * u))} s", "for a whole sample"),
        ]
        for i, (start, number, label) in enumerate(columns):
            x = LEFT + i * (width + 0.6)
            shown = ramp(t, start, start + 0.8)
            if shown <= 0.01:
                continue
            ax.plot([x, x + ease_out((t - start) / 1.2) * width], [4.75, 4.75], color=PINK,
                    lw=1.6, alpha=shown)
            ax.text(x, 3.85, number(ease_out((t - start) / 2.5)), fontsize=60, color=PINK,
                    alpha=shown, va="center")
            ax.text(x, 2.75, label, fontsize=24, color=LIGHT, alpha=shown, va="center")

    def draw(self, ax, t):
        ax.clear()
        ax.set_xlim(0, WIDTH)
        ax.set_ylim(0, HEIGHT)
        ax.set_aspect("equal")
        ax.axis("off")
        header = ramp(t, 0.0, 0.8)
        ax.text(LEFT, 7.95, "PROGRESS", fontsize=15, color=PINK, alpha=header, va="center")
        ax.text(LEFT, 7.3, "What we built", fontsize=36, color=LIGHT, alpha=header, va="center")
        self.graphs_scene(ax, t)
        self.race(ax, t)
        self.stats(ax, t)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, required=True, help="new .gif path")
    parser.add_argument("--width", type=int, default=960)
    parser.add_argument("--fps", type=float, default=12.5,
                        help="GIF delays are whole 1/100 s, so use 10, 12.5, 20 or 25")
    parser.add_argument("--still", type=float, action="append", default=[],
                        help="render only this time in seconds to a PNG beside --output")
    args = parser.parse_args()

    scene = Scene()
    figure = plt.figure(figsize=(12.8, 7.2), dpi=150, facecolor=BG)
    ax = figure.add_axes([0, 0, 1, 1])
    height = round(args.width * 9 / 16)

    def render(t):
        scene.draw(ax, t)
        figure.canvas.draw()
        image = Image.fromarray(np.asarray(figure.canvas.buffer_rgba())[..., :3])
        return image.resize((args.width, height), Image.LANCZOS)

    if args.still:
        for t in args.still:
            path = args.output.with_name(f"{args.output.stem}_{t:05.2f}s.png")
            render(t).save(path)
            print(path)
        return 0

    if args.output.exists():
        raise SystemExit(f"{args.output} already exists; choose a new path")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # One shared palette keeps the plum background from flickering between frames.
    samples = [render(x) for x in (6.4, 10.5, 23.5, DURATION - 0.01)]
    strip = Image.new("RGB", (args.width, height * len(samples)))
    for i, image in enumerate(samples):
        strip.paste(image, (0, i * height))
    palette = strip.quantize(colors=96, method=Image.Quantize.MEDIANCUT)
    total = int(round(DURATION * args.fps))

    def frames():
        for frame in range(total):
            yield render(frame / args.fps).quantize(palette=palette, dither=Image.Dither.NONE)

    stream = frames()
    next(stream).save(args.output, save_all=True, append_images=stream,
                      duration=round(1000 / args.fps), loop=0, optimize=False)
    still = args.output.with_name(args.output.stem + "_final.png")
    render(DURATION - 0.01).save(still)
    print(f"Wrote {args.output} ({total} frames, {DURATION:.0f} s) and {still}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
