#!/usr/bin/env python3
"""Render a slow, plain-language animation of the graph Hodge decomposition for slides.

The example is the repository's exon-skipping fixture,
tests/fixtures/exon_skip.tsv, the same input as
`splice-girl tests/fixtures/exon_skip.tsv`. The decomposition is recomputed
here with NumPy least squares and checked against the fixture's reviewed
answer (cycle fraction 1/9) before anything is drawn, so arrow widths and
the final percentages are the real values.

The style is black and white with Computer Modern type. Arrow width shows
how much flow an edge carries; a dashed arrow carries flow backward, matching
`splice-girl visualize`. The animation runs in nine numbered steps, about
72 seconds at the default pace; `--pace 1.3` plays it 30% slower.

Requires matplotlib, NumPy and Pillow:

    python3 scripts/make_hodge_animation.py \
        --output assets/presentation/hodge_decomposition.gif
"""

from __future__ import annotations

import argparse
import csv
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle  # noqa: E402
from PIL import Image  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "exon_skip.tsv"
WIDTH, HEIGHT = 16.0, 9.0
INK, GRAY, LIGHT = "#111111", "#6f6f6f", "#d4d4d4"

# Start and end of each numbered step, in seconds at --pace 1.
STEPS = [(3, 10), (10, 17), (17, 24), (24, 31), (31, 39),
         (39, 47), (47, 55), (55, 63), (63, 72)]
DURATION = float(STEPS[-1][1])

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

# Local layout of the three-exon skipping graph (unit scale).
NODE_X = {"E1": -2.2, "E2": 0.0, "E3": 2.2}
BOX_W, BOX_H = 1.0, 0.56


def line(p0, p1, n=60):
    return np.linspace(p0, p1, n)


def bezier(p0, control, p1, n=120):
    t = np.linspace(0.0, 1.0, n)[:, None]
    p0, control, p1 = map(np.asarray, (p0, control, p1))
    return (1 - t) ** 2 * p0 + 2 * (1 - t) * t * control + t**2 * p1


PATHS = {
    ("E1", "E2"): line((-1.7 + 0.04, 0.0), (-0.5 - 0.03, 0.0)),
    ("E2", "E3"): line((0.5 + 0.04, 0.0), (1.7 - 0.03, 0.0)),
    ("E1", "E3"): bezier((-2.2, 0.30), (0.0, 2.3), (2.2, 0.31)),
}
INCLUSION = [("E1", "E2"), ("E2", "E3")]
SKIP = ("E1", "E3")


def load_and_decompose(path: Path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    nodes = sorted({r["from"] for r in rows} | {r["to"] for r in rows})
    edges = sorted((r["from"], r["to"]) for r in rows)
    if nodes != ["E1", "E2", "E3"] or set(edges) != set(PATHS):
        raise SystemExit("this animation is laid out for the three-exon skipping fixture")
    counts = {(r["from"], r["to"]): float(r["count"]) for r in rows}
    index = {node: i for i, node in enumerate(nodes)}
    incidence = np.zeros((len(nodes), len(edges)))
    for column, (source, target) in enumerate(edges):
        incidence[index[source], column] = -1.0
        incidence[index[target], column] = 1.0
    raw = np.array([counts[e] for e in edges])
    signal = raw / raw.sum()
    potential = np.linalg.lstsq(incidence.T, signal, rcond=None)[0]
    gradient = incidence.T @ potential
    cycle = signal - gradient
    fraction = float(cycle @ cycle / (signal @ signal))
    checks = (
        np.allclose(gradient + cycle, signal, atol=1e-12)
        and np.allclose(incidence @ cycle, 0.0, atol=1e-12)
        and abs(float(gradient @ cycle)) < 1e-12
        and abs(fraction - 1 / 9) < 1e-12
    )
    if not checks:
        raise SystemExit("decomposition does not match the reviewed fixture answer (1/9)")
    as_map = lambda values: dict(zip(edges, map(float, values)))  # noqa: E731
    return {
        "signal": as_map(signal),
        "gradient": as_map(gradient),
        "cycle": as_map(cycle),
        "fraction": fraction,
    }


def smooth(x: float) -> float:
    x = min(max(x, 0.0), 1.0)
    return x * x * (3 - 2 * x)


def ramp(t: float, start: float, end: float) -> float:
    return smooth((t - start) / (end - start))


def window(t, start, end, fade_in=0.7, fade_out=0.5, last=False):
    value = ramp(t, start, start + fade_in)
    return value if last else value * (1 - ramp(t, end - fade_out, end))


def dim(t, start, end, low=0.35, fade=0.7):
    """1 outside [start, end], easing down to `low` inside it."""
    return 1 - (1 - low) * ramp(t, start, start + fade) * (1 - ramp(t, end, end + fade))


def lerp(a, b, u):
    return a + (b - a) * u


class Scene:
    def __init__(self, data):
        self.data = data
        self.vmax = max(
            abs(v) for key in ("signal", "gradient", "cycle") for v in data[key].values()
        )

    # ------------------------------------------------------------ primitives
    @staticmethod
    def place(points, center, scale):
        return np.asarray(center) + scale * np.asarray(points)

    def width(self, value, scale):
        return (0.9 + 7.0 * abs(value) / self.vmax) * scale

    def edge(self, ax, key, value, center, scale, alpha, progress=1.0):
        if alpha <= 0.01 or progress <= 0.0:
            return
        points = self.place(PATHS[key], center, scale)
        points = points[: max(2, int(round(progress * len(points))))]
        style = "-" if value >= -1e-12 else (0, (3.2, 2.2))
        ax.plot(points[:, 0], points[:, 1], color=INK, lw=self.width(value, scale), ls=style,
                alpha=alpha, solid_capstyle="butt", dash_capstyle="butt", zorder=2)
        tip, before = points[-1], points[-min(4, len(points))]
        direction = (tip - before) / np.linalg.norm(tip - before)
        normal = np.array([-direction[1], direction[0]])
        length, half = 0.22 * scale, 0.11 * scale
        head = [tip, tip - length * direction + half * normal, tip - length * direction - half * normal]
        ax.add_patch(Polygon(head, closed=True, color=INK, alpha=alpha, lw=0, zorder=3))

    def graph(self, ax, values, center, scale, alpha, progress=1.0):
        for key in PATHS:
            self.edge(ax, key, values[key], center, scale, alpha, progress)

    def nodes(self, ax, center, scale, alphas):
        for node, x in NODE_X.items():
            alpha = alphas[node]
            if alpha <= 0.01:
                continue
            cx, cy = self.place((x, 0.0), center, scale)
            ax.add_patch(FancyBboxPatch(
                (cx - BOX_W * scale / 2, cy - BOX_H * scale / 2), BOX_W * scale, BOX_H * scale,
                boxstyle=f"round,pad=0,rounding_size={0.07 * scale}",
                facecolor="white", edgecolor=INK, lw=1.3 * scale, alpha=alpha, zorder=4))
            ax.text(cx, cy, f"$E_{node[1]}$", ha="center", va="center",
                    fontsize=21 * scale, color=INK, alpha=alpha, zorder=5)

    def beads(self, ax, values, center, scale, alpha, elapsed):
        if alpha <= 0.01:
            return
        for key, value in values.items():
            points = self.place(PATHS[key], center, scale)
            along = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))])
            total = along[-1]
            offset = np.sign(value) * (0.3 + 1.4 * abs(value) / self.vmax) * scale * elapsed
            spacing = 0.62 * scale
            for k in range(int(total / spacing) + 1):
                s = (k * spacing + offset) % total
                fade = min(1.0, s / (0.3 * scale), (total - s) / (0.45 * scale))
                if fade <= 0.05:
                    continue
                x = np.interp(s, along, points[:, 0])
                y = np.interp(s, along, points[:, 1])
                ax.scatter([x], [y], s=(self.width(value, scale) + 4.0 * scale) ** 2,
                           facecolor="white", edgecolor=INK, linewidths=1.1 * scale,
                           alpha=alpha * fade, zorder=3.5)

    def circulation(self, ax, center, scale, alpha):
        if alpha <= 0.01:
            return
        mid = self.place((0.0, 0.86), center, scale)
        radius = 0.27 * scale
        theta = np.radians(np.linspace(-50, 215, 80))
        xs, ys = mid[0] + radius * np.cos(theta), mid[1] + radius * np.sin(theta)
        ax.plot(xs, ys, color=INK, lw=1.4 * scale, alpha=alpha, zorder=2)
        tip = np.array([xs[-1], ys[-1]])
        direction = np.array([-np.sin(theta[-1]), np.cos(theta[-1])])
        normal = np.array([-direction[1], direction[0]])
        length, half = 0.17 * scale, 0.085 * scale
        head = [tip + 0.5 * length * direction, tip - 0.5 * length * direction + half * normal,
                tip - 0.5 * length * direction - half * normal]
        ax.add_patch(Polygon(head, closed=True, color=INK, alpha=alpha, lw=0, zorder=3))

    # ------------------------------------------------------------ the frame
    def draw(self, ax, t):
        d = self.data
        ax.clear()
        ax.set_xlim(0, WIDTH)
        ax.set_ylim(0, HEIGHT)
        ax.set_aspect("equal")
        ax.axis("off")

        header = ramp(t, 0.0, 1.2)
        ax.text(0.6, 8.42, "Hodge decomposition of a splice graph", fontsize=29,
                color=INK, alpha=header, va="center")
        ax.plot([0.6, 15.4], [7.98, 7.98], color=INK, lw=0.8, alpha=header)
        for number, (start, end) in enumerate(STEPS, 1):
            alpha = window(t, start, end, 0.5, 0.4, last=number == len(STEPS))
            if alpha > 0.01:
                ax.text(15.4, 8.42, f"step {number} of {len(STEPS)}", fontsize=17,
                        color=GRAY, alpha=alpha, ha="right", va="center")

        intro = window(t, 0.4, STEPS[0][0], fade_in=0.8, fade_out=0.6)
        ax.text(8.0, 4.95, "How Hodge decomposition looks at splicing", fontsize=28,
                color=INK, alpha=intro, ha="center", va="center")
        ax.text(8.0, 4.2, "one gene, step by step", fontsize=19,
                color=GRAY, alpha=intro, ha="center", va="center")

        big = ((8.0, 4.6), 1.3)
        panels = [((2.7, 5.0), 0.8), ((8.0, 5.0), 0.8), ((13.3, 5.0), 0.8)]
        dim_observed = dim(t, 31.0, 47.0)
        dim_gradient = dim(t, 39.0, 47.0)

        # Steps 1-3: build the measured graph large; step 4 moves it to the left.
        move = ramp(t, 24.4, 25.8)
        oc = tuple(lerp(np.array(big[0]), np.array(panels[0][0]), move))
        os_ = lerp(big[1], panels[0][1], move)
        drawn = {INCLUSION[0]: (10.8, 12.2), INCLUSION[1]: (12.6, 14.0), SKIP: (17.6, 19.4)}
        for key, (start, end) in drawn.items():
            self.edge(ax, key, d["signal"][key], oc, os_, dim_observed, progress=ramp(t, start, end))
        self.nodes(ax, oc, os_, {n: ramp(t, 3.8 + i, 4.6 + i) * dim_observed
                                 for i, n in enumerate(NODE_X)})
        self.beads(ax, d["signal"], oc, os_, ramp(t, 20.4, 21.2) * dim_observed * 0.95, t - 20.4)
        x, y = self.place((0.0, 1.62), oc, os_)
        ax.text(x, y, r"skips $E_2$", ha="center", va="center", fontsize=20 * os_, color=GRAY,
                alpha=window(t, 19.0, 24.4, fade_in=0.8, fade_out=0.6))

        # Step 4 sets up the two empty parts; steps 5 and 6 fill them in.
        placeholder = 0.3 * ramp(t, 26.4, 27.2)
        gc, gs = panels[1]
        self.nodes(ax, gc, gs, {n: lerp(placeholder, 1.0, ramp(t, 31.2, 32.0)) * dim_gradient
                                for n in NODE_X})
        self.graph(ax, d["gradient"], gc, gs, dim_gradient, progress=ramp(t, 31.8, 33.6))
        self.beads(ax, d["gradient"], gc, gs, ramp(t, 33.8, 34.6) * dim_gradient * 0.95, t - 33.8)

        cc, cs = panels[2]
        self.nodes(ax, cc, cs, {n: lerp(placeholder, 1.0, ramp(t, 39.2, 40.0)) for n in NODE_X})
        self.graph(ax, d["cycle"], cc, cs, 1.0, progress=ramp(t, 39.8, 41.6))
        self.beads(ax, d["cycle"], cc, cs, ramp(t, 41.8, 42.6) * 0.95, t - 41.8)
        self.circulation(ax, cc, cs, ramp(t, 42.2, 43.2))
        x, y = self.place((0.0, -0.66), cc, cs)
        ax.text(x, y, "dashed arrow: flow runs backward", ha="center", va="center", fontsize=15,
                color=GRAY, alpha=ramp(t, 43.4, 44.4))

        for (pc, ps), title, subtitle, alpha in (
            (panels[0], "what we measured", "", ramp(t, 25.6, 26.4) * dim_observed),
            (panels[1], "straight-through part", "the gradient", ramp(t, 26.4, 27.2) * dim_gradient),
            (panels[2], "loop part", "the cycle", ramp(t, 27.0, 27.8)),
        ):
            ax.text(pc[0], 7.15, title, ha="center", va="center", fontsize=21, color=INK, alpha=alpha)
            if subtitle:
                ax.text(pc[0], 6.68, subtitle, ha="center", va="center", fontsize=16,
                        color=GRAY, alpha=alpha)
        signs = ramp(t, 27.6, 28.4)
        for x, symbol in ((5.35, "$=$"), (10.65, "$+$")):
            ax.text(x, panels[0][0][1] + 0.35 * panels[0][1], symbol, ha="center",
                    va="center", fontsize=36, color=INK, alpha=signs)

        # Step 8: how the flow divides between the two parts.
        grow = ramp(t, 55.6, 57.6)
        if grow > 0:
            left, total_width, y, h = 3.2, 9.6, 2.35, 0.38
            loop_share = d["fraction"]
            g_width = total_width * (1 - loop_share)
            c_width = total_width * loop_share
            shown = grow * total_width
            ax.add_patch(Rectangle((left, y), min(shown, g_width), h, facecolor=LIGHT,
                                   edgecolor=INK, lw=1.0, zorder=2))
            if shown > g_width:
                ax.add_patch(Rectangle((left + g_width, y), min(shown - g_width, c_width), h,
                                       facecolor="white", edgecolor=INK, hatch="////", lw=1.0, zorder=2))
            text = ramp(t, 57.2, 58.2)
            ax.text(left - 0.25, y + h / 2, "all the flow", ha="right", va="center",
                    fontsize=19, color=INK, alpha=text)
            ax.text(left + g_width / 2, y + h + 0.38,
                    f"straight through: {100 * (1 - loop_share):.0f}%",
                    ha="center", va="center", fontsize=19, color=INK, alpha=text)
            ax.text(left + g_width + c_width / 2, y + h + 0.38,
                    f"loop: {100 * loop_share:.0f}%",
                    ha="center", va="center", fontsize=19, color=INK, alpha=text)

        percent = f"{100 * d['fraction']:.0f}%"
        captions = [
            ("A gene is a chain of exons.",
             "Each box is one exon."),
            ("RNA-seq reads show which exons are joined together.",
             "Each arrow is one splice junction."),
            (r"Sometimes exon $E_2$ is skipped.",
             r"So there are two routes from $E_1$ to $E_3$: through $E_2$, or around it."),
            ("Hodge decomposition splits this flow into two simple parts.",
             "The split is exact, and the two parts never overlap."),
            ("Part 1: flow that runs straight through the gene.",
             "It moves from the first exon toward the last and never comes back around."),
            ("Part 2: flow that goes around a loop.",
             r"The loop only exists because $E_2$ can be skipped."),
            ("Add the two parts together and you get back exactly what we measured.",
             "So we can ask: how much of the flow goes straight through, and how much loops?"),
            ("The cycle fraction is the share of the flow that goes around loops.",
             f"In this example it is {percent}."),
            ("In the study, we compare this number between endometriosis and control samples.",
             "It describes alternative splicing routes. It is not evidence of circular RNA."),
        ]
        for number, ((start, end), (main, secondary)) in enumerate(zip(STEPS, captions), 1):
            alpha = window(t, start, end, 0.6, 0.5, last=number == len(STEPS))
            if alpha > 0.01:
                ax.text(8.0, 1.28, main, ha="center", va="center", fontsize=23, color=INK, alpha=alpha)
                ax.text(8.0, 0.66, secondary, ha="center", va="center", fontsize=17,
                        color=GRAY, alpha=alpha)
        ax.text(15.4, 0.2, "Example: the exon-skipping test fixture, 100 reads on each junction; "
                "values match splice-girl output.", ha="right", va="center", fontsize=10.5,
                color=GRAY, alpha=ramp(t, 64.0, 65.0))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, required=True, help="new .gif path")
    parser.add_argument("--frames-dir", type=Path, help="optional new folder for full-size PNG frames")
    parser.add_argument("--gif-width", type=int, default=1280)
    parser.add_argument("--fps", type=float, default=20.0,
                        help="GIF frame rate; GIF delays are whole 1/100 s, so use 10, 12.5, 20 or 25")
    parser.add_argument("--pace", type=float, default=1.0,
                        help="stretch factor for the whole timeline; 1.3 plays 30%% slower")
    parser.add_argument("--still", type=float, action="append", default=[],
                        help="render only this scene time in seconds to a PNG beside --output")
    args = parser.parse_args()

    scene = Scene(load_and_decompose(FIXTURE))
    figure = plt.figure(figsize=(12.8, 7.2), dpi=150, facecolor="white")
    ax = figure.add_axes([0, 0, 1, 1])

    def render(t: float) -> Image.Image:
        scene.draw(ax, t)
        figure.canvas.draw()
        return Image.fromarray(np.asarray(figure.canvas.buffer_rgba())[..., :3])

    if args.still:
        for t in args.still:
            path = args.output.with_name(f"{args.output.stem}_{t:05.2f}s.png")
            render(t).save(path)
            print(path)
        return 0

    if args.output.exists():
        raise SystemExit(f"{args.output} already exists; choose a new path")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.frames_dir:
        args.frames_dir.mkdir(parents=True, exist_ok=False)
    height = round(args.gif_width * 9 / 16)
    total = int(round(DURATION * args.pace * args.fps))

    def frames():
        for frame in range(total):
            image = render(frame / args.fps / args.pace)
            if args.frames_dir:
                image.save(args.frames_dir / f"frame_{frame:04d}.png", compress_level=1)
            yield image.resize((args.gif_width, height), Image.LANCZOS).convert("L")

    stream = frames()
    next(stream).save(args.output, save_all=True, append_images=stream,
                      duration=round(1000 / args.fps), loop=0, optimize=True)
    still = args.output.with_name(args.output.stem + "_final.png")
    render(DURATION - 0.01).save(still)
    print(f"Wrote {args.output} ({total} frames, {DURATION * args.pace:.0f} s) and {still}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
