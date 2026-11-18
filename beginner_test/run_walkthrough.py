#!/usr/bin/env python3
"""Run a documented beginner walkthrough in a fresh local source copy.

Uses existing Rust/Python installations. No downloads, Git mutations, global
installation changes or real-read processing. Expected errors are checked as
successful demonstrations of input protection, not hidden as tool failures.
"""

import argparse
import csv
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
from zoneinfo import ZoneInfo

ROOT = Path(__file__).resolve().parents[1]


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"sha256": digest.hexdigest(), "bytes": path.stat().st_size}


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if "id" in attributes:
            self.ids.add(attributes["id"])
        for attribute in ("href", "src"):
            if attribute in attributes:
                self.links.append(attributes[attribute])


def inspect_report(folder):
    documents = {}
    for file in folder.rglob("*.html"):
        parser = Links()
        parser.feed(file.read_text())
        documents[file.resolve()] = parser
    checked = 0
    for file, parser in documents.items():
        for link in parser.links:
            parts = urlsplit(link)
            if parts.scheme or parts.netloc:
                continue
            destination = (
                (file.parent / unquote(parts.path)).resolve() if parts.path else file
            )
            assert destination.is_file(), f"missing local report link: {link}"
            if parts.fragment and destination.suffix == ".html":
                assert unquote(parts.fragment) in documents[destination].ids, link
            checked += 1
    graphs = list(folder.rglob("*.svg"))
    for graph in graphs:
        ET.parse(graph)
    return {"local_links_checked": checked, "valid_svg_files": len(graphs)}


def walkthrough(args):
    source, output = args.source.resolve(), args.output_dir.resolve()
    if output.exists():
        raise ValueError("output folder already exists; use a new folder name")
    if (
        source == output
        or source in output.parents
        and output.relative_to(source).parts[0] == "target"
    ):
        # A folder under beginner_test is intentional and reviewable; target
        # is discouraged because it would hide the user's requested record.
        raise ValueError(
            "save this review in beginner_test/runs/<new-name>, not target"
        )

    def git(*arguments):
        return subprocess.run(
            ["git", "-C", str(source), *arguments], check=True, capture_output=True
        ).stdout

    commit = git("rev-parse", "HEAD").decode().strip()
    status = git("status", "--short").decode()
    tracked = [Path(p) for p in git("ls-files", "-z").decode().split("\0") if p]
    workspace = Path(tempfile.mkdtemp(prefix="splice-girl-beginner-"))
    fresh = workspace / "source"
    fresh.mkdir()
    copied = []
    for relative in tracked:
        original = source / relative
        if original.is_symlink():
            raise ValueError(f"symlink requires separate review: {relative}")
        if not original.is_file():
            raise ValueError(f"tracked source file missing: {relative}")
        destination = fresh / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(original, destination)
        copied.append({"path": str(relative), **fingerprint(destination)})
    assert not (fresh / "target").exists()
    output.mkdir(parents=True)
    logs, artifacts = output / "logs", output / "artifacts"
    logs.mkdir()
    artifacts.mkdir()
    install = workspace / "installed"
    env = os.environ.copy()
    env.update(
        CARGO_INSTALL_ROOT=str(install),
        CARGO_NET_OFFLINE="true",
        CARGO_TERM_COLOR="never",
    )
    env["PATH"] = str(install / "bin") + os.pathsep + env["PATH"]
    record = {
        "review_type": "assistant-run beginner walkthrough; not independent human usability sign-off",
        "started_at": datetime.now(ZoneInfo("America/New_York")).isoformat(),
        "environment": platform.platform(),
        "source_base_commit": commit,
        "source_working_tree_status": status,
        "source_copy": str(fresh),
        "isolated_install": str(install),
        "scope": "fresh copy of tracked working-tree files; existing Rust/Python; offline synthetic CLI; no real-read processing",
        "runner": fingerprint(Path(__file__)),
        "source_files": copied,
        "steps": [],
        "checks": {},
        "status": "running",
    }

    def save():
        (output / "run_manifest.json").write_text(json.dumps(record, indent=2) + "\n")

    def command(
        label, argv, expected=0, contains=(), stdin=None, cwd=None, timeout=120
    ):
        number = len(record["steps"]) + 1
        started = time.perf_counter()
        result = subprocess.run(
            argv,
            cwd=cwd or fresh,
            env=env,
            input=stdin,
            capture_output=True,
            text=True,
            errors="replace",
            timeout=timeout,
        )
        text = result.stdout + result.stderr
        passed = result.returncode == expected and all(
            word in text for word in contains
        )
        log = logs / f"{number:02d}_{label}.txt"
        log.write_text(
            f"$ {shlex.join(map(str, argv))}\nWorking folder: {cwd or fresh}\n"
            f"Typed input: {stdin!r}\nExpected exit: {expected}\nActual exit: {result.returncode}\n"
            f"\n--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
        )
        record["steps"].append(
            {
                "number": number,
                "name": label,
                "command": list(map(str, argv)),
                "stdin": stdin,
                "expected_exit": expected,
                "actual_exit": result.returncode,
                "seconds": round(time.perf_counter() - started, 4),
                "status": "PASS" if passed else "FAIL",
                "log": str(log.relative_to(output)),
            }
        )
        save()
        print(f"{number:02d} {label}: {'PASS' if passed else 'FAIL'}", flush=True)
        if not passed:
            raise ValueError(f"unexpected result at {label}; see {log}")
        return result

    try:
        command("rust_version", ["rustc", "--version"], contains=("rustc ",))
        command("cargo_version", ["cargo", "--version"], contains=("cargo ",))
        command(
            "project_folder",
            [
                sys.executable,
                "-c",
                "from pathlib import Path; print(Path.cwd()); print('Cargo.toml exists:', Path('Cargo.toml').is_file()); print('README.md exists:', Path('README.md').is_file())",
            ],
            contains=("Cargo.toml exists: True", "README.md exists: True"),
        )
        command(
            "install",
            ["cargo", "install", "--path", "."],
            contains=("Installed package",),
            timeout=180,
        )
        resolved = shutil.which("splice-girl", path=env["PATH"])
        assert Path(resolved).resolve() == (install / "bin/splice-girl").resolve()
        record["checks"]["installed_executable"] = resolved
        command(
            "installed_version",
            ["splice-girl", "--version"],
            contains=("splice-girl 0.1.0",),
        )
        command(
            "installed_help",
            ["splice-girl", "--help"],
            contains=("--start", "visualize", "stats"),
        )
        command(
            "menu_cycle",
            ["cargo", "run", "--", "--start"],
            stdin="1\ntests/fixtures/cycle.tsv\n\n\nq\n",
            contains=(
                "S P L I C E   G I R L",
                "cycle_fraction: 1.0000000000",
                "checks: PASS",
            ),
        )
        command(
            "direct_exon_skip",
            ["cargo", "run", "--", "tests/fixtures/exon_skip.tsv"],
            contains=("cycle_fraction: 0.1111111111", "checks: PASS"),
        )
        command(
            "validate_junctions",
            ["cargo", "run", "--", "validate", "examples/demo/junctions.tsv"],
            contains=("validation_status: PASS",),
        )
        command(
            "intentional_manifest_mismatch",
            [
                "cargo",
                "run",
                "--",
                "validate",
                "examples/demo/junctions.tsv",
                "examples/demo/manifest.tsv",
            ],
            expected=1,
            contains=("control_03", "orphan_01", "validation_status: FAIL"),
        )
        command(
            "compare_all_genes",
            [
                "cargo",
                "run",
                "--",
                "compare",
                "examples/demo/all_genes.tsv",
                "examples/demo/manifest.tsv",
            ],
            contains=(
                "GENE_HIGH_CYCLE",
                "GENE_LOW_CYCLE",
                "0.5555555556",
                "-0.5555555556",
            ),
        )
        graph_args = [
            "cargo",
            "run",
            "--",
            "visualize",
            "case_02",
            "GENE_DEMO",
            "examples/demo/junctions.tsv",
            "splice-girl-demo.svg",
        ]
        command("draw_graph", graph_args)
        ET.parse(fresh / "splice-girl-demo.svg")
        command("protect_graph", graph_args, expected=1, contains=("already exists",))
        gene_args = [
            "cargo",
            "run",
            "--",
            "gene",
            "GENE_DEMO",
            "examples/demo/junctions.tsv",
            "examples/demo/manifest.tsv",
            "--output",
            "demo_results",
        ]
        command("save_gene", gene_args)
        gene = table(fresh / "demo_results/sample_metrics.tsv")
        assert len(gene) == 5 and all(row["checks_status"] == "PASS" for row in gene)
        command("protect_results", gene_args, expected=1, contains=("already exists",))
        command(
            "save_comparison",
            [
                "cargo",
                "run",
                "--",
                "compare",
                "examples/demo/all_genes.tsv",
                "examples/demo/manifest.tsv",
                "--output",
                "demo_all_genes",
            ],
        )
        command("rust_tests", ["cargo", "test"], contains=("0 failed",), timeout=180)
        command(
            "browser_report",
            [
                "cargo",
                "run",
                "--",
                "report",
                "examples/demo/all_genes.tsv",
                "examples/demo/manifest.tsv",
                "--output",
                "target/demo-report",
            ],
        )
        report = fresh / "target/demo-report"
        record["checks"]["report_files"] = inspect_report(report)
        assert record["checks"]["report_files"]["valid_svg_files"] == 8
        summaries = table(report / "gene_summaries.tsv")
        record["checks"]["gene_summary_columns"] = list(summaries[0])
        expected_effects = {"GENE_HIGH_CYCLE": 5 / 9, "GENE_LOW_CYCLE": -5 / 9}
        assert len(summaries) == 2
        for row in summaries:
            assert row["case_n"] == row["control_n"] == "2"
            assert (
                abs(
                    float(row["mean_difference_case_minus_control"])
                    - expected_effects[row["gene_id"]]
                )
                < 1e-9
            )
        record["checks"]["comparison_effects"] = expected_effects
        command("release_build", ["cargo", "build", "--release"])
        command(
            "statistics",
            [
                "cargo",
                "run",
                "--release",
                "--",
                "stats",
                "examples/demo/all_genes.tsv",
                "examples/demo/manifest.tsv",
                "--genes",
                "examples/demo/stats_genes.tsv",
                "--family",
                "synthetic-demo",
                "--alternative",
                "two-sided",
                "--design",
                "independent",
                "--min-per-group",
                "2",
                "--rank-by",
                "absolute-effect",
                "--output",
                "target/demo_complete",
            ],
        )
        stats = table(fresh / "target/demo_complete/statistics.tsv")
        assert len(stats) == 2
        record["checks"]["statistics_columns"] = list(stats[0])
        for row in stats:
            assert row["case_n"] == row["control_n"] == "2"
            assert row["permutations_considered"] == "6"
            assert abs(float(row["p_value"]) - 1 / 3) < 1e-9
            assert abs(float(row["adjusted_p_value"]) - 1 / 3) < 1e-9
        record["checks"]["exact_demo_p_values"] = (
            "two raw and two adjusted values equal 1/3 to exported precision"
        )
        command(
            "missing_file",
            ["splice-girl", "analyze", "does-not-exist.tsv"],
            expected=1,
            contains=("does-not-exist.tsv", "error:"),
        )
        command(
            "zero_support",
            ["splice-girl", "analyze", "tests/fixtures/zero_support.tsv"],
            expected=1,
            contains=("positive total support",),
        )
        command(
            "duplicate_edge",
            ["splice-girl", "analyze", "tests/fixtures/duplicate_edge.tsv"],
            expected=1,
            contains=("duplicate",),
        )
        command(
            "wrong_folder",
            [str(install / "bin/splice-girl"), "analyze", "tests/fixtures/cycle.tsv"],
            cwd=workspace,
            expected=1,
            contains=("tests/fixtures/cycle.tsv", "error:"),
        )
        command(
            "menu_error_recovery",
            ["splice-girl", "--start"],
            stdin="1\ndoes-not-exist.tsv\n\nq\n",
            contains=("error:", "does-not-exist.tsv"),
        )
        record["checks"]["menu_error_recovery"] = (
            "expected interactive error reported; menu accepts q and exits 0"
        )
        if args.chrome:
            command(
                "render_offline_report",
                [
                    str(args.chrome.resolve()),
                    "--headless",
                    "--disable-gpu",
                    "--no-sandbox",
                    "--no-pdf-header-footer",
                    f"--user-data-dir={workspace / 'chrome-profile'}",
                    "--no-first-run",
                    "--no-default-browser-check",
                    "--hide-scrollbars",
                    "--window-size=1400,1100",
                    f"--screenshot={artifacts / 'report-preview.png'}",
                    (report / "report.html").as_uri(),
                ],
                timeout=45,
            )
            assert (artifacts / "report-preview.png").stat().st_size > 1000
        for name, origin in {
            "one_gene": fresh / "demo_results",
            "comparison": fresh / "demo_all_genes",
            "report": report,
            "statistics": fresh / "target/demo_complete",
            "inputs": fresh / "examples/demo",
        }.items():
            shutil.copytree(origin, artifacts / name)
        shutil.copyfile(
            fresh / "splice-girl-demo.svg", artifacts / "splice-girl-demo.svg"
        )
        record["checks"]["copied_report_files"] = inspect_report(artifacts / "report")
        record["checks"]["copied_statistics_files"] = inspect_report(
            artifacts / "statistics"
        )
        record["status"] = "PASS"
    except (
        OSError,
        ValueError,
        AssertionError,
        KeyError,
        subprocess.SubprocessError,
    ) as error:
        record.update(status="FAIL", error=str(error))
    finally:
        record["finished_at"] = datetime.now(ZoneInfo("America/New_York")).isoformat()
        record["outputs"] = [
            {"path": str(p.relative_to(output)), **fingerprint(p)}
            for p in sorted(output.rglob("*"))
            if p.is_file() and p.name != "run_manifest.json"
        ]
        save()
    print(f"Walkthrough {record['status']}: {output / 'run_manifest.json'}", flush=True)
    if record["status"] != "PASS":
        print(record.get("error"), file=sys.stderr)
        return 1
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--chrome",
        type=Path,
        help="optional local Chrome executable for an offline report screenshot",
    )
    try:
        return walkthrough(parser.parse_args())
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
