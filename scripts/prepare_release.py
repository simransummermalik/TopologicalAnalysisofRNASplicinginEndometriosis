#!/usr/bin/env python3
"""Prepare an unpublished source snapshot for local release review.

No git mutation, tagging, publishing or license choice. Only allowlisted source,
tests, guides and small demonstration inputs enter the archive; no raw data,
human reference, alignments, run folders or compiled binaries are included.
"""

import argparse
import json
import subprocess
import sys
import tarfile
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "preprocessing"))


def prepare(output):
    from hisat2_pipeline import checksum, save

    if output.exists():
        raise ValueError("output folder exists; choose a new folder")
    package = tomllib.loads((ROOT / "Cargo.toml").read_text())["package"]
    name = f"splice-girl-{package['version']}-review"
    names = [
        "Cargo.toml",
        "Cargo.lock",
        "CHANGELOG.md",
        "README.md",
        "ROADMAP.md",
        "priorities.md",
        "mermaid.md",
        "roadmap.tex",
        ".gitignore",
        "src",
        "tests",
        "scripts",
        "preprocessing",
        "gene_sets",
        "data/gene_sets",
        "docs",
        "examples/demo",
        ".github/workflows",
        "assets/splice_girl_banner.txt",
        "assets/hodgetheorynew.pdf",
    ]
    files = set()
    for entry in names + (["LICENSE"] if (ROOT / "LICENSE").is_file() else []):
        path = ROOT / entry
        candidates = path.rglob("*") if path.is_dir() else [path]
        for candidate in candidates:
            relative = candidate.relative_to(ROOT)
            if (
                candidate.is_file()
                and not candidate.is_symlink()
                and not ({"runs", "__pycache__", ".DS_Store"} & set(relative.parts))
                and candidate.suffix not in {".pyc", ".bam", ".sam", ".gz"}
            ):
                files.add(relative)
    output.mkdir(parents=True)
    archive = output / (name + ".tar.gz")
    with tarfile.open(archive, "w:gz") as bundle:
        for relative in sorted(files):
            bundle.add(
                ROOT / relative, arcname=str(Path(name) / relative), recursive=False
            )

    def git(*args):
        return subprocess.run(
            ["git", *args], cwd=ROOT, check=True, capture_output=True, text=True
        ).stdout.strip()

    record = {
        "status": "local review snapshot; unpublished",
        "version": package["version"],
        "base_commit": git("rev-parse", "HEAD"),
        "working_tree_status": git("status", "--short"),
        "license": "present; review before release"
        if (ROOT / "LICENSE").is_file()
        else "team choice pending; not ready to publish",
        "archive": checksum(archive),
        "files": [
            {**checksum(ROOT / relative), "path": str(relative)}
            for relative in sorted(files)
        ],
        "limits": "This snapshot includes uncommitted source. It is not a tagged release, a final cohort result, or a license grant. Small human-pilot artifacts remain in the repository, not this software archive.",
    }
    save(output / "release_manifest.json", record)
    print(
        json.dumps(
            {
                "archive": str(archive),
                "source_files": len(files),
                "license": record["license"],
            },
            indent=2,
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    try:
        prepare(parser.parse_args().output_dir.resolve())
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
