#!/usr/bin/env python3
"""Review the actual crate archive, test it, and install its extracted source.

Never uploads, commits, tags or publishes. Python 3.11+ and Cargo are needed
for this developer check; the installed Rust application needs neither Python
nor the source checkout at runtime.
"""

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tarfile
import tomllib
from pathlib import Path, PurePosixPath

from check_installed_cli import check as check_installed

ROOT = Path(__file__).resolve().parents[1]


def check(output):
    output = output.resolve()
    if output.exists():
        raise ValueError("output directory already exists; choose a new folder")
    output.mkdir(parents=True)
    record = {
        "status": "running",
        "scope": "unpublished crate archive acceptance",
        "commands": [],
    }
    manifest = output / "acceptance.json"

    def call(command, cwd, name):
        item = {"command": list(map(str, command)), "cwd": str(cwd), "log": name}
        record["commands"].append(item)
        with (output / name).open("x") as log:
            process = subprocess.run(
                item["command"],
                cwd=cwd,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        item["returncode"] = process.returncode
        if process.returncode:
            print(f"--- {name}: last 100 lines ---", file=sys.stderr)
            print(
                "\n".join((output / name).read_text(errors="replace").splitlines()[-100:]),
                file=sys.stderr,
            )
            raise ValueError(f"check failed; see {output / name}")

    try:
        package = tomllib.loads((ROOT / "Cargo.toml").read_text())["package"]
        if package.get("license") != "CC-BY-4.0" or not (ROOT / "LICENSE").is_file():
            raise ValueError(
                "the approved CC BY 4.0 license must be present before release review"
            )
        call(
            ["cargo", "package", "--allow-dirty", "--locked", "--offline"],
            ROOT,
            "package.log",
        )
        name = f"{package['name']}-{package['version']}"
        archive = output / f"{name}.crate"
        shutil.copy2(ROOT / "target/package" / f"{name}.crate", archive)
        record["archive"] = {
            "path": str(archive),
            "bytes": archive.stat().st_size,
            "sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
        }
        destination = output / "source"
        destination.mkdir()
        allowed_files = {
            "Cargo.toml",
            "Cargo.toml.orig",
            "Cargo.lock",
            ".cargo_vcs_info.json",
            "CRATE_README.md",
            "CHANGELOG.md",
            "LICENSE",
            "CITATION.cff",
        }
        allowed_prefixes = (
            "src/",
            "tests/",
            "examples/demo/",
            "examples/shared_graph/",
        )
        allowed_docs = {
            "docs/data_contract.md",
            "docs/math_spec.md",
            "docs/shared_graphs.md",
            "docs/configuration.md",
            "docs/user_acceptance.md",
        }
        names = []
        with tarfile.open(archive, "r:gz") as source:
            for member in source.getmembers():
                parts = PurePosixPath(member.name).parts
                if (
                    not parts
                    or parts[0] != name
                    or ".." in parts
                    or member.issym()
                    or member.islnk()
                ):
                    raise ValueError(f"unsafe crate member: {member.name}")
                relative = PurePosixPath(*parts[1:]).as_posix()
                if member.isdir():
                    continue
                if not member.isfile() or relative.endswith(".DS_Store"):
                    raise ValueError(f"unexpected crate member: {relative}")
                if not (
                    relative in allowed_files
                    or relative in allowed_docs
                    or relative == "assets/splice_girl_banner.txt"
                    or relative.startswith(allowed_prefixes)
                ):
                    raise ValueError(
                        f"unexpected payload outside software allowlist: {relative}"
                    )
                target = destination / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                with source.extractfile(member) as contents, target.open("xb") as dest:
                    dest.write(contents.read())
                names.append(relative)
        required = {
            "LICENSE",
            "CITATION.cff",
            "CRATE_README.md",
            "src/main.rs",
            "examples/shared_graph/analysis.json",
            "examples/shared_graph/INSTALLED.md",
        }
        if not required.issubset(names):
            raise ValueError(
                f"required payload files missing: {sorted(required - set(names))}"
            )
        record["files"] = sorted(names)
        call(
            [
                "cargo",
                "test",
                "--all-targets",
                "--locked",
                "--offline",
                "--target-dir",
                ROOT / "target",
            ],
            destination,
            "archive-tests.log",
        )
        check_installed(output / "installed", source=destination)
        record.update(
            status="complete",
            license=package["license"],
            citation="CITATION.cff",
            upload_performed=False,
            file_count=len(names),
        )
    except BaseException as error:
        record.update(status="failed", error=str(error))
        raise
    finally:
        manifest.write_text(json.dumps(record, indent=2) + "\n")
    print(f"Crate archive and extracted installation passed: {manifest}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    try:
        check(parser.parse_args().output_dir)
    except (OSError, ValueError, subprocess.SubprocessError, tarfile.TarError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
