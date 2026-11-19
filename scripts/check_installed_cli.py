#!/usr/bin/env python3
"""Install into an isolated root and run the shipped guide outside the checkout.

Python orchestrates this developer check; the installed application is run
with an empty PATH and uses no Python, Cargo, external tools or source files.
Nothing is published and the user's normally installed executable is untouched.
"""

import argparse
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def check(output):
    output = output.resolve()
    if output.exists():
        raise ValueError("output folder already exists; choose a new folder")
    output.mkdir(parents=True)
    record = {
        "status": "running",
        "scope": "isolated installed CLI; synthetic data only",
        "commands": [],
    }
    manifest = output / "acceptance.json"
    temporary = None

    def call(command, cwd, name, env=None):
        item = {"command": list(map(str, command)), "cwd": str(cwd), "log": name}
        record["commands"].append(item)
        with (output / name).open("x", encoding="utf-8") as log:
            process = subprocess.run(
                item["command"],
                cwd=cwd,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        item["returncode"] = process.returncode
        if process.returncode:
            raise ValueError(f"command failed; see {output / name}")

    try:
        call(
            [
                "cargo",
                "install",
                "--path",
                ROOT,
                "--root",
                output / "install",
                "--locked",
                "--offline",
                "--target-dir",
                ROOT / "target",
            ],
            ROOT,
            "install.log",
        )
        binary = (
            output
            / "install/bin"
            / ("splice-girl.exe" if os.name == "nt" else "splice-girl")
        )
        record["installed_binary_sha256"] = digest(binary)
        temporary = tempfile.TemporaryDirectory(prefix="splice-girl-installed-")
        outside = Path(temporary.name)
        env = dict(os.environ, PATH="")
        call(
            [binary, "example", "--output", "demo-inputs"], outside, "example.log", env
        )
        demo = outside / "demo-inputs"
        guide = (demo / "README.md").read_text()
        commands = re.findall(r"```sh\n(.*?)\n```", guide, re.S)
        number = 0
        for block in commands:
            for line in block.replace("\\\n", " ").splitlines():
                if not line.strip():
                    continue
                command = shlex.split(line)
                if command[0] != "splice-girl":
                    raise ValueError("unexpected executable in installed guide")
                number += 1
                call([binary, *command[1:]], demo, f"guide-{number:02d}.log", env)
        if number != 5:
            raise ValueError(f"expected 5 installed guide commands, found {number}")
        prepared = json.loads((demo / "prepared/run_manifest.json").read_text())
        sensitivity = json.loads((demo / "sensitivity/run_manifest.json").read_text())
        if (
            prepared["summary"]["eligible_genes"] != 1
            or len(sensitivity["scenarios"]) != 7
        ):
            raise ValueError("unexpected synthetic eligibility/scenario counts")
        for path in [
            demo / "statistics/report.html",
            demo / "sensitivity/primary/statistics/report.html",
        ]:
            if not path.is_file():
                raise ValueError(f"missing report: {path}")
        for path in [
            demo / "prepared/run_manifest.json",
            demo / "sensitivity/run_manifest.json",
        ]:
            metadata = json.loads(path.read_text())
            if metadata["status"] != "complete":
                raise ValueError(f"incomplete run: {path}")
            for item in metadata["outputs"]:
                if digest(path.parent / item["path"]) != item["sha256"]:
                    raise ValueError(f"output checksum mismatch: {item['path']}")
        record.update(
            status="complete",
            installed_guide_commands=number,
            application_path="empty",
            eligible_genes=1,
            sensitivity_scenarios=7,
            zero_eligible_scenarios=sum(
                row["statistics_status"] == "no_eligible_genes"
                for row in sensitivity["scenarios"]
            ),
        )
        shutil.copytree(outside, output / "artifacts")
        record["artifacts"] = "artifacts/demo-inputs"
    except BaseException as error:
        record.update(status="failed", error=str(error))
        raise
    finally:
        if temporary is not None:
            temporary.cleanup()
        manifest.write_text(json.dumps(record, indent=2) + "\n")
    print(f"Installed workflow passed with empty PATH: {manifest}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    try:
        check(parser.parse_args().output_dir)
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
