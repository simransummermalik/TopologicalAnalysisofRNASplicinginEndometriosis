"""Failed release subprocesses must expose their error and retain their status."""

import contextlib
import io
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import check_crate
import check_installed_cli


class ReleaseDiagnosticsTests(unittest.TestCase):
    def test_failed_subprocess_is_visible_and_acceptance_remains_failed(self):
        for module, log_name in [
            (check_crate, "package.log"),
            (check_installed_cli, "install.log"),
        ]:
            with self.subTest(script=module.__name__), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                (root / "Cargo.toml").write_text('[package]\nlicense = "MIT"\n')
                (root / "LICENSE").write_text("fixture license\n")
                output = root / "review"
                stderr = io.StringIO()

                def fail(command, **kwargs):
                    kwargs["stdout"].write("underlying compiler error\n")
                    return subprocess.CompletedProcess(command, 101)

                with (
                    patch.object(module, "ROOT", root),
                    patch.object(module.subprocess, "run", side_effect=fail),
                    contextlib.redirect_stderr(stderr),
                    self.assertRaises(ValueError),
                ):
                    module.check(output)

                self.assertIn("underlying compiler error", stderr.getvalue())
                self.assertIn(log_name, stderr.getvalue())
                self.assertIn("underlying compiler error", (output / log_name).read_text())
                record = json.loads((output / "acceptance.json").read_text())
                self.assertEqual(record["status"], "failed")
                self.assertEqual(record["commands"][0]["returncode"], 101)


if __name__ == "__main__":
    unittest.main()
