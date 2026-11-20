"""Offline regression tests for installation and failure reporting."""
import contextlib
import io
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import run_revision as setup
from revision_driver import run_r


class SetupTests(unittest.TestCase):
    def test_cmake_probe(self):
        with patch.object(setup.shutil, "which", return_value="cmake"), \
             patch.object(setup.subprocess, "run", return_value=subprocess.CompletedProcess(
                 [], 0, "cmake version 3.31.10\n", "")):
            self.assertEqual(setup.cmake_version(), (3, 31))

    def test_missing_cmake(self):
        with patch.object(setup.shutil, "which", return_value=None):
            self.assertIsNone(setup.cmake_version())

    def test_keep_existing_cmake(self):
        python = Path(sys.executable)
        with patch.dict(os.environ), patch.object(setup, "cmake_version", return_value=(3, 31)), \
             patch.object(setup.subprocess, "check_call") as install:
            setup.ensure_cmake(python)
            install.assert_not_called()
            self.assertEqual(os.environ["PATH"].split(os.pathsep)[0], str(python.parent))

    def test_install_cmake_without_reinstalling_requirements(self):
        python = Path(sys.executable)
        with patch.dict(os.environ), patch.object(setup, "cmake_version", side_effect=[None, (3, 31)]), \
             patch.object(setup.subprocess, "check_call") as install:
            setup.ensure_cmake(python)
            install.assert_called_once_with([str(python), "-m", "pip", "install", 
                                             "--only-binary=:all:", "cmake==3.31.10"])

    def test_abort_if_cmake_still_missing(self):
        with patch.dict(os.environ), patch.object(setup, "cmake_version", return_value=None), \
             patch.object(setup.subprocess, "check_call"), self.assertRaises(SystemExit):
            setup.ensure_cmake(Path(sys.executable))

    def test_r_failure_log_visible_and_previous_log_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "install_R.log"
            log.write_text("original error", encoding="utf-8")
            err = io.StringIO()
            with contextlib.redirect_stderr(err), self.assertRaisesRegex(ValueError, "R infrastructure"):
                run_r([sys.executable, "-c", "print('ERROR missing cmake'); raise SystemExit(1)"],
                      tmp, log)
            self.assertIn("ERROR missing cmake", err.getvalue())
            self.assertEqual(next(Path(tmp).glob("*.previous.*")).read_text(), "original error")

    def test_r_timeout_keeps_diagnostics(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "R.log"
            err = io.StringIO()
            with patch("revision_driver.subprocess.run", side_effect=subprocess.TimeoutExpired("Rscript", 1)), \
                 contextlib.redirect_stderr(err), self.assertRaises(subprocess.TimeoutExpired):
                run_r(["Rscript"], tmp, log, timeout=1)
            self.assertIn("R log tail", err.getvalue())
            self.assertTrue(log.exists())


if __name__ == "__main__":
    unittest.main()
