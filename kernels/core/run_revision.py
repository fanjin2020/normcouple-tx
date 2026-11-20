#!/usr/bin/env python3
"""Entry point; --setup installs software, never research datasets."""
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"


def cmake_version():
    exe = shutil.which("cmake")
    if exe is None:
        return None
    try:
        result = subprocess.run([exe, "--version"], capture_output=True, text=True,
                                timeout=30, check=True)
    except (OSError, subprocess.SubprocessError):
        return None
    match = re.search(r"cmake version (\d+)\.(\d+)", result.stdout)
    return tuple(map(int, match.groups())) if match else None


def ensure_cmake(python):
    # R's source-package installers inherit this PATH; pip's wheel needs no sudo.
    os.environ["PATH"] = str(python.parent) + os.pathsep + os.environ.get("PATH", "")
    version = cmake_version()
    if version is None or version < (3, 5):
        print("Installing CMake for RcppParallel in the private Python environment", flush=True)
        subprocess.check_call([str(python), "-m", "pip", "install", "--only-binary=:all:",
                               "cmake==3.31.10"])
        version = cmake_version()
    if version is None or version < (3, 5):
        raise SystemExit("CMake >=3.5 is still unavailable; inspect the pip output above.")
    print("CMAKE", shutil.which("cmake"), version, flush=True)


def setup_python():
    import venv
    if sys.version_info < (3, 10):
        raise SystemExit("Python >= 3.10 required; Python 3.11 or 3.12 recommended.")
    root = Path(os.environ.get("PAPER3_ENV_DIR", str(HERE / ".venv"))).expanduser().resolve()
    python = root / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    if not python.exists():
        venv.EnvBuilder(with_pip=True).create(root)
    requirements = HERE / "requirements.txt"
    import hashlib
    pin = hashlib.sha256(requirements.read_bytes()).hexdigest()
    receipt = root / "paper3_requirements.sha256"
    if not receipt.exists() or receipt.read_text().strip() != pin:
        subprocess.check_call([str(python), "-m", "pip", "install", "-r", str(requirements)])
        receipt.write_text(pin + "\n", encoding="ascii")
    else:
        print("REUSE Python environment", root, flush=True)
    ensure_cmake(python)
    args = [x for x in sys.argv[1:] if x != "--setup"]
    if "--install-r-packages" not in args:
        args.append("--install-r-packages")
    raise SystemExit(subprocess.call([str(python), str(HERE / "run_revision.py"), *args]))


if __name__ == "__main__":
    if "--setup" in sys.argv:
        setup_python()
    from revision_driver import main
    raise SystemExit(main())
