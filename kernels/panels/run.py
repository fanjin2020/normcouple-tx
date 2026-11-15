"""Reuse the pinned Python/R environments; never download research inputs."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PACKAGES = dict(numpy="2.2.6", pandas="2.2.3", scipy="1.15.3", h5py="3.13.0", matplotlib="3.10.1")


def probe(path):
    if not path.is_file():
        return False
    code = "import importlib.metadata as m; import numpy,pandas,scipy,h5py,matplotlib; assert all(m.version(k)==v for k,v in " + repr(PACKAGES) + ".items())"
    try:
        return subprocess.run([str(path), "-c", code], stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL, timeout=45).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        os.environ[name] = "1"
    os.environ["MPLBACKEND"] = "Agg"
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    from cli import parser
    args = parser().parse_args(argv)
    if args.status:
        from support import progress
        return progress(args.output, args.verify)
    if args.setup:
        argv.remove("--setup")
        suffix = Path("Scripts/python.exe" if os.name == "nt" else "bin/python")
        candidates = [Path(sys.executable), HERE / ".venv" / suffix]
        for pattern in ("paper3_random_marker_source_*", "paper3_strengthening_source_*",
                        "paper3_mechanism_source_*", "paper3_revision_source_*"):
            candidates.extend(s / ".venv" / suffix for s in sorted(args.project_root.glob(pattern)))
        for python in dict.fromkeys(candidates):
            if probe(python):
                print("REUSE PYTHON", python, flush=True)
                return subprocess.call([str(python), str(HERE / "run.py"), *argv])
        if args.offline:
            raise SystemExit("No pinned cached Python environment. Retry --setup without --offline to install dependencies only.")
        env = HERE / ".venv"
        if not env.exists():
            subprocess.run([sys.executable, "-m", "venv", str(env)], check=True)
        python = env / suffix
        subprocess.run([str(python), "-m", "pip", "install", "-r", str(HERE / "engine" / "requirements.txt")], check=True)
        if not probe(python):
            raise SystemExit("Pinned Python dependency check failed")
        return subprocess.call([str(python), str(HERE / "run.py"), *argv])
    sys.path.insert(0, str(HERE / "engine"))
    from driver import main as run
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
