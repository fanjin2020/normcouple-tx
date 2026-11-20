#!/usr/bin/env python3
"""One-command environment setup and execution, without modifying old outputs."""
import os
import subprocess
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS",
            "RCPP_PARALLEL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[key]="1"


def setup():
    import hashlib
    import venv
    from run_revision import ensure_cmake
    root=Path(os.environ.get("PAPER3_EXTENSION_ENV",str(HERE/".venv"))).resolve()
    python=root/("Scripts/python.exe" if os.name=="nt" else "bin/python")
    if not python.exists(): venv.EnvBuilder(with_pip=True).create(root)
    pin=hashlib.sha256((HERE/"requirements.txt").read_bytes()).hexdigest()
    receipt=root/"requirements.sha256"
    if not receipt.exists() or receipt.read_text().strip()!=pin:
        subprocess.check_call([str(python),"-m","pip","install","-r",str(HERE/"requirements.txt")])
        receipt.write_text(pin,encoding="ascii")
    ensure_cmake(python)
    args=[a for a in sys.argv[1:] if a!="--setup"]
    args.append("--install-r-packages")
    raise SystemExit(subprocess.call([str(python),str(HERE/"run_extension.py"),*args]))


if __name__=="__main__":
    if "--setup" in sys.argv: setup()
    from extension_driver import main
    raise SystemExit(main())
