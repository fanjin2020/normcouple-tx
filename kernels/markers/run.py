"""Use an existing pinned Python environment; do not alter running simulations."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PACKAGES = {'numpy': '2.2.6', 'pandas': '2.2.3', 'scipy': '1.15.3', 'h5py': '3.13.0', 'matplotlib': '3.10.1'}


def probe(path):
    code = 'import importlib.metadata as m; import numpy,pandas,scipy,h5py,matplotlib; assert all(m.version(k)==v for k,v in ' + repr(PACKAGES) + '.items())'
    try:
        return subprocess.run([str(path), '-c', code], stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL, timeout=45).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def main():
    args = list(sys.argv[1:])
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'):
        os.environ[name] = '1'
    os.environ['MPLBACKEND'] = 'Agg'
    if '--setup' in args:
        args.remove('--setup')
        p = argparse.ArgumentParser(add_help=False)
        p.add_argument('--project-root', type=Path, default=Path('/home/research/git/repository/paper1_public_evidence'))
        p.add_argument('--offline', action='store_true')
        known, _ = p.parse_known_args(args)
        suffix = Path('Scripts/python.exe' if os.name=='nt' else 'bin/python')
        choices = [Path(sys.executable), HERE/'.venv'/suffix]
        for pattern in ('paper3_strengthening_source_*','paper3_mechanism_source_*','paper3_revision_source_*'):
            choices.extend(d/'.venv'/suffix for d in sorted(known.project_root.glob(pattern)))
        for python in dict.fromkeys(choices):
            if python.is_file() and probe(python):
                print('REUSE PYTHON', python, flush=True)
                return subprocess.call([str(python), str(HERE/'run.py'), *args])
        if known.offline:
            raise SystemExit('No pinned Python runtime found in offline mode')
        env = HERE/'.venv'
        if not env.exists():
            subprocess.run([sys.executable, '-m', 'venv', str(env)], check=True)
        python = env/suffix
        subprocess.run([str(python), '-m', 'pip', 'install', *[k+'=='+v for k,v in PACKAGES.items()]], check=True)
        if not probe(python):
            raise SystemExit('Python dependency check failed')
        return subprocess.call([str(python), str(HERE/'run.py'), *args])
    from driver import main as run
    return run(args)


if __name__ == '__main__':
    raise SystemExit(main())
