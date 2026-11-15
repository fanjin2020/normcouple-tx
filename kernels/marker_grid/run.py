"""Pinned environment entry point; default 36 CPU workers, one BLAS thread each."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
PACKAGES=dict(numpy='2.2.6',pandas='2.2.3',scipy='1.15.3',h5py='3.13.0',matplotlib='3.10.1')


def probe(path):
    if not path.is_file(): return False
    code='import importlib.metadata as m; import numpy,pandas,scipy,h5py,matplotlib; assert all(m.version(k)==v for k,v in '+repr(PACKAGES)+'.items())'
    try:
        return subprocess.run([str(path),'-c',code],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=45).returncode==0
    except (OSError,subprocess.TimeoutExpired): return False


def main(argv=None):
    argv=list(sys.argv[1:] if argv is None else argv)
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'):
        os.environ[name]='1'
    os.environ['MPLBACKEND']='Agg'; os.environ['PYTHONDONTWRITEBYTECODE']='1'
    if '--setup' in argv:
        argv.remove('--setup')
        p=argparse.ArgumentParser(add_help=False)
        p.add_argument('--project-root',type=Path,default=Path('/home/research/git/repository/paper1_public_evidence'))
        p.add_argument('--offline',action='store_true')
        args,_=p.parse_known_args(argv)
        suffix=Path('Scripts/python.exe' if os.name=='nt' else 'bin/python')
        candidates=[Path(sys.executable)]
        for pattern in ('paper3_mechanism_source_*','paper3_strengthening_source_*','paper3_revision_source_*'):
            candidates.extend(s/'.venv'/suffix for s in sorted(args.project_root.glob(pattern)))
        candidates.append(HERE/'.venv'/suffix)
        for python in dict.fromkeys(candidates):
            if probe(python):
                print('REUSE PYTHON',python,flush=True)
                return subprocess.call([str(python),str(HERE/'run.py'),*argv])
        if args.offline: raise SystemExit('No pinned cached Python environment; cannot install in offline mode')
        env=HERE/'.venv'
        if not env.exists(): subprocess.run([sys.executable,'-m','venv',str(env)],check=True)
        python=env/suffix
        subprocess.run([str(python),'-m','pip','install','-r',str(HERE/'engine'/'requirements.txt')],check=True)
        if not probe(python): raise SystemExit('Python environment preflight failed')
        return subprocess.call([str(python),str(HERE/'run.py'),*argv])
    sys.path.insert(0,str(HERE/'engine'))
    from driver import main as run
    return run(argv)


if __name__=='__main__':
    raise SystemExit(main())
