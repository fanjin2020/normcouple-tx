"""Reuse installed environments only; terminate only R children owned by this run."""
import importlib.metadata
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
from common import HERE, require, read_json

PINNED = dict(numpy='2.2.6',pandas='2.2.3',scipy='1.15.3',h5py='3.13.0',matplotlib='3.10.1')


def choose_python(args, argv):
    probe = ('import importlib,importlib.metadata as m; p='+repr(PINNED)+
             '; [importlib.import_module(k) for k in p]; assert {k:m.version(k) for k in p}==p')
    candidates = [args.python] if args.python else [Path(sys.executable)]
    suffix = Path('Scripts/python.exe' if os.name=='nt' else 'bin/python')
    if not args.python:
        candidates.append(args.project_root/'.venv'/suffix)
        for pattern in ('paper3_mechanism_source_*','paper3_random_marker_source_*','paper3_revision_source_*'):
            candidates.extend(d/'.venv'/suffix for d in sorted(args.project_root.glob(pattern)))
    for candidate in dict.fromkeys(candidates):
        if candidate and candidate.is_file():
            try:
                check = subprocess.run([str(candidate),'-c',probe],capture_output=True,timeout=90)
            except (OSError,subprocess.TimeoutExpired):
                continue
            if check.returncode==0:
                # Symlinked venv interpreters share a binary but not site-packages.
                if os.path.normcase(os.path.abspath(candidate))==os.path.normcase(os.path.abspath(sys.executable)):
                    return None
                require(not os.environ.get('P3_REVIEW_REEXEC'),'Selected Python failed re-entry')
                print('REUSE PYTHON',candidate,flush=True)
                return subprocess.call([str(candidate),str(HERE/'run_coupling.py'),*argv],
                    env=dict(os.environ,P3_REVIEW_REEXEC='1'))
    raise ValueError('No compatible existing Python. No packages installed. Use --python /path/to/pinned/python')


def run_r(command, library, log, timeout=14400):
    log = Path(log)
    log.parent.mkdir(parents=True,exist_ok=True)
    env = dict(os.environ,R_LIBS_USER=str(library))
    with log.open('w',encoding='utf-8') as stream:
        child = subprocess.Popen(command,env=env,stdout=stream,stderr=subprocess.STDOUT,
                                 start_new_session=os.name!='nt')
        try:
            code = child.wait(timeout=timeout)
        except BaseException:
            if child.poll() is None:
                if os.name!='nt':
                    os.killpg(child.pid,signal.SIGTERM)
                else:
                    child.terminate()
                try:
                    child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    if os.name!='nt':
                        os.killpg(child.pid,signal.SIGKILL)
                    else:
                        child.kill()
                    child.wait()
            raise
    if code:
        print('\n'.join(log.read_text(encoding='utf-8',errors='replace').splitlines()[-35:]),flush=True)
        raise RuntimeError('R failed; see '+str(log))


def existing_r(args):
    exe = shutil.which(str(args.rscript))
    require(exe is not None,'Rscript missing. Use --rscript; packages will not be installed.')
    explicit = [p.resolve() for p in args.r_library if p.is_dir()]
    candidates = []
    for root in (args.random_results,args.project_root,args.project_root.parent):
        candidates.extend([root/'R_library',root/'runtime'/'R_library'])
        candidates.extend(sorted(root.glob('*/R_library')))
        candidates.extend(sorted(root.glob('*/runtime/R_library')))
    for key in ('R_LIBS','R_LIBS_USER'):
        candidates.extend(Path(s) for s in os.environ.get(key,'').split(os.pathsep) if s)
    candidates = list(dict.fromkeys(p.resolve() for p in candidates if p.is_dir()))
    preferred = [p for p in candidates if (p/'dtangle'/'DESCRIPTION').is_file() and
                 'Version: 2.0.9' in (p/'dtangle'/'DESCRIPTION').read_text(errors='replace')]
    libraries = list(dict.fromkeys([*explicit,*preferred,*candidates]))
    os.environ['R_LIBS'] = os.pathsep.join(map(str,libraries))
    own = args.output/'R_library'
    own.mkdir(exist_ok=True)
    path = args.output/'R_versions.json'
    run_r([exe,'--vanilla',str(HERE/'bridge.R'),'--check',str(own),str(path)],own,args.output/'R_check.log',180)
    print('Rscript',exe,'\nEXISTING R LIBRARIES',*libraries,sep='\n',flush=True)
    return dict(executable=exe,library=str(own),versions=read_json(path))


def recorded_inputs(args):
    path = args.random_results/'input_paths.json'
    if not path.exists():
        return
    entries = {x['name']:Path(x['path']) for x in read_json(path)}
    for attr,key in [('libd_dir','LIBD/counts.mtx'),('seaad_dir','SEAAD/counts.mtx'),('gene_map','gene_map')]:
        if getattr(args,attr) is None and key in entries and entries[key].exists():
            setattr(args,attr,entries[key] if attr=='gene_map' else entries[key].parent)
    key = 'external/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5'
    if args.external_cache is None and key in entries and entries[key].exists():
        args.external_cache = entries[key].parent
