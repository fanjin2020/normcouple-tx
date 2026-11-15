"""Immutable parent verification and a short, resumable follow-up workspace."""
import contextlib
import hashlib
import importlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
BASE = HERE / 'baseline'
PARENT_SIGNATURE = 'bdd118b1375add8e9fe5ed708b46587834a2312f689de1e2704e7dbc0701b211'
PARENT_SOURCE = 'f05bd8bd6b8254b6b7031c2eb14164cecab2029fa175a13f6e5862ea9f399686'
SERVER = Path('/home/research/git/repository')


def add_paths():
    for p in reversed([BASE, BASE / 'frozen' / 'engine', BASE / 'frozen']):
        sys.path.insert(0, str(p))


add_paths()
from common import require, read_json, save_json, sha, digest


def safe_file(root, name):
    require(isinstance(name, str) and name and '\\' not in name, 'Invalid manifest path')
    path = (root / name).resolve()
    require(path.is_relative_to(root.resolve()) and path.is_file(), 'Missing/unsafe file: ' + str(path))
    return path


def verify_manifest(root, name='files_sha256.json'):
    manifest = read_json(root / name)
    require(isinstance(manifest, dict) and bool(manifest), 'Empty manifest: ' + str(root))
    for rel, expected in manifest.items():
        require(sha(safe_file(root, rel)) == expected, 'Checksum mismatch: ' + rel)
    return manifest


def verify_source():
    manifest = verify_manifest(HERE, 'source_sha256.json')
    import run_coupling
    require(run_coupling.source_signature() == PARENT_SOURCE, 'Frozen parent source changed')
    return digest(manifest)


def parent_protocol(root):
    cfg = read_json(root / 'protocol.json')
    body = {k: v for k, v in cfg.items() if k != 'signature'}
    require(digest(body) == cfg.get('signature') == PARENT_SIGNATURE,
            'Expected completed bdd118b137 run, not a different experiment: ' + str(root))
    require(cfg.get('source_signature') == PARENT_SOURCE, 'Parent source mismatch')
    require(not cfg['settings']['software_fixture'], 'A software fixture is not research data')
    return cfg


def find_parent(explicit, project):
    if explicit:
        path = explicit.expanduser().absolute()
        parent_protocol(path)
        require((path / 'final').is_dir() and (path / 'clinical' / 'complete.json').is_file(),
                'Supply the FULL SERVER RUN directory, not the compact p3_review export')
        return path
    candidates = set()
    for base in {project.parent, project}:
        if not base.is_dir():
            continue
        candidates.update(p for p in base.iterdir() if p.is_dir() and p.name.startswith(('p3', 'paper3')))
    found = []
    for path in sorted(candidates):
        try:
            if (path / 'clinical' / 'complete.json').is_file():
                parent_protocol(path)
                found.append(path.absolute())
        except (OSError, ValueError, KeyError, TypeError):
            continue
    require(len(found) == 1, 'Specify --previous-results FULL_SERVER_RUN. Matching runs: ' +
            (', '.join(map(str, found)) or 'none found; only immediate project/server children were scanned'))
    return found[0]


@contextlib.contextmanager
def parent_read_lock(parent):
    require((parent / 'run.lock').is_file(), 'Parent writer lock is missing')
    with (parent / 'run.lock').open('rb') as handle:
        if os.name == 'nt':
            raise RuntimeError('Research runs require Linux flock; Windows is for software tests only')
        import fcntl
        try:
            fcntl.flock(handle, fcntl.LOCK_SH | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError('Parent run is still writing. Wait for completion; do not delete run.lock') from exc
        yield


def verify_parent(parent):
    from support import checkpoint
    cfg = parent_protocol(parent)
    manifest = verify_manifest(parent / 'final')
    status = read_json(parent / 'final' / 'analysis_scope.json')
    require(status.get('simulation_complete') and status.get('clinical_complete'), 'Parent run is incomplete')
    require(read_json(parent / 'final' / 'protocol.json') == cfg, 'Root/final parent protocols disagree')
    receipt = checkpoint(parent / 'clinical', cfg['signature'])
    require(receipt, 'Missing verified clinical receipt')
    return cfg, dict(protocol=cfg['signature'], final_manifest=digest(manifest),
                     clinical_receipt=sha(parent / 'clinical' / 'complete.json'))


def choose_python(args, parent, argv):
    from runtime import PINNED
    probe = ("import importlib,json; pins=" + repr(PINNED) +
             "; got={n:importlib.import_module(n).__version__ for n in pins}; "
             "print(json.dumps(got)); assert got==pins,(got,pins)")
    runtime = read_json(parent / 'runtime.json')
    candidates = [args.python] if args.python else [runtime.get('python'), sys.executable]
    for prefix in (args.project_root, args.project_root.parent / 'paper1-from-publicV1'):
        candidates.extend(prefix / v / 'bin' / 'python' for v in ('.venv', '.raw-reproduction-venv'))
    seen = set()
    for candidate in candidates:
        if not candidate:
            continue
        # Do not resolve venv symlinks: their containing directory selects site-packages.
        path = str(Path(candidate).expanduser().absolute())
        if path in seen or not Path(path).is_file():
            continue
        seen.add(path)
        test = subprocess.run([path, '-c', probe], text=True, capture_output=True, timeout=90)
        if test.returncode:
            print('PYTHON NOT SELECTED', path, test.stderr[-900:], flush=True)
            continue
        if Path(sys.executable).absolute() == Path(path):
            print('REUSE PYTHON', path, flush=True)
            return None
        print('REUSE PYTHON', path, flush=True)
        if os.name != 'nt':
            os.execv(path, [path, str(HERE / 'finish.py'), *argv])
        return subprocess.run([path, str(HERE / 'finish.py'), *argv], check=False).returncode
    raise RuntimeError('No pinned existing Python found. NO packages installed. Use --python /existing/venv/bin/python')


def assert_separate(output, *inputs):
    out = output.resolve()
    for root in inputs:
        p = root.resolve()
        require(not (out == p or out.is_relative_to(p) or p.is_relative_to(out)),
                'Output overlaps source/input tree: ' + str(p))


def write_manifest(root):
    manifest = {p.relative_to(root).as_posix(): sha(p) for p in sorted(root.rglob('*'))
                if p.is_file() and p.name != 'files_sha256.json'}
    save_json(root / 'files_sha256.json', manifest)
    return manifest


def export_review(output):
    cfg = read_json(output / 'protocol.json')
    root = output / 'review'
    verify_manifest(root)
    tag = 'p3n_' + cfg['signature'][:10]
    target = output.parent / (tag + '.zip')
    temporary = target.with_suffix('.zip.partial')
    with zipfile.ZipFile(temporary, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for p in sorted(root.rglob('*')):
            if p.is_file():
                archive.write(p, tag + '/' + p.relative_to(root).as_posix())
    with zipfile.ZipFile(temporary) as archive:
        require(archive.testzip() is None, 'Export archive CRC failed')
    os.replace(temporary, target)
    print('TAKE ONLY', target, flush=True)
    print('EXPORT MiB', round(target.stat().st_size / 1024**2, 1), 'SHA256', sha(target), flush=True)
    return target
