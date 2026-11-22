"""Runtime checks and isolated import layouts for the scientific kernels."""
from importlib import import_module
import platform
import shutil
import sys
from pathlib import Path
from .io import ROOT, atomic_json, digest, inside, read_json, separate, sha


def check_enhancement_environment():
    """Require the core versions recorded by the original server runs."""
    lock = ROOT / 'requirements.txt'
    pins = dict(line.split('==', 1) for line in map(str.strip, lock.read_text(encoding='utf-8').splitlines())
                if line and not line.startswith('#'))
    if set(pins) != {'numpy', 'pandas', 'scipy', 'h5py', 'matplotlib'}:
        raise ValueError('Invalid enhancement environment lock: ' + str(lock))
    expected = {'python': '3.12.13', **pins}
    actual = {'python': platform.python_version()}
    for name in pins:
        try:
            actual[name] = import_module(name).__version__
        except (ImportError, OSError, AttributeError) as error:
            actual[name] = str(error)
    differences = [f'{name}: expected {version}, found {actual[name]}'
                   for name, version in expected.items() if actual[name] != version]
    if differences:
        raise ValueError('Use the original server enhancement environment (PAPER3_PYTHON). '
                         'No packages were installed.\n' + sys.executable + '\n' + '\n'.join(differences))
    return actual


def materialize(destination):
    destination = separate(destination)
    spec = read_json(ROOT / 'kernels/layout.json')
    identity = digest(spec)
    destination.mkdir(parents=True, exist_ok=True)
    marker = destination / '.kernel-layout.json'
    if marker.exists() and read_json(marker).get('signature') != identity:
        raise ValueError('Runtime belongs to a different kernel layout; use a new run directory')
    for relative, entry in spec['files'].items():
        source = inside(ROOT, entry['source'])
        if sha(source) != entry['sha256']:
            raise ValueError('Kernel checksum differs: ' + entry['source'])
        target = inside(destination, relative)
        if target.exists():
            if sha(target) != entry['sha256']:
                raise ValueError('Runtime kernel was modified: ' + relative)
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
    atomic_json(marker, {'signature': identity, 'files': len(spec['files'])})
    return destination
