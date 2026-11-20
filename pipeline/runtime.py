"""One physical copy of each kernel, with isolated import layouts at runtime."""
import shutil
from pathlib import Path
from .io import ROOT, atomic_json, digest, inside, read_json, separate, sha


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
