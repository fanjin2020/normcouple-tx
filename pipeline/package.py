"""Verify and package the same source tree that is tested locally."""
import csv
import hashlib
import os
from pathlib import Path
import zipfile
from .io import ROOT, inside, separate, sha


def verify():
    with (ROOT / 'MANIFEST_SHA256.csv').open(encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream))
    names = set()
    for item in rows:
        if item['file'] in names:
            raise ValueError('Duplicate manifest entry: ' + item['file'])
        names.add(item['file'])
        p = inside(ROOT, item['file'])
        if not p.is_file() or sha(p) != item['sha256']:
            raise ValueError('Package checksum differs: ' + item['file'])
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and
              not any(x in ('.git', '__pycache__', '.pytest_cache') for x in p.relative_to(ROOT).parts)}
    extra = actual - names - {'MANIFEST_SHA256.csv'}
    if extra:
        raise ValueError('Unregistered files in upload folder: ' + ', '.join(sorted(extra)[:10]))
    print('Verified package files:', len(rows))
    return sorted(names | {'MANIFEST_SHA256.csv'})


def archive(destination):
    names = verify()
    manifest = ROOT / 'MANIFEST_SHA256.csv'
    manifest_hash = sha(manifest)
    with manifest.open(encoding='utf-8') as stream:
        expected = {r['file']:r['sha256'] for r in csv.DictReader(stream)}
    expected['MANIFEST_SHA256.csv'] = manifest_hash
    destination = separate(destination)
    if destination.exists():
        raise ValueError('Archive already exists: ' + str(destination))
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + '.partial')
    if partial.exists():
        raise ValueError('An interrupted archive exists: ' + str(partial))
    with zipfile.ZipFile(partial, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for name in names:
            info = zipfile.ZipInfo('paper3-normalization-coupling/' + name, (2026, 11, 15, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = (0o100755 if name.endswith('.sh') else 0o100644) << 16
            payload = inside(ROOT, name).read_bytes()
            if hashlib.sha256(payload).hexdigest() != expected.get(name):
                raise ValueError('Source changed during packaging: ' + name)
            z.writestr(info, payload)
    with zipfile.ZipFile(partial) as z:
        if z.testzip() is not None:
            raise ValueError('Archive CRC check failed')
    os.replace(partial, destination)
    print('PACKAGE', destination)
    print('SHA256', sha(destination))
    return destination
