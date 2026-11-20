"""Import exact cached inputs or fetch author-specified, checksum-pinned objects."""
import os
from pathlib import Path
import shutil
import tarfile
from .download import download
from .io import ROOT, atomic_json, inside, lock, read_json, separate, sha

DEFINITIONS = ROOT / 'publication/input_definitions'


def specification():
    spec = read_json(DEFINITIONS / 'computational_specification.json')
    items = {}
    for name, value in spec['input_sha256'].items():
        name = 'gene_identity_map.csv' if name == 'gene_map' else name
        url = None
        if name.startswith('external/'):
            url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174367/suppl/' + Path(name).name
        items[name] = {'sha256': value, 'url': url}
    for item in read_json(DEFINITIONS / 'public_input_provenance.json')['sources']:
        items['bulk/' + item['name']] = {'sha256': item['sha256'], 'url': item['url']}
    return items


def source_manifest(path=None):
    result = specification()
    if path:
        supplied = read_json(path)
        supplied = supplied.get('files', supplied)
        for name, item in supplied.items():
            if name not in result or item.get('sha256') != result[name]['sha256']:
                raise ValueError('Source manifest changes a fixed study input: ' + name)
            if item.get('url'):
                result[name]['url'] = item['url']
    return result


def verify_inputs(folder, bulk_cache=None, include_bulk=True):
    folder = Path(folder).expanduser().resolve()
    bulk_cache = Path(bulk_cache).expanduser().resolve() if bulk_cache else folder / 'bulk'
    records, missing = {}, []
    for name, item in specification().items():
        if name.startswith('bulk/'):
            if not include_bulk:
                continue
            path = bulk_cache / Path(name).name
        else:
            path = inside(folder, name)
        if not path.is_file():
            missing.append(name + ' (missing)')
        elif sha(path) != item['sha256']:
            missing.append(name + ' (checksum differs)')
        else:
            records[name] = item['sha256']
            print('VERIFIED INPUT', name, flush=True)
    if missing:
        raise ValueError('Inputs are not ready:\n' + '\n'.join(missing) +
                         '\nUse prepare-inputs with your existing input/cache directory. No fits were started.')
    return records


def stage_bulk_cache(folder, source, verified):
    """Adapt pinned existing objects to the kernels' receipt format without editing the source cache."""
    folder = separate(folder, source)
    folder.mkdir(parents=True, exist_ok=True)
    for name, item in specification().items():
        if not name.startswith('bulk/'):
            continue
        original = Path(source) / Path(name).name
        if verified.get(name) != item['sha256'] or sha(original) != item['sha256']:
            raise ValueError('Bulk cache changed after input verification: ' + name)
        target = folder / original.name
        if target.exists():
            if sha(target) != item['sha256']:
                raise ValueError('Runtime input differs: ' + name)
        else:
            temporary = target.with_name(target.name + '.importing')
            shutil.copyfile(original, temporary)
            if sha(temporary) != item['sha256']:
                raise ValueError('Bulk object changed during staging: ' + name)
            os.replace(temporary, target)
        atomic_json(target.with_name(target.name + '.download.json'), dict(
            url=item['url'], sha256=item['sha256'], bytes=target.stat().st_size,
            verification='Pinned checksum of an existing local object; not a new download'))
    return folder


def candidates(root, name):
    root = Path(root).expanduser().resolve()
    if root.is_file():
        if root.name == Path(name).name:
            yield root
        return
    bases = [root, root / 'inputs', root / 'prepared_inputs', root / 'external', root / '01_data',
             root / '01_data/paper3_GSE125583', root / 'downloads/GSE125583', root / 'bulk']
    if root.is_dir():
        # Search only known cache layouts, not every old analysis tree.
        bases += list(root.glob('*prepared_inputs'))
        bases += list(root.glob('*/downloads/GSE125583'))
        bases += list(root.glob('*gap*/cache'))
    seen = set()
    for base in bases:
        for path in (base / name, base / Path(name).name):
            if path not in seen:
                seen.add(path)
                if path.is_file():
                    yield path


def prepare(folder, caches=(), manifest=None, allow_downloads=False):
    folder = separate(folder)
    folder.mkdir(parents=True, exist_ok=True)
    items = source_manifest(manifest)
    unresolved = []
    with lock(folder / '.prepare.lock'):
        for name, item in items.items():
            dest = inside(folder, name)
            if dest.is_file():
                if sha(dest) != item['sha256']:
                    raise ValueError('Existing input differs; it was not overwritten: ' + str(dest))
                print('REUSE', name, flush=True)
                continue
            found = next((p for root in caches for p in candidates(root, name) if sha(p) == item['sha256']), None)
            if found:
                dest.parent.mkdir(parents=True, exist_ok=True)
                temporary = dest.with_name(dest.name + '.importing')
                if temporary.exists():
                    if sha(temporary) != item['sha256']:
                        shutil.copyfile(found, temporary)
                else:
                    try:
                        os.link(found, temporary)
                    except OSError:
                        shutil.copyfile(found, temporary)
                if sha(temporary) != item['sha256']:
                    raise ValueError('Cached input changed while importing: ' + str(found))
                os.replace(temporary, dest)
                print('IMPORTED CACHE', name, flush=True)
            elif allow_downloads and item['url']:
                download(item['url'], dest, item['sha256'])
            else:
                unresolved.append(name)
        atomic_json(folder / 'input_sources.json', {'files': items})
        if unresolved:
            raise ValueError('Still missing:\n' + '\n'.join(unresolved) +
                             '\nThe exact LIBD/SEA-AD banks and gene map need a verified cache or deposited URLs. '
                             'The framework does not replace them with a different cell selection.')
        atomic_json(folder / 'inputs_manifest.json', [dict(relative_path=n, sha256=v['sha256'],
                    bytes=inside(folder, n).stat().st_size) for n, v in items.items()])
    return folder


def pack_inputs(folder, destination):
    folder = Path(folder).expanduser().resolve()
    destination = separate(destination, folder)
    hashes = verify_inputs(folder)
    if destination.exists():
        raise ValueError('Input archive already exists: ' + str(destination))
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + '.partial')
    if temporary.exists():
        raise ValueError('An interrupted archive exists: ' + str(temporary))
    with tarfile.open(temporary, 'w:gz' if str(destination).endswith('.gz') else 'w') as archive:
        for name in sorted(hashes):
            archive.add(inside(folder, name), arcname=name, recursive=False)
    for name, expected in hashes.items():
        if sha(inside(folder, name)) != expected:
            raise ValueError('Input changed during packaging: ' + name)
    os.replace(temporary, destination)
    atomic_json(destination.with_suffix(destination.suffix + '.json'), {
        'archive': destination.name, 'sha256': sha(destination),
        'files': {name: {'sha256': value} for name, value in hashes.items()},
        'publication': 'Local package only. Dataset access permissions must be checked before public deposition.'})
    return destination


def unpack_inputs(archive_path, folder):
    folder = separate(folder)
    items = specification()
    folder.mkdir(parents=True, exist_ok=True)
    with lock(folder / '.prepare.lock'), tarfile.open(archive_path, 'r:*') as archive:
        names = set()
        for member in archive:
            if member.name in names or member.name not in items or not member.isfile():
                raise ValueError('Unexpected, duplicate or non-regular archive member: ' + member.name)
            names.add(member.name)
            dest = inside(folder, member.name)
            if dest.exists():
                if sha(dest) != items[member.name]['sha256']:
                    raise ValueError('Existing input differs: ' + str(dest))
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            temporary = dest.with_name(dest.name + '.importing')
            with archive.extractfile(member) as source, temporary.open('wb') as target:
                shutil.copyfileobj(source, target, 8 * 1024 * 1024)
            if sha(temporary) != items[member.name]['sha256']:
                raise ValueError('Input archive checksum differs: ' + member.name)
            os.replace(temporary, dest)
        if names != set(items):
            raise ValueError('Input archive is incomplete: ' + ', '.join(sorted(set(items) - names)))
    return folder
