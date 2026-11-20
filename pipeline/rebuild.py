"""Reconstruct the pinned reference banks from author-distributed raw counts."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from .download import download
from .inputs import DEFINITIONS, specification
from .io import ROOT, atomic_json, lock, read_json, separate, sha

CLASSES = ['Astro', 'EndoMural', 'Excit', 'Inhib', 'Micro', 'OPC', 'Oligo']


def column(group, name):
    node = group[name.decode() if isinstance(name, bytes) else name]
    if isinstance(node, h5py.Dataset):
        value = node[()]
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in value]) if value.dtype.kind in 'OSU' else value
    categories, codes = column(node, 'categories'), node['codes'][()].astype(int)
    if ((codes < -1) | (codes >= len(categories))).any():
        raise ValueError('Invalid categorical codes')
    value = np.full(len(codes), '', dtype=object)
    value[codes >= 0] = categories[codes[codes >= 0]]
    return value


def seaad_donor(path, seed, mapping):
    donor = path.name.split('_SEAAD_')[0]
    with h5py.File(path, 'r') as f:
        obs, var = f['obs'], f['var']
        ids = column(obs, obs.attrs.get('_index', '_index')).astype(str)
        genes = column(var, var.attrs.get('_index', '_index')).astype(str)
        if set(column(obs, 'Donor ID')) != {donor} or set(column(obs, 'Brain Region')) != {'DFC'}:
            raise ValueError('Source donor or anatomy mismatch')
        labels = column(obs, 'Subclass').astype(str)
        keep = (np.char.lower(column(obs, 'Used in analysis').astype(str)) == 'true') & (np.char.lower(column(obs, 'Neurotypical reference').astype(str)) != 'true')
        if set(labels[keep]) - set(mapping) or len(set(ids)) != len(ids) or len(set(genes)) != len(genes):
            raise ValueError('Unmapped subclass or duplicate source identifiers')
        broad = np.array([mapping.get(c, '') for c in labels])
        pools = {c: np.flatnonzero(keep & (broad == c)) for c in CLASSES}
        if any(len(p) < 10 for p in pools.values()):
            raise ValueError('Fewer than ten eligible cells in a class')
        fraction = min(1., 2000 / max(map(len, pools.values())))
        rng = np.random.default_rng(seed)
        selected = []
        for pool in pools.values():
            n = min(len(pool), max(10, int(np.floor(len(pool) * fraction))))
            selected.extend(np.sort(rng.choice(pool, n, replace=False)))
        selected = np.array(sorted(selected))
        matrix = f['layers/UMIs']
        encoding = matrix.attrs.get('encoding-type', '')
        if encoding not in ('csr_matrix', b'csr_matrix') or tuple(matrix.attrs['shape']) != (len(ids), len(genes)):
            raise ValueError('Expected author cell-by-gene CSR UMI layer')
        values, indices, indptr = [], [], [0]
        for i in selected:
            lo, hi = int(matrix['indptr'][i]), int(matrix['indptr'][i + 1])
            v, j = matrix['data'][lo:hi], matrix['indices'][lo:hi]
            if not (np.isfinite(v).all() and (v >= 0).all() and np.equal(v, np.floor(v)).all() and ((j >= 0) & (j < len(genes))).all()):
                raise ValueError('Invalid raw UMI values or feature indices')
            values.append(v)
            indices.append(j)
            indptr.append(indptr[-1] + len(v))
        x = sparse.csr_matrix((np.concatenate(values), np.concatenate(indices), np.array(indptr)), shape=(len(selected), len(genes))).T.tocsc()
        meta = pd.DataFrame(dict(cell_id=[donor + '::' + i for i in ids[selected]], donor_id=donor,
            cell_type=broad[selected], source_id='SEAAD_DFC_2026_06_22', author_subclass=labels[selected]))
    if (np.asarray(x.sum(0)).ravel() <= 0).any():
        raise ValueError('Selected zero-library cell')
    return x, meta, genes


def execute(args):
    output = separate(args.output, args.source_dir)
    output.mkdir(parents=True, exist_ok=True)
    sources = read_json(DEFINITIONS / 'raw_reference_sources.json')
    expected = {n: v['sha256'] for n, v in specification().items() if n.startswith(('LIBD/', 'SEAAD/')) or n == 'gene_identity_map.csv'}
    with lock(output / '.rebuild.lock'):
        if any((output / n).exists() for n in expected):
            if all((output / n).is_file() and sha(output / n) == h for n, h in expected.items()):
                print('RECONSTRUCTION ALREADY VERIFIED', output)
                return
            raise ValueError('Use an empty destination; existing partial or different banks were not overwritten')
        paths = []
        args.source_dir.mkdir(parents=True, exist_ok=True)
        for item in sources['files']:
            path = args.source_dir / item['name']
            if not path.is_file() and args.allow_downloads:
                download(item['url'], path, item['sha256'])
            if not path.is_file() or sha(path) != item['sha256']:
                raise ValueError('Missing or different raw source: ' + str(path))
            paths.append(path)
            print('VERIFIED RAW', path.name, flush=True)
        with tempfile.TemporaryDirectory(prefix='.rebuild_', dir=output) as temp:
            temp = Path(temp)
            libd, seaad = temp / 'LIBD', temp / 'SEAAD'
            libd.mkdir()
            seaad.mkdir()
            with zipfile.ZipFile(paths[0]) as archive:
                for name in ('sce_DLPFC_annotated/se.rds', 'sce_DLPFC_annotated/assays.h5'):
                    # Only two exact, non-user-controlled member names are extracted.
                    with archive.open(name) as src, (temp / Path(name).name).open('wb') as dst:
                        shutil.copyfileobj(src, dst, 8 * 1024 * 1024)
            subprocess.run([args.rscript, str(ROOT / 'pipeline/rebuild_libd.R'), str(temp), str(libd),
                            *[str(p.resolve()) for p in args.r_library]], check=True)
            blocks, metadata, genes = [], [], None
            for i, path in enumerate(paths[1:]):
                x, m, g = seaad_donor(path, sources['seed'] + 35400 + i, sources['mapping'])
                if genes is not None and not np.array_equal(g, genes):
                    raise ValueError('Source feature order differs')
                genes = g
                blocks.append(x)
                metadata.append(m)
                print('RECONSTRUCTED', path.name, flush=True)
            mmwrite(seaad / 'counts.mtx', sparse.hstack(blocks, format='csc'))
            pd.concat(metadata, ignore_index=True).to_csv(seaad / 'cells.csv', index=False)
            pd.DataFrame({'gene_id': genes}).to_csv(seaad / 'genes.csv', index=False)
            shutil.move(libd / 'gene_identity_map.csv', temp / 'gene_identity_map.csv')
            actual = {name: sha(temp / name) for name in expected}
            if actual != expected:
                raise ValueError('Raw reconstruction differs from the study input hashes: ' + str({k: v for k, v in actual.items() if v != expected[k]}))
            for name in expected:
                dest = output / name
                dest.parent.mkdir(parents=True, exist_ok=True)
                (temp / name).replace(dest)
        atomic_json(output / 'reconstruction.json', dict(raw_sources=sources, output_sha256=actual,
                    status='All seven bank files match the computational specification byte for byte'))
        print('RECONSTRUCTION VERIFIED', output)
