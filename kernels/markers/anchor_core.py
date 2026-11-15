"""Donor-level preparation and inference for the fixed-panel empirical anchor."""
import hashlib
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats

CLASSES = ['Astro', 'EndoMural', 'Excit', 'Inhib', 'Micro', 'OPC', 'Oligo']
CELL_MAP = {'ASC': 'Astro', 'PER.END': 'EndoMural', 'EX': 'Excit',
            'INH': 'Inhib', 'MG': 'Micro', 'OPC': 'OPC', 'ODC': 'Oligo'}
META = 'GSE174367_snRNA-seq_cell_meta.csv.gz'
COUNTS = 'GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5'
HASHES = {META: 'ab1a029deb43196e2bb1fea5907d838885750cfff3850c600c07588c7c7cdb2b',
          COUNTS: '6ba98a1af8772af08c8cdd5e5e63eebb0890daed48cac2e1b8961dcc67069b77'}
BASE = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174367/suppl/'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.' + str(os.getpid()) + '.partial')
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False) + '\n', encoding='utf-8')
    os.replace(tmp, path)


def metadata(path):
    d = pd.read_csv(path, dtype={'Barcode': str, 'SampleID': str})
    fields = ['Diagnosis', 'Batch', 'Age', 'Sex', 'PMI', 'RIN']
    require({'Barcode', 'SampleID', 'Cell.Type', *fields} <= set(d), 'Metadata columns missing')
    require(d.Barcode.notna().all() and d.Barcode.is_unique, 'Duplicate or missing cell barcode')
    require(d.SampleID.notna().all(), 'Missing SampleID')
    require(set(d['Cell.Type']) <= set(CELL_MAP), 'Unmapped deposited cell type')
    require(set(d.Diagnosis) == {'AD', 'Control'}, 'Unexpected diagnosis coding')
    for sid, f in d.groupby('SampleID'):
        for col in fields:
            require(f[col].nunique(dropna=False) == 1, 'Inconsistent donor metadata: ' + sid + '/' + col)
    d['cell_type'] = d['Cell.Type'].map(CELL_MAP)
    donors = d[['SampleID', *fields]].drop_duplicates().sort_values('SampleID').reset_index(drop=True)
    n = pd.crosstab(d.SampleID, d.cell_type).reindex(index=donors.SampleID, columns=CLASSES, fill_value=0)
    return d, donors, n


def design(donors, nuclei, min_cells, model):
    fields = ['Diagnosis', 'Age', 'Sex', 'Batch']
    if model == 'extended':
        fields += ['PMI', 'RIN']
    d = donors.copy()
    d['nuclei'] = d.SampleID.map(nuclei)
    d['included'] = d[fields].notna().all(axis=1) & d.nuclei.ge(min_cells)
    d['reason'] = ['' if ok else 'too_few_nuclei_or_missing_declared_covariate' for ok in d.included]
    m = d.loc[d.included].copy()
    if len(m) == 0:
        return m, None, d, {'status': 'UNAVAILABLE', 'reason': 'No complete samples'}
    x = pd.DataFrame({'Intercept': np.ones(len(m)), 'AD': m.Diagnosis.eq('AD').to_numpy(float)}, index=m.index)
    for col in ('Sex', 'Batch'):
        dummy = pd.get_dummies(m[col].astype(str), prefix=col, drop_first=True, dtype=float)
        x = pd.concat([x, dummy], axis=1)
    for col in ('Age', *(['PMI', 'RIN'] if model == 'extended' else [])):
        v = pd.to_numeric(m[col], errors='raise')
        # A constant nuisance variable adds no column space and is explicitly recorded.
        if v.std(ddof=0) > 0:
            x[col] = (v - v.mean()) / v.std(ddof=0)
    a = x.to_numpy(float)
    require(np.isfinite(a).all(), 'Nonfinite design')
    rank = int(np.linalg.matrix_rank(a))
    df = len(m) - rank
    counts = m.Diagnosis.value_counts().to_dict()
    info = dict(n=len(m), groups=counts, rank=rank, columns=list(x), residual_df=df,
                model=model, constant_covariates=[c for c in fields if m[c].nunique() == 1])
    if min(counts.get('AD', 0), counts.get('Control', 0)) < 4 or rank != a.shape[1] or df < 4:
        info.update(status='UNAVAILABLE', reason='Require >=4 donors/group, full rank and >=4 residual df')
    else:
        nuisance = np.delete(a, 1, axis=1)
        z = a[:, 1] - nuisance @ np.linalg.lstsq(nuisance, a[:, 1], rcond=None)[0]
        zz = float(z @ z)
        h = np.sum(a * np.linalg.pinv(a).T, axis=1)
        info.update(status='READY', exposure_information=zz,
                    exposure_r2=1 - zz / float(np.sum((a[:, 1] - a[:, 1].mean())**2)),
                    max_leverage=float(h.max()))
    return m, x, d, info


def pseudobulk(path, cells):
    """Exact barcode matching; all annotated nuclei, no all-class donor exclusion."""
    with h5py.File(path, 'r') as f:
        g = f['matrix']
        barcodes = pd.Index(g['barcodes'].asstr()[:])
        genes = pd.Index(g['features']['id'].asstr()[:]).str.replace(r'\.\d+$', '', regex=True)
        symbols = g['features']['name'].asstr()[:] if 'name' in g['features'] else genes.to_numpy()
        if 'feature_type' in g['features']:
            require(set(g['features']['feature_type'].asstr()[:]) == {'Gene Expression'}, 'Non-RNA features')
        require(barcodes.is_unique and genes.is_unique, 'Ambiguous H5 axes')
        pos = barcodes.get_indexer(cells.Barcode)
        require((pos >= 0).all(), 'Metadata barcode missing from H5')
        a = sparse.csc_matrix((g['data'][:], g['indices'][:], g['indptr'][:]), shape=tuple(g['shape'][:]))
    require(a.shape == (len(genes), len(barcodes)), 'Matrix shape disagrees with axes')
    require(np.isfinite(a.data).all() and (a.data >= 0).all() and (a.data == np.floor(a.data)).all(),
            'Require deposited nonnegative integer counts, not transformed expression')
    require(a.data.max(initial=0) < 2**53, 'Count exceeds exact integer range')
    a = a[:, pos].astype(np.int64)
    keys = cells.SampleID + '|' + cells.cell_type
    group_names = pd.Index(sorted(keys.unique()))
    incidence = sparse.csc_matrix((np.ones(len(cells), dtype=np.int64),
                  (np.arange(len(cells)), group_names.get_indexer(keys))), shape=(len(cells), len(group_names)))
    out = (a @ incidence).toarray()
    require(int(out.sum()) == int(a.sum()), 'Pseudobulk conservation check failed')
    return genes, np.asarray(symbols), out, group_names, dict(matrix_cells=len(barcodes),
        annotated_cells=len(cells), unused_unannotated_cells=len(barcodes) - len(cells),
        genes=len(genes), pseudobulks=len(group_names), retained_integer_umi=int(out.sum()))


def hc3(y, x):
    y, x = np.asarray(y, float), np.asarray(x, float)
    require(y.ndim == 1 and len(y) == len(x) and np.isfinite(y).all(), 'Invalid response')
    rank = np.linalg.matrix_rank(x)
    if rank != x.shape[1] or len(y) - rank < 2:
        return dict(available=False, reason='rank_or_df')
    pinv = np.linalg.pinv(x)
    beta = pinv @ y
    h = np.sum(x * pinv.T, axis=1)
    if h.max() >= 1 - 1e-10:
        return dict(available=False, reason='unit_leverage')
    resid = y - x @ beta
    se = np.sqrt(np.sum((pinv[1] * resid / (1 - h))**2))
    if not np.isfinite(se) or se <= 0:
        return dict(available=False, reason='nonpositive_se')
    df = len(y) - rank
    critical = stats.t.ppf(.975, df)
    return dict(available=True, beta=float(beta[1]), se=float(se), df=int(df),
                lower=float(beta[1] - critical * se), upper=float(beta[1] + critical * se),
                p=float(2 * stats.t.sf(abs(beta[1] / se), df)))


def bh(values):
    p = np.asarray(values, float)
    q = np.full(len(p), np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    if not len(ok):
        return q
    order = ok[np.argsort(p[ok], kind='stable')]
    q[order] = np.minimum(1, np.minimum.accumulate((p[order] * len(order) / np.arange(1, len(order)+1))[::-1])[::-1])
    return q
