"""Frozen public inputs and donor-isolated preparation for the anchor benchmark."""
from pathlib import Path
import re
import tarfile

import numpy as np
import pandas as pd

from .imaging_audit import COMMIT, FILES, UPSTREAM, matched_inputs
from .io import ROOT, sha
from .download import download

DEFINITIONS = ROOT / 'publication/input_definitions'
PB_SHA = '5ea4c5ea4214fc64d43ee6f6fca3ff2bd38c203f536ba307dc8568ababdad16e'
LIBRARY = 'library_type.csv'
LIBRARY_PATH = 'processed-data/01_SPEAQeasy/' + LIBRARY
LIBRARY_SHA = '302e873fed465262b27062be954571bbdb477b7f4f356b5ba35f3dd8ce066daf'
CLASSES = ('Astro', 'EndoMural', 'Excit', 'Inhib', 'Micro', 'OPC', 'Oligo')
PUBLIC = {name: (folder + name, checksum) for name, (folder, checksum) in FILES.items()}
PUBLIC[LIBRARY] = (LIBRARY_PATH, LIBRARY_SHA)


def checked(public_root, library, allow):
    public_root = Path(public_root)
    names = {}
    for name, (relative, checksum) in PUBLIC.items():
        path = Path(library) if name == LIBRARY and library else public_root / name
        if path.is_file():
            if sha(path) != checksum:
                raise ValueError('Pinned upstream file differs: ' + str(path))
        elif allow and not (name == LIBRARY and library):
            download(f'https://raw.githubusercontent.com/{UPSTREAM}/{COMMIT}/{relative}', path, checksum)
        else:
            raise ValueError('Missing public input: ' + str(path) + '. Supply an existing cache or --allow-downloads.')
        names[name] = path
    return names


def panels():
    frame = pd.read_csv(DEFINITIONS / 'fixed_marker_panels.csv')
    marker = frame[(frame.reference == 'LIBD') & (frame.target == 'GSE174367')].copy()
    ids = pd.read_csv(DEFINITIONS / 'complete_test_gene_families.csv')
    test = ids[(ids.reference_source == 'LIBD') & (ids.role == 'test') & ids.included].gene_id.tolist()
    if (len(test) != 100 or len(set(test)) != 100 or marker.gene_id.duplicated().any()
            or set(test) & set(marker.gene_id) or set(marker.marker_class) != set(range(1, 8))):
        raise ValueError('Invalid fixed LIBD external panels')
    return test, marker


def paired_counts(paths, test, marker):
    """Transfer the observed, equal-depth assay contrast, not an estimated disease effect."""
    counts, libraries, _, _ = matched_inputs(paths['DLPFC_bulk_counts.txt'].parent)
    library = pd.read_csv(paths[LIBRARY]).rename(columns={'Dataset': 'run'})
    if library.run.duplicated().any():
        raise ValueError('The public run identifiers are not unique')
    libraries['run'] = libraries.library.str.split('_', n=1).str[0]
    libraries = libraries.merge(library, on='run', validate='many_to_one', how='left')
    if libraries[['round', 'library_type']].isna().any().any():
        raise ValueError('Missing author library labels')
    shared = counts.index.intersection(test)
    if len(shared) < 80:
        raise ValueError('Too few fixed test genes in the study counts')
    counts = counts.astype(float)
    if not np.isfinite(counts).all().all() or counts.lt(0).any().any() or counts.sum().le(0).any():
        raise ValueError('Invalid bulk counts')
    marker_in = counts.index.intersection(marker.gene_id)
    by_block, effects = [], []
    for (donor, block, run_round), rows in libraries.groupby(['donor', 'block', 'round'], sort=True):
        if rows.library_type.duplicated().any():
            raise ValueError('Ambiguous within-block and round protocol pair')
        if set(rows.library_type) != {'polyA', 'RiboZeroGold'}:
            continue
        slots = rows.set_index('library_type').library
        before = counts[slots.polyA]
        after = counts[slots.RiboZeroGold]
        b, a = before / before.sum() * 1e6, after / after.sum() * 1e6
        log_ratio = np.log((a.loc[shared] + .5) / (b.loc[shared] + .5))
        # Algebraically substitute the test rows only, after bringing the paired
        # libraries to the same depth. All markers remain at their baseline values.
        by_block.append(dict(donor=donor, block=block, round=int(run_round),
            test_genes=len(shared), marker_genes=len(marker_in),
            polyA_reads=int(before.sum()), ribozero_reads=int(after.sum()),
            test_share_polyA=b.loc[shared].sum()/1e6, test_share_ribozero=a.loc[shared].sum()/1e6,
            total_relative_change=(a.loc[shared].sum() - b.loc[shared].sum())/1e6,
            median_absolute_test_log_ratio=float(np.median(np.abs(log_ratio)))))
        effects.append(pd.Series(log_ratio, name=donor))
    table = pd.DataFrame(by_block)
    if table.empty or table.donor.nunique() < 4:
        raise ValueError('Need at least four donors with verified assay pairs')
    # Blocks and rounds are not independent donors.
    shift = pd.DataFrame(effects).groupby(level=0).median().median().reindex(test, fill_value=0).to_numpy()
    return table, shift, len(shared)


def training_inputs(path, test, marker):
    if sha(path) != PB_SHA:
        raise ValueError('Expected the checksum-pinned GSE174367 donor-by-class pseudobulk')
    with np.load(path, allow_pickle=False) as z:
        raw, genes, groups = (z[k] for k in ('counts', 'genes', 'groups'))
    g = pd.Index(genes)
    if raw.shape != (58721, 126) or not g.is_unique or not np.issubdtype(raw.dtype, np.integer) or (raw < 0).any():
        raise ValueError('Invalid pseudobulk matrix')
    if int(raw.sum()) != 593946945:
        raise ValueError('Unexpected retained UMI total')
    pairs = [tuple(str(k).split('|')) for k in groups]
    if len(set(pairs)) != len(pairs) or any(len(k) != 2 or k[1] not in CLASSES for k in pairs):
        raise ValueError('Invalid donor/class inventory')
    # Freeze the split before looking at any holdout outcomes.
    donors = sorted({k[0] for k in pairs})
    train, validation = donors[::3], [x for x in donors if x not in donors[::3]]
    if len(train) != 6 or len(validation) != 12:
        raise ValueError('The protocol requires 6 reference and 12 held-out donors')
    order = {k: i for i, k in enumerate(pairs)}
    columns = np.array([[order[(donor, c)] for c in CLASSES] for donor in donors])
    if (raw.sum(0) == 0).any():
        raise ValueError('An empty donor/class library cannot form a profile')
    p = raw[:, columns].transpose(1, 2, 0).astype(float)
    p /= p.sum(-1, keepdims=True)
    ref = p[[donors.index(k) for k in train]]
    # A fixed analysis universe keeps the fit economical. All remaining reads
    # still contribute to the observed totals through a multinomial remainder.
    abundance = ref.mean((0, 1))
    common = np.flatnonzero((ref > 0).mean((0, 1)) >= .75)
    common = common[np.lexsort((genes[common], -abundance[common]))[:2048]]
    marker_index, test_index = g.get_indexer(marker.gene_id), g.get_indexer(test)
    if (marker_index < 0).any() or (test_index < 0).any():
        raise ValueError('Frozen marker and test IDs must be present')
    universe = np.unique(np.r_[common, marker_index, test_index])
    mark = np.searchsorted(universe, marker_index)
    test_rows = np.searchsorted(universe, test_index)
    candidate = np.setdiff1d(common, np.r_[marker_index, test_index])
    score = np.log2(1 + ref[:, :, candidate] * 1e6).std((0, 1), ddof=1)
    order = np.lexsort((genes[candidate], -abundance[candidate], score))[:256]
    if len(order) != 256:
        raise ValueError('Too few independent anchors')
    selected = candidate[order]
    w = 1 / (.25 + score[order])
    w /= w.mean()
    reference = ref.mean(0)[:, universe]
    pool = np.searchsorted(universe, selected)
    profile = p[[donors.index(k) for k in validation]][:, :, universe]
    anchor = pd.DataFrame(dict(gene_id=genes[selected], weight=w,
                              reference_frequency=abundance[selected], log_cpm_sd=score[order]))
    return dict(train=train, held=validation, ref=reference, profile=profile, groups=CLASSES,
                genes=genes[universe], mark=mark, owners=marker.marker_class.to_numpy(int)-1,
                test=test_rows, anchor=pool, weights=w, panel=anchor, abundance=abundance[selected])


def clinical_count_capacity(tar, selected=None, mapping=None):
    """Original clinical family, with optional exact one-to-one anchor mapping."""
    from .inputs import specification
    from .io import read_json
    from .runtime import materialize
    import sys
    import tempfile

    tar = Path(tar)
    pinned = specification()['bulk/GSE125583_RAW.tar']['sha256']
    if sha(tar) != pinned:
        raise ValueError('The archived GSE125583 count object differs')
    source = pd.read_csv(ROOT / 'publication/data/clinical_coefficient_pairs.csv.gz', dtype={'gene_id': str})
    tests = source.groupby('reference').gene_id.agg(set)
    if len(tests) != 2 or any(t != tests.iloc[0] for t in tests) or len(tests.iloc[0]) != 16161:
        raise ValueError('Unexpected original clinical test family')
    meta = pd.read_csv(DEFINITIONS / 'clinical_bulk_cohort.csv').set_index('sample_id')
    if len(meta) != 288 or not meta.index.is_unique:
        raise ValueError('Unexpected original cohort')
    with tarfile.open(tar, 'r') as stream:
        ids = [re.fullmatch(r'(GSM\d+)_sample\d+\.tsv\.gz', name).group(1) for name in stream.getnames()]
    if len(set(ids)) != 289 or not set(meta.index) <= set(ids):
        raise ValueError('Unexpected deposited sample family')
    with tempfile.TemporaryDirectory(prefix='normcouple_') as tmp:
        runtime = materialize(Path(tmp) / 'runtime')
        sys.path.insert(0, str(runtime / 'common_support/baseline'))
        try:
            from bulk_data import read_counts
            raw, _ = read_counts(tar, ids)
        finally:
            sys.path.pop(0)
    if not set(tests.iloc[0]) <= set(raw.index):
        raise ValueError('Clinical test family is missing from the raw assay')
    frame = raw.loc[:, meta.index]
    total = frame.sum()
    non_test = frame.loc[~frame.index.isin(tests.iloc[0])]
    cap = non_test.sum() / total
    value = pd.DataFrame(dict(group=meta.group, non_test_share=cap, non_test_counts=non_test.sum(),
                              detected_non_test=(non_test > 0).sum()))
    if mapping:
        digest = read_json(DEFINITIONS/'public_input_provenance.json')['annotation_provenance']['snapshot_sha256']
        if sha(mapping) != digest:
            raise ValueError('The pinned NCBI mapping snapshot differs')
        lookup = pd.read_csv(mapping, dtype={'entrez':str})
        lookup = lookup.loc[lookup.one_to_one.astype(str).str.lower().eq('true')]
        anchors = lookup.loc[lookup.ensembl.isin(selected),'entrez']
        if anchors.duplicated().any():
            raise ValueError('Ambiguous anchor mapping')
        present = non_test.loc[non_test.index.isin(anchors)]
        value['mapped_disjoint_anchors'] = len(present)
        value['anchor_counts'] = present.sum()
        value['detected_anchors'] = (present > 0).sum()
    rows = []
    for group, subset in value.groupby('group'):
        for field in value.columns.drop('group'):
            rows.append(dict(group=group, samples=len(subset), metric=field,
                             min=subset[field].min(), p05=subset[field].quantile(.05),
                             median=subset[field].median(), max=subset[field].max()))
    return pd.DataFrame(rows)
