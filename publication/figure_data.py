"""Paired quantities used in the publication figures (rates are fractions)."""
import numpy as np
import pandas as pd


def subset_points(data):
    keys = ['reference_source', 'target_source']
    metrics = ['proportion_rmse', 'fdp']
    df = data.loc[data.method.eq('marker_matched') & data.outcome.eq('null') &
                  (data.support_rule.eq('positive_only') | data.base_subset.eq('zero'))]
    if df.empty or df.duplicated(keys + ['base_subset', 'target_donor']).any():
        raise ValueError('Missing or repeated donor/subset observations')
    if not np.isfinite(df[metrics]).all().all():
        raise ValueError('Missing RMSE or FDR observation')
    for _, route in df.groupby(keys):
        rosters = route.groupby('base_subset').target_donor.agg(frozenset)
        if len(set(rosters)) != 1:
            raise ValueError('Marker subsets do not share the same target donors')
    expected = {'zero', 'legacy'} | {f'{kind}_{i:03}' for kind in ('random','matched') for i in range(1,4)}
    if any(set(route.base_subset) != expected for _, route in df.groupby(keys)):
        raise ValueError('Fig6 requires the control, predefined, three random and three matched subsets in each setting')
    means = df.groupby(keys + ['base_subset'], as_index=False)[metrics].mean()
    random = df[df.base_subset.str.startswith('random_')]
    if random.empty or set(map(tuple, random[keys].to_numpy())) != set(map(tuple, df[keys].to_numpy())):
        raise ValueError('Missing random subsets for a reference-target setting')
    per_donor = random.groupby(keys + ['target_donor'], as_index=False)[metrics].mean()
    average = per_donor.groupby(keys, as_index=False)[metrics].mean()
    average['base_subset'] = 'random_mean'
    return pd.concat([means, average], ignore_index=True)


CLASSES = ('Both', 'TMM only', 'MRN only', 'Neither')


def fusiform_doses(finite, probes):
    """Validate the matched clinical family before combining S13 and S20 in Fig9."""
    f = probes.loc[probes.cohort.eq('GSE125583')].copy()
    if f.probe.isna().any() or f.probe.duplicated().any():
        raise ValueError('Fig9: duplicated or missing clinical probe key')
    keyed = f.set_index('probe')
    ranks = ('low', 'middle', 'high')
    singles = [f'single_{rank}_{d}' for rank in ranks for d in ('up','down')]
    full = ['test_family_down','test_family_up']
    if not set(singles+full) <= set(keyed.index):
        raise ValueError('Fig9 requires single-gene and full-family probes in both directions')
    one = keyed.loc[singles]
    ids = one.selected_ids
    shares = one.AD_count_share_pct.to_numpy().reshape(3,2)
    if (not one.selected_genes.eq(1).all() or ids.isna().any() or ids.nunique()!=3 or
            not np.array_equal(ids.iloc[::2].to_numpy(),ids.iloc[1::2].to_numpy()) or
            not np.isfinite(shares).all() or not np.allclose(shares[:,0],shares[:,1]) or
            (np.diff(shares[:,0])<0).any()):
        raise ValueError('Fig9: singleton rank, identity or baseline count share differs between directions')
    for group in ('AD','control'):
        q = keyed.loc[full,[f'{group}_effective_{v}_pct' for v in ('q25','median','q75')]].to_numpy()
        if not np.isfinite(q).all() or (np.diff(q,axis=1)<0).any():
            raise ValueError('Fig9: invalid effective-denominator interval')
    up = keyed.loc['test_family_up']
    for ref in ('LIBD','SEAAD'):
        curve = finite.loc[finite.reference.eq(ref) & finite.method.eq('tmm')]
        if curve.empty or curve.log_count_step.duplicated().any():
            raise ValueError('Fig9: missing or repeated TMM dose')
        end = curve.loc[curve.log_count_step.idxmax()]
        n = up.selected_genes
        crossed = up[ref+'_family_BH_crossings']
        if (not n>0 or end.paired_available!=n or end.threshold_crossings!=crossed or
                not np.isclose(end.threshold_crossing_fraction,crossed/n,rtol=1e-10,atol=1e-13) or
                not np.isclose(up.AD_input_change_pct,100*np.expm1(end.log_count_step))):
            raise ValueError(f'Fig9: S13 and S20 disagree on the dose, family or BH crossings for {ref}')
    return f


def clinical_status(g, summary):
    if g.empty or g.gene_id.isna().any() or g.gene_id.duplicated().any():
        raise ValueError('Expected unique paired clinical genes')
    q = g[['tmm_q', 'mrn_q']]
    if not (q.ge(0) & q.le(1)).all().all():
        raise ValueError('Incomplete or invalid paired clinical q values')
    t, m = g.tmm_q.le(.05), g.mrn_q.le(.05)
    status = np.select([t & m, t, m], CLASSES[:3], default=CLASSES[3])
    counts = pd.Series(status).value_counts().reindex(CLASSES, fill_value=0)
    expect = [len(g), counts['Both'] + counts['MRN only'], counts['MRN only'], counts['TMM only']]
    if list(summary[['planned_genes', 'discoveries', 'gained', 'lost']]) != expect:
        raise ValueError('Paired BH counts disagree with the full-family stability table')
    return status, counts


def imaging_rows(frame):
    keys = ['reference', 'cell_type']
    if (frame.empty or frame.duplicated(keys).any() or
            set(frame.panel) != {'reference_transfer'} or set(frame.reference) != {'LIBD','SEAAD'}):
        raise ValueError('Expected paired reference-transfer imaging estimates')
    counts = frame.pivot(index='cell_type', columns='reference', values='donors')
    if counts.isna().any().any() or not counts.LIBD.eq(counts.SEAAD).all():
        raise ValueError('Imaging donor support differs across references')
    blocks = frame.pivot(index='cell_type', columns='reference', values='blocks')
    if blocks.isna().any().any() or not blocks.LIBD.eq(blocks.SEAAD).all() or (blocks < counts).any().any():
        raise ValueError('Imaging tissue-block support differs')
    if not np.allclose(frame.mrn_mae - frame.tmm_mae, frame.mrn_minus_tmm_mae, rtol=1e-8, atol=1e-12):
        raise ValueError('Imaging difference is not MRN MAE minus TMM MAE')
    return frame.pivot(index='cell_type', columns='reference', values='mrn_minus_tmm_mae').mul(100).join(
        counts.LIBD.rename('donors')).join(blocks.LIBD.rename('blocks'))
