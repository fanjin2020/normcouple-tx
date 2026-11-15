"""Re-report completed evidence without changing any fit or scientific outcome."""
import shutil
import numpy as np
import pandas as pd
from scipy import stats
from next_core import require, save_json, read_json
from enhancement_report import FIELDS, CONTEXT
from report import read_nonempty as read, validate_diagnostic

NUMERICAL_FLOOR = 1e-6


def effective_counts(detail, monte_carlo):
    keys = ['unit_id', 'target_donor', 'base_subset', *CONTEXT]
    require(not detail.duplicated(keys).any() and not monte_carlo.duplicated(keys).any(), 'Duplicate donor/subset rows')
    counts = detail.merge(monte_carlo[keys + [f+'__count' for f in FIELDS]],
                          on=keys, validate='one_to_one', how='outer', indicator=True)
    require(counts['_merge'].eq('both').all(), 'Missing Monte Carlo counts')
    records = []
    for key, g in counts.groupby(CONTEXT, dropna=False):
        row = dict(zip(CONTEXT, key), planned_donors=g.target_donor.nunique(),
                   planned_subset_records=int(g.planned_records.sum()),
                   feasible_subset_records=int(g.feasible_records.sum()),
                   feasible_donors=g.loc[g.feasible_records.gt(0), 'target_donor'].nunique())
        for f in FIELDS:
            n = g[f+'__count']
            valid = n.gt(0) & np.isfinite(g[f])
            row[f+'__contributing_donors'] = int(g.loc[valid, 'target_donor'].nunique())
            row[f+'__contributing_donor_subset_pairs'] = int(valid.sum())
            row[f+'__usable_repetition_subset_pairs'] = int(n.sum())
            row[f+'__min_repetitions_per_contributing_pair'] = int(n[valid].min()) if valid.any() else 0
        records.append(row)
    return pd.DataFrame(records)


def reliability(frame):
    result = frame.copy()
    result['original_status'] = result.status
    result['original_score'] = result.score
    finite = np.isfinite(result.score)
    base_ok = result.availability.gt(0) & finite
    warn = result.get('local_smoothness_warning', pd.Series('', index=result.index)).astype('string').fillna('').str.lower()
    no_warning = warn.isin(['', 'false', '0', '0.0'])
    agreement = pd.to_numeric(result.get('derivative_agreement', np.nan), errors='coerce')
    rank = result.get('stable_nuisance_rank', pd.Series(False, index=result.index)).astype('string').fillna('false').str.lower().isin(['true', '1', '1.0'])
    status_ok = result.status.eq('LOCAL_SENSITIVITY')
    explicit_nonsmooth = result.status.eq('FINITE_DIFFERENCE_SENSITIVITY_NONSMOOTH')
    result['reliability'] = 'UNAVAILABLE'
    result.loc[status_ok & ~base_ok, 'reliability'] = 'UNAVAILABLE_BASE_FIT_OR_SCORE'
    smooth = status_ok & base_ok & rank & no_warning & np.isfinite(agreement) & (agreement <= .1)
    result.loc[(status_ok | explicit_nonsmooth) & base_ok & ~smooth, 'reliability'] = 'NONSMOOTH_OR_RANK_UNSTABLE'
    result.loc[smooth, 'reliability'] = 'LOCAL_DERIVATIVE_QUALIFIED'
    numerical = smooth & result.score.abs().le(NUMERICAL_FLOOR)
    result.loc[numerical, 'reliability'] = 'NUMERICAL_INVARIANCE'
    result['derivative_qualified'] = smooth
    result['score'] = result.score.where(smooth)
    result.loc[numerical, 'score'] = 0.
    return result


def qualified_diagnostics(frame):
    qualified = reliability(frame)
    heldout, performance = validate_diagnostic(qualified)
    planned = []
    for (method, outcome), g in qualified.groupby(['method', 'outcome']):
        planned.append(dict(method=method, outcome=outcome, planned_records=len(g),
            original_finite_scores=int(np.isfinite(g.original_score).sum()),
            qualified_scores=int(g.derivative_qualified.sum()),
            nonsmooth_records=int(g.reliability.eq('NONSMOOTH_OR_RANK_UNSTABLE').sum()),
            numerical_invariance_records=int(g.reliability.eq('NUMERICAL_INVARIANCE').sum()),
            status='QUALIFIED_SUBSET_ONLY' if g.derivative_qualified.any() else 'NO_QUALIFIED_DERIVATIVE',
            interpretation='Post-hoc reliability audit of the original held-out scores; not independent validation'))
    return heldout, performance, pd.DataFrame(planned)


def bh(p):
    p = np.asarray(p, float)
    require(np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all(), 'Invalid p-values')
    ix = np.argsort(p, kind='stable')
    q = np.minimum.accumulate((p[ix]*len(p)/np.arange(1, len(p)+1))[::-1])[::-1]
    out = np.empty(len(p))
    out[ix] = np.minimum(q, 1.)
    return out


def validate_associations(frame, expected_genes, scale=False):
    keys = ['reference', 'method'] + (['fraction_scale'] if scale else [])
    require(not frame.duplicated(keys+['gene_id']).any(), 'Duplicated clinical gene/model keys')
    expected = set(map(str, expected_genes))
    checks = []
    for key, g in frame.groupby(keys):
        require(set(g.gene_id.astype(str)) == expected, 'Clinical gene family changed')
        require(np.allclose(bh(g.p), g.q, atol=2e-12, rtol=1e-9), 'Clinical BH recomputation disagrees')
        available = g.available.astype(str).str.lower().eq('true')
        require(np.isfinite(g.loc[available, ['beta', 'se']]).all().all() and (g.loc[available, 'se'] > 0).all(),
                'Invalid available clinical estimates')
        require((g.loc[~available, 'p'] == 1).all(), 'Unavailable genes changed BH family handling')
        checks.append(dict(zip(keys, key), genes=len(g), available=int(available.sum()),
                           discoveries=int((available & g.q.le(.05)).sum()), bh_verified=True))
    require(len(checks) == (20 if scale else 10), 'Missing clinical reference/method combinations')
    return pd.DataFrame(checks)


def clinical_exports(parent, dest):
    source = parent / 'clinical'
    genes = read_json(source / 'protocol.json')['test_genes']
    require(len(genes) == 16161, 'Unexpected parent clinical family')
    association = read(source / 'disease_associations.csv.gz')
    scales = read(source / 'fraction_scale_associations.csv.gz')
    checks = validate_associations(association, genes)
    scale_checks = validate_associations(scales, genes, scale=True)
    stability = read(source / 'method_stability.csv')
    require(not stability.duplicated(['reference', 'method']).any() and len(stability) == 10,
            'Incomplete clinical method summary')
    for reference, group in association.groupby('reference'):
        base = group[group.method.eq('tmm')].set_index('gene_id')
        for method, g in group.groupby('method'):
            other = g.set_index('gene_id').reindex(base.index)
            a = base.available & base.q.le(.05)
            b = other.available & other.q.le(.05)
            row = stability[stability.reference.eq(reference) & stability.method.eq(method)].iloc[0]
            for name, value in [('discoveries', b.sum()), ('gained', (b & ~a).sum()),
                                ('lost', (a & ~b).sum()), ('threshold_crossings', (a ^ b).sum())]:
                require(int(row[name]) == int(value), 'Clinical summary differs from gene-level table: '+name)
    checks.to_csv(dest / 'clinical_gene_checks.csv', index=False)
    scale_checks.to_csv(dest / 'clinical_fraction_scale_checks.csv', index=False)
    for name in ('disease_associations.csv.gz', 'fraction_scale_associations.csv.gz',
                 'method_stability.csv', 'clinical_fraction_scale.csv', 'analysis_cohort.csv',
                 'cohort_exclusions.csv', 'design_diagnostics.csv', 'fraction_yields.csv'):
        shutil.copy2(source / name, dest / ('clinical_'+name if not name.startswith('clinical_') else name))


def make_corrections(parent, output, signature):
    from support import checkpoint, seal
    dest = output / 'corrections'
    if checkpoint(dest, signature):
        return
    dest.mkdir(exist_ok=True)
    old = parent / 'final'
    effective = effective_counts(read(old/'donor_subset_results.csv'), read(old/'within_unit_monte_carlo.csv'))
    effective.to_csv(dest/'effective_sample_counts.csv', index=False)
    original = read(old/'method_summary.csv')
    result = original.merge(effective, on=CONTEXT, validate='one_to_one')
    require(len(result) == len(original), 'Effective-count report lost settings')
    result.to_csv(dest/'method_summary_with_effective_counts.csv', index=False)
    source = read(old/'diagnostic_heldout.csv.gz')
    heldout, performance, coverage = qualified_diagnostics(source)
    heldout.to_csv(dest/'diagnostic_reliability.csv.gz', index=False, compression='gzip')
    performance.to_csv(dest/'diagnostic_qualified_performance.csv', index=False)
    coverage.to_csv(dest/'diagnostic_qualification_counts.csv', index=False)
    for name in ('diagnostic_performance.csv', 'diagnostic_setting_performance.csv', 'perturbation_support.csv'):
        shutil.copy2(old/name, dest/('original_'+name))
    clinical_exports(parent, dest)
    save_json(dest/'reporting_policy.json', dict(
        parent_fits_changed=False, primary_outcomes_changed=False,
        numerical_score_floor=NUMERICAL_FLOOR, derivative_agreement_max=.1,
        original_diagnostic_performance='Retained, including non-smooth scores; not replaced by the qualified subset',
        qualified_diagnostic_performance='Post-hoc reliability audit, not an independently validated predictor',
        averaging='Original equal-subset, then equal-target-donor weighting preserved',
        effective_counts='Per-metric contributing donors and repetition-subset pairs; not independent biological replicates',
        clinical='Fixed test family; gene tables exported from verified parent receipt; BH independently recalculated'))
    seal(dest, [p.name for p in dest.iterdir() if p.is_file() and p.name != 'complete.json'], signature)
