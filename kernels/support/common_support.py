"""Equal-budget follow-up on the intersection of positive sample/class support."""
import hashlib
import json
from pathlib import Path
import signal
import tempfile
import time
import numpy as np
import pandas as pd
from next_core import BASE, require, read_json
import revision_core as C
import revision_driver as D
from extension_driver import marker_properties
from extension_stats import evaluate
from experiment import METHODS as DT_METHODS, read_predictions, exclusion_anchor
from enhancement_experiment import METHODS, nnls_prediction
from matched_markers import allocate, variants, split_variant
from reference_allocation import frozen_weights
from support import checkpoint, seal, atomic_json
import selection

JOBS, OPTIONS = [], {}
CASES = [('null', 0.), ('up025', .25)]


def intersection_budget(observed, owners, weights, subsets, requested):
    observed, owners, weights = np.asarray(observed), np.asarray(owners), np.asarray(weights)
    requested = np.asarray(requested)
    require(requested.shape == (len(observed), 7) and (requested >= 0).all() and
            np.array_equal(requested, np.rint(requested)), 'Invalid requested budgets')
    require(np.isfinite(observed).all() and (observed >= 0).all(), 'Invalid marker counts')
    require(np.isfinite(weights).all() and (weights >= 0).all(), 'Invalid reference weights')
    mask = np.ones(requested.shape, bool)
    for name, subset in subsets.items():
        if name == 'zero':
            continue
        ix = np.asarray(subset, int)
        require(len(ix) > 0 and len(np.unique(ix)) == len(ix), 'Empty or duplicated selected markers')
        for k in range(1, 8):
            eligible = ix[(owners[ix] == k) & (weights[ix] > 0)]
            require(len(eligible) > 0, 'No positive reference weight for '+name+' class '+str(k))
            mask[:, k-1] &= (observed[:, eligible] > 0).any(axis=1)
    return requested.astype(np.int64) * mask, mask


def make_perturbations(job, base):
    observed = base['draw']['counts'][:, job['panels']['balanced']]
    owners = C.E.official_markers(job).cell_type.to_numpy(int)
    weights = frozen_weights(job, owners, 'reference')
    legacy, _ = selection.stress(job, base, 'legacy')
    requested = np.column_stack([legacy[:, owners == k].sum(1) for k in range(1, 8)])
    common, mask = intersection_budget(observed, owners, weights, job['stress_subsets'], requested)
    group = base['group'].astype(bool)
    require(not common[~group].any(), 'Control samples received a budget')
    by_class = []
    for k in range(7):
        den = int(requested[:, k].sum())
        num = int(common[:, k].sum())
        by_class.append(dict(cell_class=k+1, requested_counts=den, retained_counts=num,
            retained_fraction=num/den if den else np.nan,
            requested_sample_classes=int((requested[:, k] > 0).sum()),
            retained_sample_classes=int((common[:, k] > 0).sum()),
            common_positive_support_samples=int(mask[group, k].sum()), case_samples=int(group.sum())))
    additions, records = {}, []
    budget_hash = hashlib.sha256(common.astype('<i8').tobytes()).hexdigest()
    for variant in variants(job['stress_subsets']):
        subset, rule = split_variant(variant)
        selected = np.asarray(job['stress_subsets'][subset], int)
        if subset == 'zero':
            value = np.zeros_like(observed, dtype=np.int64)
            used = np.zeros_like(common)
        else:
            seed = int(C.M.digest(['p3-common-support-v1', base['key'], subset])[:16], 16)
            value, unsupported = allocate(weights, selected, owners, common, observed, rule, seed)
            require(value is not None and not unsupported, 'Common-support allocation should be feasible')
            used = common
        actual = np.column_stack([value[:, owners == k].sum(1) for k in range(1, 8)])
        require(np.array_equal(actual, used), 'Unequal sample/class budgets')
        require(rule != 'positive_only' or not value[observed == 0].any(), 'Positive-support violation')
        require(not value[~group].any(), 'Control samples received additions')
        requested_total, retained_total = int(requested.sum()), int(used.sum())
        selected_mass = float(observed[group][:, selected].sum())
        records.append(dict(subset_id=variant, base_subset=subset, subset_kind=subset.split('_')[0],
            support_rule=rule, intervention_feasible=True, requested_budget_total=requested_total,
            common_budget_total=int(common.sum()), applied_budget_total=retained_total,
            retained_budget_fraction=float(common.sum()/requested_total) if requested_total else np.nan,
            no_requested_budget=requested_total == 0, zero_common_budget=int(common.sum()) == 0,
            dropped_sample_classes=int(((requested > 0) & ~mask).sum()),
            common_budget_sha256=budget_hash, selected_markers=len(selected),
            selected_baseline_case_mass=selected_mass,
            effective_relative_amplitude=retained_total/selected_mass if selected_mass else np.nan,
            added_on_observed_zero=int(value[observed == 0].sum()),
            allocation_seed='shared across support rules within subset; no outcome access'))
        additions[variant] = value
    return additions, records, by_class


def folder_for(output, task):
    unit, rep = task
    return Path(output)/'support_batches'/unit['unit_id']/('r%05d' % rep)


def stop_worker(signum, frame):
    raise KeyboardInterrupt('Worker interrupted; unsealed work will be recomputed')


def initialize_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGTERM, stop_worker)


def run_batch(task):
    unit, rep = task
    opt, job = OPTIONS, JOBS[unit['job_index']]
    folder = folder_for(opt['output'], task)
    previous = checkpoint(folder, opt['signature'])
    if previous:
        return dict(previous, resumed=True)
    folder.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    meta = {k: unit[k] for k in ('unit_id', 'reference_source', 'target_source', 'target_donor')}
    base = C.draw_base(job, unit, rep, opt)
    A, my, anchor = marker_properties(job, base)
    fy = base['ref']['pooled_yield']
    additions, stress_records, class_records = make_perturbations(job, base)
    axes = read_json(Path(job['reference_folder'])/'axes.json')
    marker_ix, test_ix = np.array(axes['marker_rows'])-1, np.array(axes['test_rows'])-1
    metrics = []
    data_by_case = {name: C.intervene(job, base, dict(name=name, effect=effect, pattern='up', panel='frozen'), opt)
                    for name, effect in CASES}
    with tempfile.TemporaryDirectory(prefix='work_', dir=folder) as tmp:
        work = Path(tmp)
        C.write_bin(work/'A.bin', A)
        entries, payload, expected = [], {}, set()
        for variant, added in additions.items():
            subset, rule = split_variant(variant)
            for name, effect in CASES:
                data = data_by_case[name]
                raw = base['draw']['counts'][:, job['norm_genes']].T.copy()
                raw[marker_ix] += added.T
                marker = raw[marker_ix].T.copy()
                raw[test_ix] += data['additions'].T
                key = variant+'_'+name
                C.write_bin(work/(key+'.bin'), raw)
                C.write_bin(work/(key+'_total.bin'), data['total']+added.sum(1))
                entries.append(dict(key=key, n=unit['n'], raw=key+'.bin', A='A.bin', total=key+'_total.bin',
                    group=base['group'].tolist(), marker_anchor=anchor,
                    non_test_anchor=exclusion_anchor(job, base, A), steps=[0.]))
                identity = dict(meta, rep=rep, subset_id=variant, base_subset=subset,
                    subset_kind=subset.split('_')[0], support_rule=rule, outcome=name, effect=effect,
                    common_budget_total=int(stress_records[0]['common_budget_total']))
                payload[key] = (identity, data, marker)
                expected.update((key, m, 0.) for m in DT_METHODS)
        C.write_json(work/'jobs.json', dict(axes=axes, entries=entries))
        r = opt['r_info']
        D.run_r([r['executable'], '--vanilla', str(BASE/'bridge.R'), str(work), r['library']],
                r['library'], folder/'bridge.log', opt['r_timeout'])
        predictions, r_records = read_predictions(work, expected)
        for key, (identity, data, marker) in payload.items():
            for method in ['oracle', *METHODS]:
                if method == 'oracle':
                    props, reason = data['p'], ''
                elif method in DT_METHODS:
                    q, _, reason, _ = predictions[key, method, 0.]
                    props = C.E.convert_available_rows(q, my if method == 'marker_matched' else fy)
                else:
                    denom = predictions[key, 'tmm', 0.][1] if method == 'nnls_linear' else marker.sum(1)
                    reference = A if method == 'nnls_linear' else A/A.sum(0)
                    q, reason = nnls_prediction(reference, marker, denom)
                    props = C.E.convert_available_rows(q, fy if method == 'nnls_linear' else my)
                row, arrays, fit = evaluate(data, props, 'full_six', reason)
                valid = arrays['available'].astype(bool)
                hit = valid & (arrays['q'] <= .05)
                null = data['truth'] == 0
                row.update(intervention_feasible=True, fit_status='ESTIMABLE' if fit['basic'] else 'UNAVAILABLE',
                    planned_tests=len(null), valid_tests=int(valid.sum()), discoveries=int(hit.sum()),
                    false_discoveries=int((hit & null).sum()), true_discoveries=int((hit & ~null).sum()),
                    coefficient_error_rms=float(np.sqrt(np.mean((arrays['beta'][valid]-data['truth'][valid])**2)))
                    if valid.any() else np.nan)
                metrics.append(dict(identity, method=method, **row))
        r_records.to_csv(folder/'r_diagnostics.csv.gz', index=False, compression='gzip')
    table = pd.DataFrame(metrics)
    require(len(table) == len(additions)*len(CASES)*(len(METHODS)+1) and
            not table.duplicated(['subset_id', 'outcome', 'method']).any(), 'Incomplete follow-up batch')
    table.to_csv(folder/'metrics.csv.gz', index=False, compression='gzip')
    pd.DataFrame([dict(meta, rep=rep, **v) for v in stress_records]).to_csv(folder/'perturbations.csv', index=False)
    pd.DataFrame([dict(meta, rep=rep, **v) for v in class_records]).to_csv(folder/'budget_by_class.csv', index=False)
    record = dict(unit=unit['unit_id'], rep=rep, rows=len(table), elapsed_seconds=time.monotonic()-start,
                  finished_unix=time.time(), common_budget=int(stress_records[0]['common_budget_total']))
    atomic_json(folder/'batch.json', record)
    seal(folder, ['metrics.csv.gz', 'perturbations.csv', 'budget_by_class.csv', 'r_diagnostics.csv.gz',
                  'bridge.log', 'batch.json'], opt['signature'], **record)
    return dict(record, resumed=False)


def summarize(output, cfg):
    from finish import tasks
    from enhancement_report import donor_tables, direct_controls, subset_comparisons, CONTEXT, FIELDS
    from next_report import effective_counts
    dest = output/'common_support'
    if checkpoint(dest, cfg['signature']):
        return
    dest.mkdir(exist_ok=True)
    folders = [folder_for(output, t) for t in tasks(cfg)]
    require(all(checkpoint(p, cfg['signature']) for p in folders), 'Follow-up has incomplete batches')
    def combine(name):
        return pd.concat([pd.read_csv(p/name, keep_default_na=False, na_values=['']) for p in folders], ignore_index=True)
    metrics = combine('metrics.csv.gz')
    expected = len(folders)*len(cfg['settings']['subset_ids'])*len(CASES)*(len(METHODS)+1)
    require(len(metrics) == expected and not metrics.duplicated(['unit_id', 'rep', 'subset_id', 'outcome', 'method']).any(),
            'Incomplete common-support result grid')
    require(metrics.intervention_feasible.all(), 'Common-support design unexpectedly infeasible')
    summary = donor_tables(metrics, dest)
    direct_controls(metrics, dest)
    # Average dependent subsets inside each repetition before computing Monte Carlo SE.
    unit_rep = metrics.groupby(['unit_id', 'target_donor', *CONTEXT, 'rep'], dropna=False)[FIELDS].mean().reset_index()
    agg = unit_rep.groupby(['unit_id', 'target_donor', *CONTEXT], dropna=False)[FIELDS].agg(['mean', 'sem', 'count'])
    agg.columns = ['__'.join(c) for c in agg.columns]
    agg.reset_index().to_csv(dest/'unit_monte_carlo.csv', index=False)
    counts = effective_counts(pd.read_csv(dest/'donor_subset_results.csv', keep_default_na=False, na_values=['']),
                              pd.read_csv(dest/'within_unit_monte_carlo.csv', keep_default_na=False, na_values=['']))
    summary.merge(counts, on=CONTEXT, validate='one_to_one').to_csv(dest/'method_summary.csv', index=False)
    comparisons = subset_comparisons(metrics)
    keys = ['unit_id', 'target_donor', 'reference_source', 'target_source', 'comparison', 'base_subset',
            'subset_kind', 'support_rule', 'outcome', 'method']
    columns = [f+'_difference' for f in FIELDS]
    paired = comparisons.groupby(keys, dropna=False)[columns].agg(['mean', 'sem', 'count'])
    paired.columns = ['__'.join(c) for c in paired.columns]
    paired.reset_index().to_csv(dest/'paired_subset_differences.csv', index=False)
    budgets = combine('budget_by_class.csv')
    perturb = combine('perturbations.csv')
    budgets.to_csv(dest/'budget_by_class.csv.gz', index=False, compression='gzip')
    perturb.to_csv(dest/'perturbations.csv.gz', index=False, compression='gzip')
    mask = perturb.subset_id.eq('zero')
    per_draw = perturb.loc[mask, ['unit_id', 'reference_source', 'target_source', 'target_donor', 'rep',
        'requested_budget_total', 'common_budget_total', 'retained_budget_fraction', 'no_requested_budget',
        'zero_common_budget', 'dropped_sample_classes']]
    dose = per_draw.groupby(['unit_id', 'reference_source', 'target_source', 'target_donor']).agg(
        draws=('rep', 'size'), no_requested_budget_draws=('no_requested_budget', 'sum'),
        zero_common_budget_draws=('zero_common_budget', 'sum'),
        original_counts=('requested_budget_total', 'sum'), retained_counts=('common_budget_total', 'sum'),
        mean_draw_retained_fraction=('retained_budget_fraction', 'mean'),
        mean_dropped_sample_classes=('dropped_sample_classes', 'mean')).reset_index()
    dose['pooled_retained_fraction'] = dose.retained_counts/dose.original_counts.replace(0, np.nan)
    dose.to_csv(dest/'dose_retention.csv', index=False)
    # Include zero-budget draws in the primary estimand, and label a secondary positive-dose stratum.
    positive = metrics[metrics.common_budget_total.gt(0)]
    pos = positive.groupby(['unit_id', 'target_donor', 'base_subset', *CONTEXT], dropna=False)[FIELDS].mean()
    pos.reset_index().to_csv(dest/'positive_dose_only_secondary.csv', index=False)
    metrics.to_csv(dest/'followup_metrics.csv.gz', index=False, compression='gzip')
    atomic_json(dest/'scope.json', dict(
        requested_batches=len(folders), completed_batches=len(folders), rows=len(metrics),
        rep_start=cfg['settings']['rep_start'], repetitions=cfg['settings']['repetitions'],
        primary='All new draws, including draws with zero retained budget',
        estimand='Common-positive-support, reduced equal-budget intervention; not the original intervention',
        support='Same sample/class mask across all seven prespecified sets and both allocation rules',
        randomization='New repetition IDs; same donor resources and generators; no independent biological cohort',
        uncertainty='Within-unit Monte Carlo SE; target donors reused across reference analyses',
        selection='Fixed reference-selected sets; support uses unperturbed marker counts, never fitted outcomes'))
    seal(dest, [p.name for p in dest.iterdir() if p.is_file() and p.name != 'complete.json'], cfg['signature'])
