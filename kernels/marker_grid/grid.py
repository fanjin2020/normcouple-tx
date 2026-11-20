"""Prespecified marker-stress grid using the unchanged Paper 3 generator."""
import json
import os
import shutil
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

import revision_core as C
import revision_driver as D
from extension_driver import marker_properties
from extension_stats import evaluate
from mrn_tools import actual_span_hc3
from theory import coefficient_bound, log_input_bound

HERE = Path(__file__).resolve().parent
JOBS, OPTIONS = [], {}
METHODS = ['oracle', 'observed', 'tmm', 'marker_matched', 'marker_denominator_only']


def stress_grid():
    points = [dict(grid_id='zero', fraction=0., amplitude=0., selection='gene_id')]
    for fraction in (.05, .10, .20):
        for amplitude in (.25, .5, 1.):
            points.append(dict(grid_id=f'f{int(fraction*100):02d}_a{int(amplitude*100):03d}',
                               fraction=fraction, amplitude=amplitude, selection='gene_id'))
    points.append(dict(grid_id='f10_a050_random', fraction=.1, amplitude=.5, selection='fixed_random'))
    return points


def selected_rows(job, fraction, selection):
    marker = job['panels']['balanced']
    owners = C.E.official_markers(job).cell_type.to_numpy(int)
    genes = job['source']['genes'][marker]
    selected = []
    if fraction == 0:
        return np.array([], dtype=int)
    for k in range(1, 8):
        ix = np.flatnonzero(owners == k)
        C.require(len(ix) > 0, 'Missing marker class')
        ordered = ix[np.argsort(genes[ix])]
        if selection == 'fixed_random':
            ordered = C.M.rng_for(2026101401, 'grid_marker_order', job['source']['label'],
                                  genes[ordered].tolist()).permutation(ordered)
        selected.extend(ordered[:max(1, int(np.ceil(len(ix) * fraction)))] )
    return np.asarray(selected, dtype=int)


def stress_draw(job, base, point):
    selected = selected_rows(job, point['fraction'], point['selection'])
    marker = job['panels']['balanced']
    expectation = np.zeros((len(base['group']), len(marker)))
    expectation[:, selected] = point['amplitude'] * base['draw']['counts'][:, marker[selected]] * base['group'][:, None]
    # Preserve the reviewed central condition exactly, including its RNG stream.
    if point['selection'] == 'gene_id':
        rng = C.M.rng_for(2026101303, 'marker_stress', base['key'])
    else:
        rng = C.M.rng_for(2026101401, 'grid_marker_stress', base['key'], point['selection'])
    return rng.poisson(expectation), selected


def cases_for(job):
    cases = []
    for point in stress_grid():
        for outcome in C.conditions(True):
            cases.append(dict(point, outcome=outcome['name'], effect=outcome['effect'],
                              case_id=point['grid_id'] + '_' + outcome['name']))
    if 'target_source' not in job:
        cases.append(dict(stress_grid()[0], outcome='up025', effect=.25, case_id='zero_up025'))
    return cases


def run_batch(unit, start):
    job, options = JOBS[unit['job_index']], OPTIONS
    folder = Path(options['output']) / 'units' / unit['unit_id'] / f'chunk_{start:05d}'
    if D.valid_checkpoint(folder, options['signature']):
        print('RESUME', unit['unit_id'], start, flush=True)
        return dict(resumed=True, unit=unit['unit_id'], start=start)
    folder.mkdir(parents=True, exist_ok=True)
    clock = time.monotonic()
    reps = list(range(start, min(start + options['batch_size'], options['repetitions'])))
    print('START', unit['unit_id'], reps[0], reps[-1], flush=True)
    records, genes, theory, invariants, bridge_rows = [], [], [], [], []
    with tempfile.TemporaryDirectory(prefix='working_', dir=folder) as td:
        work = Path(td)
        entries, payload, meta = [], {}, {}
        for rep in reps:
            base = C.draw_base(job, unit, rep, options)
            A, marker_yield, anchor = marker_properties(job, base)
            key = f'r{rep}'
            C.write_bin(work / (key + '_raw.bin'), base['draw']['counts'][:, job['norm_genes']].T)
            C.write_bin(work / (key + '_A.bin'), A)
            original = base['draw']['mass'].sum(1)
            C.write_bin(work / (key + '_total.bin'), original)
            stress_by_point = {}
            for point in stress_grid():
                stress, selected = stress_draw(job, base, point)
                stress_by_point[point['grid_id']] = (stress, selected)
                C.write_bin(work / (key + '_' + point['grid_id'] + '_stress.bin'), stress.T)
            outcomes = {}
            for case in cases_for(job):
                name = case['outcome']
                if name not in outcomes:
                    condition = dict(name=name, effect=case['effect'], pattern='up', panel='frozen')
                    outcomes[name] = C.intervene(job, base, condition, options)
                    C.write_bin(work / (key + '_' + name + '_add.bin'), outcomes[name]['additions'].T)
                    if rep == 0 and name in ('null', 'up05'):
                        C.verify_baseline(job, unit, rep, options, outcomes[name], case['effect'])
            conditions = []
            marker_raw = base['draw']['counts'][:, job['panels']['balanced']]
            for case in cases_for(job):
                stress, selected = stress_by_point[case['grid_id']]
                data = outcomes[case['outcome']]
                payload[rep, case['case_id']] = data
                conditions.append(dict(name=case['case_id'], additions=key + '_' + case['outcome'] + '_add.bin',
                                       stress=key + '_' + case['grid_id'] + '_stress.bin'))
                total = original + stress.sum(1)
                check = log_input_bound(marker_raw + stress, total, total + data['additions'].sum(1))
                C.require(check['max_violation'] < 1e-8, 'Finite log-input bound failed')
                theory.append(dict(rep=rep, case_id=case['case_id'], kind='log_input',
                                   eligible=True, max_violation=check['max_violation'],
                                   actual_max=float(check['actual'].max()), bound_max=float(check['bound'].max())))
                meta[rep, case['case_id']] = dict(realized_fraction=len(selected)/len(marker_raw[0]),
                    selected_markers=len(selected), marker_library_mean=float((marker_raw+stress).sum(1).mean()),
                    marker_mass_fraction=float(np.mean((marker_raw+stress).sum(1)/(total+data['additions'].sum(1)))))
            entries.append(dict(key=key, n=unit['n'], raw=key+'_raw.bin', A=key+'_A.bin',
                                total=key+'_total.bin', marker_anchor=anchor, conditions=conditions))
            payload[rep, 'yields'] = marker_yield, base['ref']['pooled_yield']
        C.write_json(work/'jobs.json', dict(reference_folder=job['reference_folder'], entries=entries))
        r = options['r_info']
        D.run_r([r['executable'], '--vanilla', str(HERE/'grid_bridge.R'), str(work), r['library']],
                r['library'], folder/'bridge.log', options['r_timeout'])
        outputs = pd.read_csv(work/'outputs.csv', keep_default_na=False)
        C.require(not outputs.duplicated(['key','method']).any(), 'Duplicate bridge outputs')
        predictions = {}
        for row in outputs.itertuples():
            p = np.loadtxt(work/row.file, delimiter=',', ndmin=2)
            C.require(p.shape == (unit['n'],7), 'Prediction shape mismatch')
            predictions[row.key, row.method] = (p, row.reason)
        for rep in reps:
            null_predictions = {}
            marker_yield, full_yield = payload[rep, 'yields']
            for case in cases_for(job):
                data = payload[rep, case['case_id']]
                fitted = {}
                for method in METHODS:
                    reason = ''
                    if method == 'oracle':
                        p = data['p']
                    else:
                        q, reason = predictions[f'r{rep}_' + case['case_id'], method]
                        p = C.E.convert_available_rows(q, marker_yield if method=='marker_matched' else full_yield)
                    row, arr, fit = evaluate(data, p, 'full_six', reason)
                    fitted[method] = fit
                    attrs = dict(rep=rep, **case, method=method)
                    records.append(dict(attrs, **meta[rep, case['case_id']], **row))
                    if method.startswith('marker_'):
                        if case['outcome']=='null':
                            null_predictions[case['grid_id'],method] = p
                        elif case['outcome']=='up05':
                            prev = null_predictions[case['grid_id'],method]
                            finite = bool(np.isfinite(prev).all() and np.isfinite(p).all())
                            change = float(abs(p-prev).max()) if finite else np.nan
                            C.require(not finite or change < 1e-8, 'MRN paired invariance failed')
                            invariants.append(dict(rep=rep, grid_id=case['grid_id'], method=method,
                                                   finite=finite, max_change=change))
                    if case['grid_id']=='zero' and method in ('tmm','marker_matched','marker_denominator_only'):
                        # Preserve per-gene projections needed for a paired variance audit.
                        if fit['basic']:
                            audit = actual_span_hc3(data['y'], data['group'], p)
                            C.require(np.allclose(audit['beta'],fit['beta'],rtol=1e-9,atol=1e-9), 'Reusable API differs from frozen fit')
                            C.require(np.allclose(audit['se'],fit['se'],rtol=1e-9,atol=1e-9), 'Reusable HC3 differs from frozen fit')
                            residual_var = np.sum(audit['residual']**2,0) / audit['df']
                            hc3_max = audit['hc3_contributions'].max(0)
                        else:
                            residual_var = hc3_max = np.full(len(data['truth']),np.nan)
                        f = pd.DataFrame(arr)
                        f['gene_id'] = job['source']['genes'][job['test']]
                        f['residual_variance'] = residual_var
                        f['max_sample_hc3_contribution'] = hc3_max
                        for k,v in attrs.items(): f[k] = v
                        genes.append(f)
                for comparator in ('tmm','marker_denominator_only'):
                    check = coefficient_bound(data['y'],data['group'],data['truth'],fitted[comparator],fitted['marker_matched'])
                    theory.append(dict(rep=rep, case_id=case['case_id'], kind='coefficient', comparator=comparator, **check))
        outputs.to_csv(folder/'bridge_diagnostics.csv.gz', index=False, compression='gzip')
    pd.DataFrame(records).to_csv(folder/'replicates.csv.gz',index=False,compression='gzip')
    pd.concat(genes,ignore_index=True).to_csv(folder/'power_gene_fits.csv.gz',index=False,compression='gzip')
    C.write_csv(folder/'theory_checks.csv',theory)
    C.write_csv(folder/'invariance.csv',invariants)
    C.write_json(folder/'batch.json',dict(unit=unit,start=start,reps=reps,elapsed_seconds=time.monotonic()-clock))
    names=['replicates.csv.gz','power_gene_fits.csv.gz','theory_checks.csv','invariance.csv','batch.json','bridge_diagnostics.csv.gz','bridge.log']
    D.seal(folder,names,options['signature'])
    print('DONE',unit['unit_id'],start,'seconds',round(time.monotonic()-clock,1),flush=True)
    return dict(resumed=False,unit=unit['unit_id'],start=start)


def selection_manifest(jobs):
    rows=[]
    for i,job in enumerate(jobs):
        for point in stress_grid():
            ix=selected_rows(job,point['fraction'],point['selection'])
            rows.append(dict(job_index=i,reference=job['source']['label'],**point,
                             marker_count=len(job['panels']['balanced']),selected_count=len(ix),
                             selected_gene_ids=job['source']['genes'][job['panels']['balanced'][ix]].tolist()))
    return rows
