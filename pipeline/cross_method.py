"""Paired, fixed-response sensitivity experiment with an official second estimator."""
import concurrent.futures
import importlib
import json
import multiprocessing
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np
import pandas as pd

from .io import ROOT, atomic_json, digest, lock, read_json, receipt, separate, sha, verify_receipt

MUSIC_COMMIT = 'f21fe67f5670d5e9fca0ad7550abaae3423eb59c'
METHODS = ('dtangle_tmm', 'dtangle_marker', 'music_marker', 'linear_svr_marker')
PROBES = ('baseline', 'test_counts', 'marker_counts')
GLOBAL = {}


def paired_metrics(before, after, y, p0, p1):
    """Change in an estimable OLS weight vector, not a clinical FDR score."""
    w0, w1 = before['weights'], after['weights']
    out = {'max_fraction_change': float(np.max(np.abs(p1 - p0)))}
    if not (np.isfinite(w0).all() and np.isfinite(w1).all()):
        out.update(weight_l2=np.nan, beta_rms=np.nan, beta_max=np.nan,
                   weight_identity_error=np.nan, se_ratio_median=np.nan)
        return out
    delta = (w1 - w0) @ y
    out.update(weight_l2=float(np.linalg.norm(w1 - w0)),
               beta_rms=float(np.sqrt(np.mean(delta**2))), beta_max=float(np.max(np.abs(delta))),
               weight_identity_error=float(np.max(np.abs(delta - (after['beta'] - before['beta'])))),
               se_ratio_median=float(np.median(after['se'] / before['se'])))
    if out['weight_identity_error'] > 1e-8 * (1 + out['beta_max']):
        raise ValueError('Fixed-response coefficient identity failed')
    return out


def _setup_runtime(work):
    from .runtime import materialize
    runtime = materialize(work / '.runtime')
    sys.path.insert(0, str(runtime / 'mechanism'))
    c = importlib.import_module('revision_core')
    d = importlib.import_module('extension_driver')
    s = importlib.import_module('extension_stats')
    return c, d, s


def _reference(job, work, verified):
    from scipy.io import mmwrite
    source = job['source']
    folder = work / 'references' / source['label']
    folder.mkdir(parents=True, exist_ok=True)
    marker = folder / 'source_receipt.json'
    identity = digest(dict(reference=source['label'], donors=job['refs'],
                           inputs={k:v for k,v in verified.items() if k.startswith(source['label'] + '/')}))
    if marker.is_file():
        record = read_json(marker)
        if record.get('identity') != identity:
            raise ValueError('Reference inputs or donors changed')
        verify_receipt(folder, record['files'])
    else:
        cells = np.concatenate([np.concatenate(source['masks'][d]) for d in job['refs']])
        x = source['x'][:, cells].tocsc()
        keep = np.asarray(x.sum(0)).ravel() > 0
        cells, x = cells[keep], x[:, keep]
        with tempfile.TemporaryDirectory(prefix='.source_', dir=folder) as temp:
            temp = Path(temp)
            mmwrite(temp / 'reference.mtx', x, field='integer')
            source['meta'].iloc[cells][['donor_id', 'cell_type']].to_csv(temp / 'cells.csv', index=False)
            pd.DataFrame({'gene_id': source['genes']}).to_csv(temp / 'genes.csv', index=False)
            for name in ('reference.mtx', 'cells.csv', 'genes.csv'):
                (temp / name).replace(folder / name)
        atomic_json(marker, {'identity': identity, 'files': receipt(folder, ['reference.mtx', 'cells.csv', 'genes.csv'])})
    return folder


def _valid_batch(folder, signature):
    path = folder / 'receipt.json'
    if not path.is_file():
        return False
    record = read_json(path)
    if record.get('protocol') != signature:
        raise ValueError('Batch belongs to another protocol: ' + str(folder))
    verify_receipt(folder, record['files'])
    return True


def _batch(task):
    unit, reps = task
    c, d, s, jobs, options, work, args, signature = (GLOBAL[k] for k in
        ('c', 'd', 's', 'jobs', 'options', 'work', 'args', 'signature'))
    job = jobs[unit['job_index']]
    folder = work / 'batches' / (unit['unit_id'] + '_' + str(reps[0]).zfill(4))
    if _valid_batch(folder, signature):
        return len(reps), True
    folder.mkdir(parents=True, exist_ok=True)
    marker, test, norm = job['panels']['balanced'], job['test'], job['norm_genes']
    norm_index = pd.Index(norm)
    rows, contrasts = [], []
    with tempfile.TemporaryDirectory(prefix='compute_', dir=folder) as temporary:
        temporary = Path(temporary)
        entries, payload = [], {}
        for rep in reps:
            base = c.draw_base(job, unit, rep, options)
            # Condition on the four observed reference donors, equally weighted,
            # for all methods and probes; no fictitious extra reference donors.
            base['ref'] = job['base']
            base['weights'] = np.ones(4, dtype=int)
            base['allocation'] = '1;1;1;1'
            a, marker_yield, anchor = d.marker_properties(job, base)
            stress, selected = d.stress_additions(job, base)
            probe_truth = np.zeros(len(test))
            probe_truth[job['injections']['frozen']] = .5
            expectation = options['nuclei'] * base['scale'][None, :] * base['group'][:, None] * probe_truth
            extra = c.M.rng_for(2026111701, 'fixed_response_test_probe', base['key']).poisson(expectation)
            raw_name = f'{rep}_raw.bin'
            c.write_bin(temporary / raw_name, base['draw']['counts'][:, norm].T)
            for effect in (0., .5):
                condition = dict(name='null' if effect == 0 else 'up05', effect=effect, pattern='up', panel='frozen')
                data = c.intervene(job, base, condition, options)
                key = f'r{rep}_{condition["name"]}'
                for tag, value in [('A', a), ('total', base['draw']['mass'].sum(1)),
                                   ('added', data['additions'].T), ('extra', extra.T), ('stress', stress.T)]:
                    c.write_bin(temporary / f'{key}_{tag}.bin', value)
                entries.append(dict(key=key, n=unit['n'], raw=raw_name, anchor=anchor))
                payload[key] = data, marker_yield, base['ref']['pooled_yield'], rep, effect
        plan = dict(libraries=[str(x.resolve()) for x in args.r_library], entries=entries,
                    reference=str(job['music_reference']), classes=c.M.CT,
                    norm_genes=job['source']['genes'][norm].tolist(),
                    marker_rows=(norm_index.get_indexer(marker) + 1).tolist(),
                    test_rows=(norm_index.get_indexer(test) + 1).tolist(),
                    marker_owner=c.E.official_markers(job).cell_type.to_numpy(int).tolist(),
                    music_commit=MUSIC_COMMIT,
                    validate_official=bool(reps[0] == 0 and unit['target_donor'] == job['targets'][0]))
        atomic_json(temporary / 'plan.json', plan)
        with (folder / 'engine.log').open('w', encoding='utf-8') as log:
            subprocess.run([args.rscript, str(ROOT / 'pipeline/cross_method.R'), str(temporary),
                            *[str(p.resolve()) for p in args.r_library]],
                           check=True, stdout=log, stderr=subprocess.STDOUT)
        audit = read_json(temporary / 'engine.json')
        for key, (data, marker_yield, full_yield, rep, effect) in payload.items():
            fixed_y = data['y'].copy()
            for method in METHODS:
                results = {}
                for probe in PROBES:
                    q = np.loadtxt(temporary / f'{key}_{probe}_{method}.csv', delimiter=',', ndmin=2)
                    props = q if method == 'music_marker' else c.E.convert_available_rows(
                        q, marker_yield if method in ('dtangle_marker', 'linear_svr_marker') else full_yield)
                    row, arrays, fit = s.evaluate(data, props, rule='full_six')
                    if not np.array_equal(data['y'], fixed_y):
                        raise ValueError('Outcome changed within a paired comparison')
                    meta = {k: unit[k] for k in ('reference_source', 'target_source', 'target_donor')}
                    meta.update(rep=rep, effect=effect, method=method)
                    rows.append(dict(meta, probe=probe, **row))
                    results[probe] = (props, arrays, fit)
                    if probe != 'baseline':
                        p0, a0, f0 = results['baseline']
                        comparison = paired_metrics(f0, fit, data['y'], p0, props)
                        comparison['bh_crossings'] = int(np.sum((a0['q'] <= .05) != (arrays['q'] <= .05)))
                        comparison['family_size'] = len(job['test'])
                        comparison['max_q_change'] = float(np.max(np.abs(a0['q'] - arrays['q'])))
                        contrasts.append(dict(meta, probe=probe, **comparison))
        pd.DataFrame(rows).to_csv(folder / 'metrics.csv', index=False)
        pd.DataFrame(contrasts).to_csv(folder / 'contrasts.csv', index=False)
        atomic_json(folder / 'engine.json', audit)
    atomic_json(folder / 'receipt.json', dict(protocol=signature, datasets=len(reps) * 2,
        files=receipt(folder, ['metrics.csv', 'contrasts.csv', 'engine.json'])))
    return len(reps), False


def summarize(work, signature):
    from scipy.stats import t
    frames = {}
    for name in ('metrics', 'contrasts'):
        pieces = []
        for folder in sorted((work / 'batches').iterdir()):
            if _valid_batch(folder, signature):
                pieces.append(pd.read_csv(folder / (name + '.csv')))
        frames[name] = pd.concat(pieces, ignore_index=True)
    keys = ['reference_source', 'target_source', 'effect', 'method', 'probe']
    identity = keys + ['target_donor', 'rep']
    final = work / 'final'
    final.mkdir(exist_ok=True)
    output = []
    for name, frame in frames.items():
        if frame.duplicated(identity).any():
            raise ValueError('Duplicate method-dataset record')
        frame.to_csv(final / (name + '.csv.gz'), index=False, compression='gzip')
        metrics = (['fdp', 'power', 'coverage', 'availability', 'proportion_rmse', 'mean_se'] if name == 'metrics'
                   else ['max_fraction_change', 'weight_l2', 'beta_rms', 'beta_max', 'se_ratio_median',
                         'bh_crossings', 'max_q_change', 'weight_identity_error'])
        for values, group in frame.groupby(keys, sort=True):
            donors = group.groupby('target_donor', sort=True)
            for metric in metrics:
                dmeans = donors[metric].mean()
                variances = donors[metric].var(ddof=1)
                counts = donors[metric].count()
                valid = counts.gt(0)
                n = int(valid.sum())
                mean = float(dmeans[valid].mean()) if n else np.nan
                se = float(np.sqrt((variances[valid] / counts[valid]).sum()) / n) if n and counts[valid].min() > 1 else np.nan
                critical = float(t.ppf(.975, counts[valid].min() - 1)) if n and counts[valid].min() > 1 else np.nan
                output.append(dict(zip(keys, values), metric=metric, donor_equal_mean=mean,
                    conditional_mc_se=se, conditional_mc_low=mean-critical*se, conditional_mc_high=mean+critical*se,
                    min_donor_mean=float(dmeans[valid].min()) if n else np.nan,
                    max_donor_mean=float(dmeans[valid].max()) if n else np.nan,
                    min_value=float(group[metric].min()), max_value=float(group[metric].max()),
                    donors=n, datasets=len(group)))
    pd.DataFrame(output).to_csv(final / 'summary.csv', index=False)
    paired = paired_precision(frames['metrics'])
    paired.to_csv(final / 'paired_precision.csv', index=False)
    atomic_json(final / 'receipt.json', dict(protocol=signature,
        uncertainty='Conditional Monte Carlo variation with observed donors and references fixed; donor ranges reported separately.',
        files=receipt(final, ['metrics.csv.gz', 'contrasts.csv.gz', 'summary.csv', 'paired_precision.csv'])))
    return final


def paired_precision(frame):
    """Monte Carlo uncertainty of within-dataset contrasts with donors fixed."""
    from scipy.stats import t
    identity = ['reference_source', 'target_source', 'target_donor', 'rep', 'effect']
    records = []
    for method in METHODS:
        current = frame[frame.method.eq(method)]
        base = current[current.probe.eq('baseline')]
        for probe in PROBES[1:]:
            pair = base.merge(current[current.probe.eq(probe)], on=identity,
                              validate='one_to_one', suffixes=('_base', '_probe'))
            if len(pair) != len(base):
                raise ValueError('Incomplete paired precision contrast')
            for keys, group in pair.groupby(['reference_source', 'target_source', 'effect']):
                for metric in ('fdp', 'power', 'coverage', 'proportion_rmse'):
                    values = group[metric + '_probe'] - group[metric + '_base']
                    groups = values.groupby(group.target_donor)
                    counts, means, variances = groups.count(), groups.mean(), groups.var(ddof=1)
                    good = counts.gt(0)
                    n = int(good.sum())
                    mean = float(means[good].mean()) if n else np.nan
                    se = (float(np.sqrt((variances[good]/counts[good]).sum())/n)
                          if n and counts[good].min() > 1 else np.nan)
                    half = float(t.ppf(.975, counts[good].min()-1)*se) if n else np.nan
                    records.append(dict(reference_source=keys[0], target_source=keys[1], effect=keys[2],
                        method=method, probe=probe, metric=metric, donor_equal_change=mean,
                        conditional_mc_se=se, conditional_mc_low=mean-half, conditional_mc_high=mean+half,
                        half_width=half, precision_target=.02 if metric=='fdp' else np.nan,
                        target_met=bool(half<=.02) if metric=='fdp' and np.isfinite(half) else None,
                        donors=n, paired_datasets=len(group)))
    return pd.DataFrame(records)


def execute(args):
    if min(args.workers, args.repetitions, args.batch_size) < 1:
        raise ValueError('Workers, repetitions and batch size must be positive')
    if os.name != 'posix':
        raise ValueError('The raw-count benchmark uses Linux fork workers; analytical functions are portable')
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[name] = '1'
    work = separate(args.work, args.inputs)
    work.mkdir(parents=True, exist_ok=True)
    with lock(work / '.writer.lock'):
        from .inputs import verify_inputs
        from .run import python_environment
        environment = dict(python=sys.version, packages=python_environment(),
            R=json.loads(subprocess.check_output([args.rscript, str(ROOT / 'pipeline/cross_method.R'),
                '--environment', *[str(p.resolve()) for p in args.r_library]], text=True)))
        verified = verify_inputs(args.inputs, include_bulk=False)
        c, d, s = _setup_runtime(work)
        folders = {label: args.inputs / label for label in ('LIBD', 'SEAAD')}
        sources, jobs = c.full_inputs(folders, args.inputs / 'gene_identity_map.csv')
        external = c.S.read_external(args.inputs / 'external/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5',
            args.inputs / 'external/GSE174367_snRNA-seq_cell_meta.csv.gz', next(iter(sources.values()))['genes'],
            work / 'external_inputs')
        jobs += c.S.make_external_jobs(jobs, external, work)
        jobs = [c.prepare_job(j, work) for j in jobs]
        for j in jobs:
            j['music_reference'] = _reference(j, work, verified)
        for reference in sorted({j['music_reference'] for j in jobs}):
            subprocess.run([args.rscript, str(ROOT / 'pipeline/cross_method.R'), '--prepare', str(reference),
                            *[str(p.resolve()) for p in args.r_library]], check=True)
        units = c.plan_units(jobs, smoke=args.smoke)
        for unit in units:
            unit.pop('conditions', None)
        reps = min(args.repetitions, 2) if args.smoke else args.repetitions
        options = dict(c.SEEDS)
        protocol = dict(schema=2, repetitions=reps, batch_size=args.batch_size, n=80, effects=[0., .5], probes=list(PROBES), methods=list(METHODS),
            reference='Four original reference donors, equally weighted and fixed across replicates and methods.',
            test_probe='Independent Poisson additions, effect 0.5 on the fixed 20 test genes in exposed samples; outcome held fixed.',
            marker_probe='Original 10% per-class marker selection, 50% expected count increase in exposed samples; outcome held fixed.',
            test_family='Original disjoint 100-gene family; unavailable tests have p=1.',
            precision='Fixed repetitions before fitting; paired conditional Monte Carlo intervals. FDR-change half-width target 0.02 is reported, never an outcome-dependent stopping rule.',
            music='Official music_prop body with memoized music_basis, fixed disjoint marker panel, native nuclear proportions; checked against the uncached API.',
            music_commit=MUSIC_COMMIT, seeds=dict(options, probe_seed=2026111701), inputs=verified,
            environment=environment,
            units=units, software_smoke=bool(args.smoke),
            reference_hashes={str(p.relative_to(work)): sha(p) for p in sorted((work / 'references').rglob('*')) if p.is_file()},
            linear_svr='e1071 1.7.17 LIBSVM linear nu-SVR; nu selected from .25,.5,.75 by reconstruction error, no disease labels; nonnegative closure then marker-yield conversion.',
            code={p: sha(ROOT / p) for p in ('pipeline/cross_method.py', 'pipeline/cross_method.R', 'pipeline/estimators.R', 'kernels/layout.json')})
        signature = digest(protocol)
        path = work / 'protocol.json'
        if path.exists():
            if not args.resume or read_json(path) != protocol:
                raise ValueError('Existing work requires --resume with an identical protocol')
        else:
            atomic_json(path, protocol)
        GLOBAL.update(c=c, d=d, s=s, jobs=jobs, options=options, work=work, args=args, signature=signature)
        tasks = [(unit, list(range(start, min(reps, start + args.batch_size))))
                 for unit in units for start in range(0, reps, args.batch_size)]
        start = time.monotonic()
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context('fork')) as pool:
            futures = [pool.submit(_batch, task) for task in tasks]
            for done, future in enumerate(concurrent.futures.as_completed(futures), 1):
                count, reused = future.result()
                print('CROSS_METHOD', f'{done}/{len(tasks)}', 'reused' if reused else 'computed',
                      count * 2, 'datasets', round(time.monotonic()-start, 1), 's', flush=True)
        final = summarize(work, signature)
        atomic_json(work / 'complete.json', dict(protocol=signature, batches=len(tasks), datasets=len(units)*reps*2,
                    final=str(final)))
        print('COMPLETE', final, flush=True)
