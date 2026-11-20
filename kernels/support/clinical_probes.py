"""A fixed-response multistep sensitivity check on the same 288 clinical donors."""
import json
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
from scipy import stats
from next_core import BASE, require, read_json, save_json, sha, safe_file
import bulk_data
import revision_core as C
import revision_driver as D
from extension_stats import fit_hc3
from experiment import METHODS, INVARIANT, read_predictions, evaluate_probe
from support import checkpoint, seal
from next_report import read, bh


def association(y, x, cov, p):
    fit = fit_hc3(y, x, p, covariates=cov)
    available = fit['basic'] & fit['good']
    pv = np.ones(y.shape[1])
    pv[available] = 2*stats.t.sf(abs(fit['beta'][available]/fit['se'][available]), fit['df'])
    return fit, available, pv, bh(pv)


def verify_baseline(genes, fit, valid, pv, qv, old):
    old = old.copy()
    old['gene_id'] = old.gene_id.astype(str)
    old = old.set_index('gene_id').reindex(list(map(str, genes)))
    require(len(old) == len(genes) and old.available.notna().all(), 'Clinical baseline gene identities changed')
    require(np.array_equal(valid, old.available.to_numpy(bool)), 'Clinical baseline availability changed')
    differences = {}
    for name, value in [('beta', np.where(valid, fit['beta'], np.nan)),
                        ('se', np.where(valid, fit['se'], np.nan)), ('p', pv), ('q', qv)]:
        actual = old[name].to_numpy(float)
        require(np.allclose(value, actual, atol=1e-8, rtol=1e-7, equal_nan=True),
                'Baseline '+name+' no longer matches the sealed parent. No probe result may replace it.')
        finite = np.isfinite(actual) & np.isfinite(value)
        differences[name+'_maximum_absolute_difference'] = float(np.max(abs(value[finite]-actual[finite]))) if finite.any() else None
    return differences


def run(args, parent, cfg, r):
    root = args.output/'clinical_multistep'
    if checkpoint(root, cfg['signature']):
        return
    root.mkdir(exist_ok=True)
    oldroot = parent/'clinical'
    original = read_json(oldroot/'protocol.json')
    print('CLINICAL: verify cached public counts and fixed sample/gene axes', flush=True)
    raw, metadata, _, provenance = bulk_data.prepare(args, root/'inputs')
    require(provenance == original['source'], 'Public clinical source identity differs from parent')
    meta, x, cov, _ = bulk_data.cohort(metadata)
    require(meta.sample_id.tolist() == original['sample_ids'] and len(meta) == 288, 'Clinical donor order changed')
    raw = raw.loc[:, meta.sample_id]
    genes = list(map(str, original['test_genes']))
    test_rows = raw.index.get_indexer(genes)
    require((test_rows >= 0).all() and len(genes) == 16161, 'Clinical test genes missing')
    old_assoc = read(oldroot/'disease_associations.csv.gz')
    yields = read(oldroot/'fraction_yields.csv')
    original_rna = read(oldroot/'rna_predictions.csv')
    original_props = read(oldroot/'composition_predictions.csv')
    hs = cfg['settings']['clinical_derivative_steps']
    probes = cfg['settings']['clinical_probe_steps']
    steps = sorted({0., *probes, *(s for h in hs for s in (-h, h, -h/2, h/2))})
    diagnostics, approximations, finite_records, baselines, inputs, invariance = [], [], [], [], [], []
    y_common = None
    for reference in ['LIBD', 'SEAAD']:
        previous = oldroot/reference
        source_job = read_json(previous/'jobs.json')
        require(len(source_job['entries']) == 1, 'Unexpected clinical job count')
        axes, oldentry = source_job['axes'], source_job['entries'][0]
        require(axes['normalization_genes'] == len(raw) and np.array_equal(np.array(axes['test_rows'])-1, test_rows),
                'Clinical count/test axis changed')
        require(oldentry['n'] == len(meta) and np.array_equal(oldentry['group'], x), 'Clinical exposure changed')
        nmarkers = len(axes['marker_rows'])
        marker_rows = np.array(axes['marker_rows'])-1
        for field, shape in [('raw', (len(raw), len(meta))), ('A', (nmarkers, 7)), ('total', (len(meta),))]:
            path = safe_file(previous, oldentry[field])
            require(path.stat().st_size == int(np.prod(shape))*8, 'Unexpected clinical binary size')
            value = np.fromfile(path, dtype='<f8').reshape(shape, order='F')
            require(np.isfinite(value).all() and (value >= 0).all(), 'Invalid clinical input matrix')
            if field == 'raw':
                require(np.array_equal(value, raw.to_numpy(float)), 'Clinical binary differs from deposited counts')
            elif field == 'total':
                require(np.array_equal(value, raw.sum(0).to_numpy(float)), 'Clinical totals differ from deposited counts')
            inputs.append(dict(reference=reference, input=field, bytes=path.stat().st_size, sha256=sha(path),
                               identity_check='public-count equality' if field != 'A' else 'fixed parent reference; baseline predictions checked'))
        work = root/reference
        work.mkdir(exist_ok=True)
        entry = dict(oldentry, steps=steps)
        for field in ('raw', 'A', 'total'):
            # Copy persistent parent inputs. The parent is never edited.
            shutil.copy2(safe_file(previous, oldentry[field]), work/Path(oldentry[field]).name)
            entry[field] = Path(oldentry[field]).name
        C.write_json(work/'jobs.json', dict(axes=axes, entries=[entry]))
        print('CLINICAL', reference, len(steps), 'steps; outcomes and covariates held fixed', flush=True)
        D.run_r([r['executable'], '--vanilla', str(BASE/'bridge.R'), str(work), r['library']],
                r['library'], work/'bridge.log', args.r_timeout)
        key = entry['key']
        pp, rr = read_predictions(work, {(key, m, t) for m in METHODS for t in steps})
        denominator = pp[key, 'tmm', 0.][1]
        require(np.isfinite(denominator).all() and (denominator > 0).all(), 'Invalid baseline denominator')
        response = np.log2(raw.iloc[test_rows].to_numpy(float).T/denominator[:, None]*1e6+1)
        if y_common is None:
            y_common = response
        else:
            require(np.array_equal(y_common, response), 'Clinical outcome changed between references')
        marker = raw.iloc[marker_rows].to_numpy(float).T
        for method in METHODS:
            q, d, reason, scale = pp[key, method, 0.]
            yg = yields[yields.reference.eq(reference) & yields.method.eq(method)].set_index('cell_class')
            yield_values = yg.reindex(C.M.CT).yield_value.to_numpy(float)
            require(np.isfinite(yield_values).all() and (yield_values > 0).all(), 'Missing fixed clinical yields')
            p = C.E.convert_available_rows(q, yield_values)
            for frame, value in [(original_rna, q), (original_props, p)]:
                old = frame[frame.reference.eq(reference) & frame.method.eq(method)].set_index('sample_id')
                require(old.index.is_unique and set(old.index) == set(meta.sample_id), 'Prediction sample identities changed')
                old = old.reindex(meta.sample_id)[C.M.CT].to_numpy(float)
                require(np.allclose(value, old, atol=1e-9, rtol=1e-7, equal_nan=True),
                        'Baseline clinical fractions no longer match sealed parent')
            fit, valid, pv, qv = association(y_common, x, cov, p)
            old = old_assoc[old_assoc.reference.eq(reference) & old_assoc.method.eq(method)]
            differences = verify_baseline(genes, fit, valid, pv, qv, old)
            baselines.append(dict(reference=reference, method=method, **differences,
                                  available_genes=int(valid.sum()), discoveries=int((valid & (qv <= .05)).sum())))
            props, denoms = {}, {0.: d}
            for step in steps:
                if step == 0.:
                    continue
                qt, dt, _, _ = pp[key, method, step]
                props[step] = C.E.convert_available_rows(qt, yield_values)
                denoms[step] = dt
                if method in INVARIANT:
                    require(np.allclose(q, qt, atol=1e-10, rtol=0, equal_nan=True), 'Complete-exclusion invariance failed')
                    invariance.append(dict(reference=reference, method=method, step=step,
                        max_prediction_difference=float(np.nanmax(abs(q-qt))) if np.isfinite(q).any() else np.nan))
            for h in hs:
                chosen = {s: props[s] for s in {-h, h, -h/2, h/2, *probes}}
                denom_chosen = {s: denoms[s] for s in [0., *chosen]}
                score, changes = evaluate_probe(dict(group=x, y=y_common), p, chosen, denom_chosen, marker, scale, h, fit, cov)
                diagnostics.append(dict(reference=reference, method=method, derivative_step=h, **score))
                approximations.extend(dict(reference=reference, method=method, derivative_step=h, **v) for v in changes)
            for step in probes:
                fp, vp, _, qp = association(y_common, x, cov, props[step])
                both = valid & vp
                b0, bt = valid & (qv <= .05), vp & (qp <= .05)
                change = fp['beta'][both]-fit['beta'][both]
                finite_records.append(dict(reference=reference, method=method, log_count_step=step,
                    baseline_available=int(valid.sum()), perturbed_available=int(vp.sum()), paired_available=int(both.sum()),
                    coefficient_change_rms=float(np.sqrt(np.mean(change**2))) if both.any() else np.nan,
                    median_absolute_standardized_change=float(np.median(abs(change)/fit['se'][both])) if both.any() else np.nan,
                    baseline_discoveries=int(b0.sum()), probe_discoveries=int(bt.sum()),
                    threshold_crossings=int((b0 ^ bt).sum()), threshold_crossing_fraction=float((b0 ^ bt).mean()),
                    interpretation='Finite perturbation of tested counts in covariate construction only; fixed outcome and BH family'))
        rr.to_csv(work/'prediction_manifest.csv', index=False)
    for name, records in [('derivative_steps', diagnostics), ('linear_approximation', approximations),
                          ('finite_perturbations', finite_records), ('baseline_identity', baselines),
                          ('input_identity', inputs), ('invariance', invariance)]:
        pd.DataFrame(records).to_csv(root/(name+'.csv'), index=False)
    save_json(root/'scope.json', dict(samples=288, tested_genes=16161, derivative_steps=hs,
        positive_probe_steps=probes, number_of_r_steps=len(steps),
        purpose='Step-size reliability and fixed-response finite sensitivity; not clinical FDR validation',
        selection='All fixed step sizes reported; no best-step selection',
        public_counts_reverified=True, parent_fractions_and_associations_reverified=True))
    names = [p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and
             p.name != 'complete.json' and p.suffix != '.bin']
    seal(root, names, cfg['signature'], samples=288, genes=16161)
