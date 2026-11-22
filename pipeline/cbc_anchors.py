"""Paired anchor benchmark with empirical doses, disjoint donors and compact outputs."""
import multiprocessing as mp
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.special import softmax
from scipy.stats import rankdata, t

from . import cbc_inputs as inp
from .download import download
from .io import ROOT, atomic_json, digest, lock, read_json, receipt, separate, sha, verify_receipt

METHODS = ('tmm', 'complete', 'mrn', 'anchor_equal', 'anchor_weighted', 'oracle')
LABELS = ('TMM', 'Complete exclusion', 'MRN', 'Equal anchors', 'Weighted anchors', 'Oracle')
DIRECTIONS = ('Astro', 'Micro')
PROBES = ('baseline', 'measured_test', 'marker')
N = 80
REPS = 60
BATCH = 10
SEED = 2026112107
GAMMA = 0.94339016910671658
LIMITS = dict(fdr_margin=.02, coverage_margin=.02, fdr_upper=.07, availability=.99)
STATE = {}


def tmm_denominator(counts, total):
    """Weighted TMM with edgeR's reference, trimming and geometric centering rules."""
    a, d = np.asarray(counts, float), np.asarray(total, float)
    if a.ndim != 2 or d.shape != (len(a),) or not np.isfinite(a).all() or (a < 0).any() or (d <= 0).any():
        raise ValueError('TMM needs nonnegative counts and strictly positive library sizes')
    a = a[:, np.any(a > 0, axis=0)]
    if a.shape[1] == 0 or len(d) == 1:
        return d.copy()
    x = a / d[:, None]
    upper = np.quantile(x, .75, axis=1)
    j = np.argmax(np.sqrt(a).sum(1)) if np.median(upper) < 1e-20 else np.argmin(abs(upper - upper.mean()))
    with np.errstate(divide='ignore', invalid='ignore'):
        m = np.log2(x / x[j])
        avg = (np.log2(x) + np.log2(x[j])) / 2
        var = (1-x)/a + (1-x[j])/a[j]
    valid = np.isfinite(m) & np.isfinite(avg) & (avg > -1e10)
    mc = rankdata(np.where(valid, m, np.nan), axis=1, nan_policy='omit')
    ac = rankdata(np.where(valid, avg, np.nan), axis=1, nan_policy='omit')
    n = valid.sum(1, keepdims=True)
    ml, al = np.floor(.3*n) + 1, np.floor(.05*n) + 1
    keep = valid & (mc >= ml) & (mc <= n + 1 - ml) & (ac >= al) & (ac <= n + 1 - al)
    w = np.divide(1., var, out=np.zeros_like(var), where=keep & (var > 0))
    log_f = np.divide((np.where(keep, m, 0) * w).sum(1), w.sum(1),
                      out=np.zeros(len(d)), where=w.sum(1) > 0)
    log_f[np.max(np.abs(np.where(valid, m, 0)), axis=1) < 1e-6] = 0
    return d * np.exp2(log_f - log_f.mean())


def anchor_denominators(counts, freq, weights):
    """Disjoint, fixed-reference weighted-median scale; zeros give no positive ratio."""
    counts = np.asarray(counts, float)
    weights = np.asarray(weights, float)
    freq = np.asarray(freq, float)
    if (counts.ndim != 2 or counts.shape[1] != len(freq) or weights.shape != freq.shape
            or not np.isfinite(counts).all() or (counts < 0).any() or (freq <= 0).any()
            or (weights <= 0).any() or not np.isfinite(weights).all()):
        raise ValueError('Invalid fixed anchors')
    support = (counts > 0).sum(1) >= (len(freq) + 1) // 2
    ratio = np.divide(counts, freq, where=counts > 0, out=np.full_like(counts, np.inf))
    order = ratio.argsort(axis=1, kind='stable')
    sorted_w = np.take_along_axis((counts > 0) * weights, order, axis=1)
    cumulative = sorted_w.cumsum(1)
    k = (cumulative >= cumulative[:, -1:] / 2).argmax(1)
    d = np.take_along_axis(ratio, order, axis=1)[np.arange(len(counts)), k]
    return np.where(support, d, np.nan)


def marker_prediction(markers, den, ref, owners, marker_scaled=False):
    den, reference = np.asarray(den, float), np.asarray(ref, float)
    if not np.isfinite(den).all() or (den <= 0).any():
        return np.full((len(markers), len(reference)), np.nan)
    scale = 1e6
    class_mass = reference.sum(1)
    if marker_scaled:
        scale *= class_mass.mean()
        reference = reference / class_mass[:, None]
    a = np.log2(1 + markers / den[:, None] * scale)
    b = np.log2(1 + reference * scale)
    score = np.column_stack([(a[:, owners == k] - b[k, owners == k]).mean(1)
                             for k in range(len(reference))]) * np.log(2) / GAMMA
    if marker_scaled:
        score -= np.log(class_mass)
    return softmax(score, axis=1)


def prediction(counts, total, model):
    mark, test, anchors = (model[k] for k in ('mark', 'test', 'anchor'))
    markers = counts[:, mark]
    anchor_counts = counts[:, anchors]
    include = np.ones(counts.shape[1], dtype=bool)
    include[test] = False
    dd = dict(tmm=tmm_denominator(counts, total),
              complete=tmm_denominator(counts[:, include], total - counts[:, test].sum(1)),
              mrn=markers.sum(1),
              anchor_equal=anchor_denominators(anchor_counts, model['abundance'], np.ones(len(anchors))),
              anchor_weighted=anchor_denominators(anchor_counts, model['abundance'], model['weights']))
    rf = model['ref'][:, mark]
    pred = {k: marker_prediction(markers, d, rf, model['owners'], k == 'mrn') for k, d in dd.items()}
    return pred, dd


def paired_draw(model, donor_index, rep, direction):
    profile = model['profile'][donor_index]
    rng = np.random.default_rng(np.random.SeedSequence([SEED, donor_index, rep, DIRECTIONS.index(direction)]))
    exposure = np.r_[np.zeros(N//2), np.ones(N//2)]
    mix = rng.dirichlet(np.array([.14,.05,.32,.15,.08,.06,.20]) * 45, size=N)
    mix[:, model['groups'].index(direction)] *= np.exp(exposure)
    mix /= mix.sum(1, keepdims=True)
    expected = mix @ profile
    pr = np.column_stack([expected, np.maximum(1-expected.sum(1), 0)])
    pr /= pr.sum(1, keepdims=True)
    depth = np.rint(rng.lognormal(np.log(50000), .38, N)).astype(int)
    raw = np.array([rng.multinomial(d, p) for d,p in zip(depth,pr)])[:, :-1]
    test = model['test']
    nuisance = .25 * mix @ np.log2(1+profile[:, test]*1e6)
    noise = rng.normal(0,.25,(N,len(test)))
    truth = np.zeros(len(test)); truth[:20] = .25
    base = dict(p=mix, group=exposure, nuisance=nuisance, noise=noise)
    scenarios = dict(null=dict(base, truth=np.zeros(len(test)), y=nuisance+noise),
                     sparse=dict(base, truth=truth, y=nuisance+noise+exposure[:, None]*truth))
    active = exposure == 1
    # The donor-weighted technical fold vector is fixed across all simulations.
    measured = raw.copy()
    measured[np.ix_(active,test)] = np.rint(raw[np.ix_(active,test)] * model['fold'][None,:]).astype(int)
    stress = raw.copy()
    subset = model['stress']
    stress[np.ix_(active,subset)] = np.rint(stress[np.ix_(active,subset)] * 1.5).astype(int)
    return raw, depth, scenarios, {
        'baseline':(raw, depth),
        'measured_test':(measured, depth + (measured-raw).sum(1)),
        'marker':(stress, depth + (stress-raw).sum(1)),
    }


def init_worker(runtime, model, out, identity):
    sys.path.insert(0, str(runtime / 'mechanism'))
    from extension_stats import evaluate
    STATE.update(eval=evaluate, model=model, out=out, identity=identity)


def batch_file(out, task):
    d, direction, reps = task
    return out / '.batches' / f'{d:02d}_{direction}_{reps[0]:04d}'


def audit_batch(folder, identity):
    path = folder / 'receipt.json'
    if not path.exists():
        return False
    done = read_json(path)
    if done.get('identity') != identity:
        raise ValueError('Batch belongs to a different protocol: ' + str(folder))
    verify_receipt(folder, done['files'])
    return True


def batch(task):
    donor, direction, reps = task
    model, out, identity = (STATE[k] for k in ('model','out','identity'))
    folder = batch_file(out, task)
    if audit_batch(folder, identity):
        return str(folder)
    folder.mkdir(parents=True, exist_ok=True)
    rows = []
    for rep in reps:
        baseline_raw, _, outcomes, probes = paired_draw(model, donor, rep, direction)
        fits = {}; d0 = {}
        for probe, (raw, total) in probes.items():
            pred, effective = prediction(raw,total,model)
            pred['oracle'] = outcomes['null']['p']
            if probe == 'baseline':
                d0 = effective
            for method in METHODS:
                if probe == 'measured_test' and method not in ('tmm','oracle'):
                    if not np.array_equal(effective[method], d0[method], equal_nan=True):
                        raise ValueError('Test counts entered a protected denominator')
                for outcome, data in outcomes.items():
                    row, _, fit = STATE['eval'](data, pred[method])
                    key = method, outcome
                    if probe == 'baseline':
                        fits[key] = fit
                    base = fits[key]
                    available = base['good'] & fit['good']
                    error = np.nan
                    beta_rms = np.nan
                    if available.any():
                        delta = (fit['weights'] - base['weights']) @ data['y'][:,available]
                        error = np.max(np.abs(delta - (fit['beta'] - base['beta'])[available]))
                        beta_rms = np.sqrt(np.mean(delta**2))
                        if error > 1e-9:
                            raise ValueError('Eq. (5) and the actual fits disagree')
                    change = effective[method] / d0[method] - 1 if method != 'oracle' else np.zeros(N)
                    active = data['group'].astype(bool)
                    def med(a):
                        return np.median(a[np.isfinite(a)]) if np.isfinite(a).any() else np.nan
                    rows.append(dict(donor=model['held'][donor], direction=direction, rep=rep,
                        probe=probe, outcome=outcome, method=method, fdr=row['fdp'], power=row['power'],
                        coverage=row['delivered_coverage'], availability=row['availability'],
                        fraction_rmse=row['proportion_rmse'], beta_rms=beta_rms, identity_error=error,
                        exposed_den_change=med(change[active]), control_den_change=med(change[~active]),
                        anchor_reads=np.median(baseline_raw[:,model['anchor']].sum(1)),
                        anchor_detected_min=np.min((baseline_raw[:,model['anchor']] > 0).sum(1))))
    pd.DataFrame(rows).to_csv(folder/'metrics.csv', index=False)
    atomic_json(folder/'receipt.json', dict(identity=identity, files=receipt(folder,['metrics.csv'])))
    return str(folder)


def interval(values):
    a = np.asarray(values, float)
    a = a[np.isfinite(a)]
    if not a.size:
        return np.nan, np.nan, np.nan
    m = float(a.mean())
    if a.size == 1:
        return m, np.nan, np.nan
    e = float(t.ppf(.975, len(a)-1)*a.std(ddof=1)/np.sqrt(len(a)))
    return m,m-e,m+e


def summarize(frame):
    keys = ['direction','probe','outcome','method']
    metrics = ['fdr','power','coverage','availability','fraction_rmse','beta_rms',
               'exposed_den_change','control_den_change','anchor_reads']
    donors = frame.groupby([*keys,'donor'], sort=False)[metrics].mean().reset_index()
    rows = []
    for group, part in donors.groupby(keys,sort=False):
        row = dict(zip(keys,group), donors=len(part))
        for metric in metrics:
            mean, low, high = interval(part[metric])
            row.update({metric:mean, metric+'_lo':low, metric+'_hi':high})
        rows.append(row)
    summary = pd.DataFrame(rows)
    contrasts = []
    for group, part in donors.groupby(keys[:-1],sort=False):
        for method, comparator in [('anchor_equal','mrn'),('anchor_weighted','mrn'),('anchor_weighted','anchor_equal')]:
            a = part[part.method == method].set_index('donor')
            b = part[part.method == comparator].set_index('donor')
            if set(a.index) != set(b.index):
                raise ValueError('Paired contrast lost donors')
            row = dict(zip(keys[:-1],group), method=method, comparator=comparator, donors=len(a))
            for metric in ('fdr','power','coverage'):
                mean, low, high = interval(a[metric] - b[metric])
                row.update({metric:mean,metric+'_lo':low,metric+'_hi':high})
            contrasts.append(row)
    return summary, pd.DataFrame(contrasts)


def decision(summary, contrasts, smoke):
    a = summary.loc[summary.method == 'anchor_weighted']
    c = contrasts.loc[(contrasts.method == 'anchor_weighted') & (contrasts.comparator == 'mrn')]
    p = c.loc[(c.direction == 'Micro') & (c.probe == 'marker') & (c.outcome == 'sparse')]
    gates = dict(fdr_absolute=bool(a.fdr_hi.le(LIMITS['fdr_upper']).all()),
                 fdr_vs_mrn=bool(c.fdr_hi.le(LIMITS['fdr_margin']).all()),
                 coverage_vs_mrn=bool(c.coverage_lo.ge(-LIMITS['coverage_margin']).all()),
                 availability=bool(a.availability.ge(LIMITS['availability']).all()),
                 primary_power=bool(len(p)==1 and p.power_lo.gt(0).all()))
    return dict(eligible=not smoke and all(gates.values()), criteria=gates, limits=LIMITS,
        primary=dict(direction='Micro',probe='marker',outcome='sparse',reference='mrn'))


def execute(args):
    from .inputs import specification
    from .runtime import materialize
    from .publication import export_published

    if args.workers < 1:
        raise ValueError('Workers must be positive')
    bulk = Path(args.bulk_cache)/'GSE125583_RAW.tar'
    public = inp.checked(args.imaging_cache, args.library_types, args.allow_downloads)
    if not bulk.is_file():
        spec = specification()['bulk/GSE125583_RAW.tar']
        if not args.allow_downloads:
            raise ValueError('Missing GSE125583_RAW.tar; supply --bulk-cache or enable pinned downloads')
        download(spec['url'],bulk,spec['sha256'])
    out = separate(args.output,args.pseudobulk,args.imaging_cache,args.bulk_cache)
    if out.exists() and not args.resume and any(out.iterdir()):
        raise ValueError('Use a new output or the same one with --resume')
    out.mkdir(parents=True,exist_ok=True)
    with lock(out/'.writer.lock'):
        test, markers = inp.panels()
        observed, shift, mapped = inp.paired_counts(public,test,markers)
        model = inp.training_inputs(args.pseudobulk,test,markers)
        model['fold'] = np.exp(np.clip(shift, -np.log(4), np.log(4)))
        model['stress'] = np.concatenate([model['mark'][model['owners']==k][:max(1,round(sum(model['owners']==k)/10))]
                                          for k in range(7)])
        runtime = materialize(out/'.runtime')
        sys.path.insert(0,str(runtime/'mechanism'))
        try:
            from revision_core import M
            if tuple(M.CT) != model['groups']:
                raise ValueError('Frozen marker ownership uses a different class order')
        finally:
            sys.path.pop(0)
        mapping = ROOT/'kernels/coupling/gene_id_mapping_snapshot.csv.gz'
        if not mapping.is_file():
            raise ValueError('Bundled clinical gene-identity snapshot is missing')
        n_reps = 2 if args.smoke else REPS
        d_count = 2 if args.smoke else len(model['held'])
        protocol = dict(seed=SEED, reps=n_reps, samples=N, donors=model['held'][:d_count],
            reference_donors=model['train'], directions=DIRECTIONS, probes=PROBES, test_genes=len(test),
            marker_genes=len(markers), anchor_size=len(model['anchor']), mapped_observed_test_genes=mapped,
            universe_genes=len(model['genes']), observed_blocks=observed.block.nunique(),
            observed_donors=observed.donor.nunique(), observed_pairs=len(observed),
            doses='Equal-depth RiboZeroGold versus poly(A) in matched LIBD blocks and processing rounds; log ratio of CPM + 0.5, donor medians, fourfold cap.',
            checksums=dict(pseudobulk=inp.PB_SHA, clinical=sha(bulk), mapping=sha(mapping),
                           **{key:sha(path) for key,path in public.items()}),
            recipe=receipt(ROOT,['pipeline/cbc_anchors.py','pipeline/cbc_inputs.py','kernels/mechanism/extension_stats.py',
                                 'publication/input_definitions/fixed_marker_panels.csv',
                                 'publication/input_definitions/complete_test_gene_families.csv']),
            libraries=dict(numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),
            formats=tuple(dict.fromkeys(args.formats)), article_replay=not args.skip_replay,
            formal_result=not args.smoke)
        identity = digest(protocol)
        prev = out/'protocol.json'
        if prev.exists() and digest(read_json(prev)) != identity:
            raise ValueError('Inputs or analysis code changed; use a new output directory')
        atomic_json(prev,protocol)
        finish = out/'final'/'receipt.json'
        if finish.exists():
            prior = read_json(finish)
            if prior['identity'] != identity:
                raise ValueError('Output belongs to a different protocol')
            verify_receipt(out/'final',prior['files'])
            print('Verified completed CBC upgrade',out/'final',flush=True)
            return
        plan = [(i,d,list(range(j,min(j+BATCH,n_reps)))) for i in range(d_count) for d in DIRECTIONS
                for j in range(0,n_reps,BATCH)]
        pending = [q for q in plan if not audit_batch(batch_file(out,q),identity)]
        print(f'CBC anchors: {len(plan)-len(pending)}/{len(plan)} batches verified; {len(pending)} pending',flush=True)
        for var in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
            os.environ[var] = '1'
        if pending:
            ctx = mp.get_context('spawn')
            with ctx.Pool(min(args.workers,len(pending)),init_worker,(runtime,model,out,identity)) as pool:
                for k,_ in enumerate(pool.imap_unordered(batch,pending),len(plan)-len(pending)+1):
                    if k % 4 == 0 or k == len(plan):
                        print(f'CBC anchor batches {k}/{len(plan)}',flush=True)
        frames=[]
        for task in plan:
            folder=batch_file(out,task)
            if not audit_batch(folder,identity):
                raise ValueError('Missing expected batch: ' + str(folder))
            frames.append(pd.read_csv(folder/'metrics.csv',keep_default_na=False,na_values=['']))
        records=pd.concat(frames,ignore_index=True)
        expected = len(DIRECTIONS)*d_count*n_reps*len(PROBES)*2*len(METHODS)
        if len(records)!=expected or records.duplicated(['direction','donor','rep','probe','outcome','method']).any():
            raise ValueError('Missing or duplicate held-out comparisons')
        s,c=summarize(records)
        if len(s)!=len(DIRECTIONS)*len(PROBES)*2*len(METHODS) or len(c)!=len(DIRECTIONS)*len(PROBES)*2*3:
            raise ValueError('Missing summarized null or sparse panels')
        capacity = inp.clinical_count_capacity(bulk,model['panel'].gene_id,mapping)
        f = out/'final'
        data = f/'data'
        data.mkdir(parents=True,exist_ok=True)
        s.to_csv(data/'Anchor_Performance.csv',index=False)
        c.to_csv(data/'Anchor_Paired_Contrasts.csv',index=False)
        capacity.to_csv(data/'Clinical_Anchor_Capacity.csv',index=False)
        observed.to_csv(data/'Observed_Bulk_Pairs.csv',index=False)
        model['panel'].to_csv(data/'Fixed_Anchor_Panel.csv',index=False)
        model_doses = pd.DataFrame(dict(gene_id=test))
        model_doses['technical_log_ratio'] = shift
        model_doses['applied_count_multiplier'] = model['fold']
        model_doses.to_csv(data/'Observed_Probe_Dose.csv',index=False)
        from publication.anchor_figures import draw
        status=decision(s,c,args.smoke)
        qc=dict(protocol, datasets=d_count*n_reps*len(DIRECTIONS)*2,
                valid_identity_checks=int(records.identity_error.notna().sum()),
                maximum_identity_error=float(records.identity_error.max()),
                minimum_simulated_anchors_detected=int(records.anchor_detected_min.min()),
                status=status)
        atomic_json(f/'summary.json',qc)
        draw(s,observed,f/'figures',args.formats,formal=not args.smoke)
        if not args.skip_replay:
            if not (f/'article'/'export_manifest.json').is_file():
                export_published(f/'article',args.formats,False)
        files=receipt(f,[p.relative_to(f).as_posix() for p in f.rglob('*') if p.is_file()
                         and p.name!='receipt.json'])
        atomic_json(finish,dict(identity=identity,files=files))
        print('CBC upgrade complete:', f, flush=True)
