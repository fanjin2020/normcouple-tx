"""Fixed-sample, fixed-response OLS/HC3 audit of supplied composition designs."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .cross_method import _setup_runtime, paired_metrics
from .io import atomic_json, lock, receipt, separate, sha


def aligned(frame, ids, label):
    if 'sample_id' not in frame or frame.sample_id.isna().any() or frame.sample_id.duplicated().any():
        raise ValueError(label + ' needs unique nonmissing sample_id values')
    if set(frame.sample_id) != set(ids):
        raise ValueError(label + ' sample IDs differ; no automatic intersection is performed')
    return frame.set_index('sample_id').loc[ids]


def execute(args):
    output = separate(args.output, args.outcomes, args.design, args.fractions)
    if output.exists() and any(output.iterdir()):
        raise ValueError('Use an empty sensitivity output directory')
    output.mkdir(parents=True, exist_ok=True)
    with lock(output / '.writer.lock'):
        c, _, s = _setup_runtime(output)
        outcomes = pd.read_csv(args.outcomes, dtype={'sample_id': str})
        if 'sample_id' not in outcomes or outcomes.sample_id.duplicated().any():
            raise ValueError('Outcome table needs one unique sample_id per row')
        ids = outcomes.sample_id.tolist()
        yframe = aligned(outcomes, ids, 'Outcomes')
        design = aligned(pd.read_csv(args.design, dtype={'sample_id': str}), ids, 'Design')
        fractions = pd.read_csv(args.fractions, dtype={'sample_id': str, 'method': str})
        if 'method' not in fractions or fractions.method.isna().any():
            raise ValueError('Fraction table needs a nonmissing method column')
        classes = [x for x in fractions.columns if x not in ('sample_id', 'method')]
        if len(classes) < 2 or fractions.method.nunique() < 2 or args.baseline not in set(fractions.method):
            raise ValueError('Need two cell classes, two methods and a present baseline method')
        required = [args.exposure, *args.covariates]
        if set(required) - set(design.columns) or len(set(required)) != len(required):
            raise ValueError('Exposure and covariates must name distinct, present numeric design columns')
        y = yframe.to_numpy(float)
        x = design[args.exposure].to_numpy(float)
        cov = design[args.covariates].to_numpy(float) if args.covariates else None
        if y.shape[1] < 1 or not np.isfinite(y).all() or not np.isfinite(x).all() or (cov is not None and not np.isfinite(cov).all()):
            raise ValueError('Outcomes and design must be finite numeric values on the intended analysis scale')
        fits, proportions, qvalues = {}, {}, {}
        summaries, genes, comparisons = [], [], []
        for method, frame in fractions.groupby('method', sort=True):
            props = aligned(frame, ids, method)[classes].to_numpy(float)
            if not np.isfinite(props).all() or (props < 0).any() or not np.allclose(props.sum(1), 1, atol=1e-6, rtol=0):
                raise ValueError('Invalid fractions for ' + method)
            # Close accepted rounded fractions to avoid a spurious extra rank.
            props = props / props.sum(1, keepdims=True)
            fit = s.fit_hc3(y, x, props, cov, rule='adaptive')
            available = fit['good'] & fit['basic']
            p = np.ones(y.shape[1])
            p[available] = 2 * stats.t.sf(np.abs(fit['beta'][available] / fit['se'][available]), fit['df'])
            q = c.M.bh(p)
            fits[method], proportions[method], qvalues[method] = fit, props, q
            summaries.append(dict(method=method, samples=len(ids), family_size=y.shape[1],
                available_tests=int(available.sum()), discoveries=int((q <= .05).sum()) if available.any() else np.nan,
                composition_rank=fit['composition_rank'], nuisance_rank=fit['nuisance_rank'],
                residual_df=fit['df'], exposure_information=fit['exposure_residual_ss'],
                maximum_leverage=fit['leverage'], reason=fit['reason']))
            for j, gene in enumerate(yframe.columns):
                genes.append(dict(method=method, gene_id=gene, available=bool(available[j]),
                    coefficient=fit['beta'][j], hc3_se=fit['se'][j], p=p[j], q=q[j]))
        baseline = fits[args.baseline]
        for method in fits:
            if method == args.baseline:
                continue
            fit = fits[method]
            change = paired_metrics(baseline, fit, y, proportions[args.baseline], proportions[method])
            estimable = bool(np.any(baseline['good'] & fit['good']) and baseline['basic'] and fit['basic'])
            change['bh_status_changes'] = int(np.sum((qvalues[args.baseline] <= .05) != (qvalues[method] <= .05))) if estimable else np.nan
            change['family_size'] = y.shape[1]
            comparisons.append(dict(baseline=args.baseline, method=method, **change))
        pd.DataFrame(genes).to_csv(output / 'gene_statistics.csv.gz', index=False, compression='gzip')
        pd.DataFrame(summaries).to_csv(output / 'method_summary.csv', index=False)
        pd.DataFrame(comparisons).to_csv(output / 'design_comparison.csv', index=False)
        atomic_json(output / 'audit.json', dict(
            inputs={name: sha(getattr(args, name)) for name in ('outcomes', 'design', 'fractions')},
            baseline=args.baseline, exposure=args.exposure, covariates=args.covariates,
            fraction_processing='Accepted row sums within 1e-6 are normalized to one before fitting.',
            inference='OLS with HC3 and residual-df t reference. BH over the complete supplied gene family; unavailable p=1.',
            interpretation='Sensitivity conditional on supplied designs. Threshold crossings are not validated changed effects or estimates of clinical FDR.',
            output=receipt(output, ['gene_statistics.csv.gz', 'method_summary.csv', 'design_comparison.csv'])))
    print('SENSITIVITY COMPLETE', output)
