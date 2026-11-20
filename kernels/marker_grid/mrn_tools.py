"""Reusable MRN inputs and actual-span HC3, without the simulation engine.

Arrays: marker_counts=(samples, markers), reference=(markers, cell types).
Reference entries are whole-library gene fractions, not logged expression.
These utilities neither select markers nor claim FDR control.
"""
import numpy as np


def mrn_inputs(marker_counts, reference, full_yield, anchor, scale=1e6, pseudocount=1.):
    """Return matched log inputs and reference marker yield for dtangle.

    ``anchor`` is a fixed reference-derived marker fraction, not estimated from
    targets. Keep it, marker IDs/order, reference and pseudocount fixed in a
    paired test-gene intervention. Zero marker libraries are errors, not zeros.
    """
    counts, a = np.asarray(marker_counts, float), np.asarray(reference, float)
    y = np.asarray(full_yield, float)
    if counts.ndim != 2 or a.ndim != 2 or counts.shape[1] != a.shape[0] or y.shape != (a.shape[1],):
        raise ValueError('MRN array dimensions do not match')
    if not all(np.isfinite(v).all() for v in (counts, a, y)) or (counts < 0).any() or (a < 0).any() or (y <= 0).any():
        raise ValueError('MRN requires finite nonnegative counts/reference and positive yields')
    if not (np.isfinite(anchor) and 0 < anchor <= 1 and np.isfinite(scale) and scale > 0 and np.isfinite(pseudocount) and pseudocount > 0):
        raise ValueError('Invalid anchor, scale or pseudocount')
    d, fraction = counts.sum(1), a.sum(0)
    if (d <= 0).any() or (fraction <= 0).any() or (fraction > 1 + 1e-8).any():
        raise ValueError('Empty marker library or invalid reference marker fraction')
    c = scale * anchor
    return dict(target_log2=np.log2(c * counts / d[:, None] + pseudocount),
                reference_log2=np.log2(c * (a / fraction).T + pseudocount),
                marker_yield=y * fraction, denominator=d)


def yield_convert(rna_fractions, reference_yield):
    """Convert RNA fractions to cell fractions using explicitly matched yields."""
    q, y = np.asarray(rna_fractions, float), np.asarray(reference_yield, float)
    if q.ndim != 2 or y.shape != (q.shape[1],) or not np.isfinite(q).all() or not np.isfinite(y).all():
        raise ValueError('Invalid fraction/yield dimensions or missing values')
    if (q < 0).any() or (y <= 0).any() or not np.allclose(q.sum(1), 1, atol=1e-6, rtol=0):
        raise ValueError('RNA fractions must sum to one; yields must be positive')
    p = q / y
    return p / p.sum(1, keepdims=True)


def column_space(a):
    a = np.asarray(a, float)
    if not a.size:
        return np.empty((a.shape[0], 0)), np.array([])
    u, s, _ = np.linalg.svd(a, full_matrices=False)
    rank = int((s > max(1e-12, s[0] * 1e-10)).sum()) if len(s) else 0
    return u[:, :rank], s


def actual_span_hc3(y, exposure, proportions=None, covariates=None):
    """Fit the exposure coefficient using the actual nuisance column space.

    Rank-deficient composition columns do not alone imply that exposure is
    unidentifiable. Returns availability and design diagnostics, not a joint
    model for uncertainty in estimated proportions. HC3 uses t(df) inference
    externally. Missing/nonfinite inputs are rejected; samples are not dropped.
    """
    y, x = np.asarray(y, float), np.asarray(exposure, float)
    if y.ndim != 2 or x.shape != (len(y),) or not np.isfinite(y).all() or not np.isfinite(x).all():
        raise ValueError('Expected finite sample-by-gene outcomes and exposure vector')
    n, g = y.shape
    parts = [np.ones((n, 1))]
    rank = 0
    for array, centered in ((covariates, False), (proportions, True)):
        if array is not None:
            a = np.asarray(array, float)
            if a.ndim != 2 or len(a) != n or not np.isfinite(a).all():
                raise ValueError('Invalid nuisance array')
            if centered:
                a, _ = column_space(a - a.mean(0))
                rank = a.shape[1]
            parts.append(a)
    u, _ = column_space(np.column_stack(parts))
    z = x - u @ (u.T @ x)
    zz, total = float(z @ z), float(np.sum((x - x.mean()) ** 2))
    result = dict(beta=np.full(g, np.nan), se=np.full(g, np.nan), available=np.zeros(g, bool),
                  composition_rank=rank, nuisance_rank=u.shape[1], exposure_residual_ss=zz,
                  exposure_r2=1 - zz / total if total > 0 else np.nan,
                  exposure_estimable=bool(total > 0 and zz > max(1e-12, total * 1e-10)),
                  df=n - u.shape[1] - 1, reason='', residualized_exposure=z)
    if not result['exposure_estimable']:
        result['reason'] = 'exposure_not_estimable'
        return result
    hat = np.sum(u * u, 1) + z * z / zz
    result['leverage'] = float(hat.max())
    if result['df'] < 2 or hat.max() >= 1 - 1e-9:
        result['reason'] = 'HC3_not_estimable'
        return result
    w = z / zz
    beta = w @ y
    residual = y - u @ (u.T @ y) - z[:, None] * beta
    contributions = (w[:, None] * residual / (1 - hat[:, None])) ** 2
    se = np.sqrt(contributions.sum(0))
    available = np.isfinite(beta) & np.isfinite(se) & (se > 1e-12)
    result.update(beta=beta, se=se, available=available, weights=w,
                  hc3_contributions=contributions, residual=residual)
    if not available.all():
        result['reason'] = 'gene_SE_unavailable'
    return result
