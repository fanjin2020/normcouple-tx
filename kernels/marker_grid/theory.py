"""Finite deterministic bounds, not a universal prediction or FDR theorem."""
import numpy as np


def log_input_bound(counts, denominator_a, denominator_b, scale=1e6, pseudocount=1.):
    x = np.asarray(counts, float)
    da = np.asarray(denominator_a, float).reshape(-1, 1)
    db = np.asarray(denominator_b, float).reshape(-1, 1)
    if x.ndim != 2 or da.shape != (len(x), 1) or db.shape != da.shape:
        raise ValueError('Mismatched input dimensions')
    if not all(np.isfinite(a).all() for a in (x, da, db)) or (x < 0).any() or (da <= 0).any() or (db <= 0).any():
        raise ValueError('Finite nonnegative counts and positive denominators required')
    if not (np.isfinite(scale) and scale > 0 and np.isfinite(pseudocount) and pseudocount > 0):
        raise ValueError('Positive scale and pseudocount required')
    actual = np.abs(np.log2(scale * x / db + pseudocount) - np.log2(scale * x / da + pseudocount))
    sensitivity = scale * x / (scale * x + pseudocount * np.minimum(da, db))
    bound = sensitivity * np.abs(np.log(db / da)) / np.log(2.)
    return dict(actual=actual, bound=bound, max_violation=float(np.max(actual - bound)))


def coefficient_bound(y, exposure, truth, fit_a, fit_b):
    """Cauchy-Schwarz coefficient bound conditional on two estimable designs.

    For z_j=(I-P_j)x and w_j=z_j/(z_j'z_j):
      delta beta = (w_b-w_a)'(y-x*tau).
    No Lipschitz assumption about the upstream deconvolution algorithm is made.
    """
    y, x, tau = np.asarray(y), np.asarray(exposure), np.asarray(truth)
    za, zb = fit_a['residualized_exposure'], fit_b['residualized_exposure']
    qa, qb = float(za @ za), float(zb @ zb)
    if not (fit_a['basic'] and fit_b['basic'] and qa > 1e-12 and qb > 1e-12):
        return dict(eligible=False, reason='one_or_both_fits_not_estimable')
    dw = zb / qb - za / qa
    residual_truth = y - x[:, None] * tau
    actual = fit_b['beta'] - fit_a['beta']
    exact = dw @ residual_truth
    bound = np.linalg.norm(dw) * np.linalg.norm(residual_truth, axis=0)
    # Two triangle inequalities; choose the tighter valid orientation.
    dz, dq = np.linalg.norm(zb - za), abs(qb - qa)
    weight_bound = min(dz / qb + np.linalg.norm(za) * dq / (qa * qb),
                       dz / qa + np.linalg.norm(zb) * dq / (qa * qb))
    scale = 1 + np.max(np.abs(actual))
    if np.max(np.abs(actual - exact)) > 1e-8 * scale or np.max(np.abs(actual) - bound) > 1e-8 * scale:
        raise ValueError('Conditional coefficient identity/bound failed')
    if np.linalg.norm(dw) > weight_bound + 1e-8:
        raise ValueError('Residualized-exposure bound failed')
    return dict(eligible=True, identity_error=float(np.max(np.abs(actual - exact))),
                max_violation=float(np.max(np.abs(actual) - bound)),
                coefficient_change_max=float(np.max(np.abs(actual))), bound_max=float(bound.max()),
                weight_change=float(np.linalg.norm(dw)), weight_bound=float(weight_bound),
                z_change=float(dz), residual_ss_a=qa, residual_ss_b=qb)


def variance_decomposition(beta_error, leakage, noise):
    """Complete paired repetitions for one gene; sample covariance uses ddof=1."""
    a, l, e = (np.asarray(v, float) for v in (beta_error, leakage, noise))
    valid = np.isfinite(a) & np.isfinite(l) & np.isfinite(e)
    a, l, e = a[valid], l[valid], e[valid]
    if len(a) < 2:
        return dict(n=int(len(a)))
    covariance = np.cov(l, e, ddof=1)[0, 1]
    remainder = a.var(ddof=1) - l.var(ddof=1) - e.var(ddof=1) - 2 * covariance
    if abs(remainder) > 1e-7 * (1 + a.var(ddof=1)):
        raise ValueError('FWL variance decomposition failed')
    return dict(n=len(a), bias=float(a.mean()), leakage_mean=float(l.mean()), noise_mean=float(e.mean()),
                error_variance=float(a.var(ddof=1)), leakage_variance=float(l.var(ddof=1)),
                noise_variance=float(e.var(ddof=1)), twice_covariance=float(2 * covariance),
                variance_identity_error=float(abs(remainder)))
