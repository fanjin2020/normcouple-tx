"""HC3 in the actual nuisance column space, with auditable denominators."""
import numpy as np
from scipy import stats
from revision_core import M


def column_space(a):
    a = np.asarray(a, float)
    if not a.size:
        return np.empty((a.shape[0], 0)), np.array([])
    u, s, _ = np.linalg.svd(a, full_matrices=False)
    rank = int((s > max(1e-12, s[0] * 1e-10)).sum()) if len(s) else 0
    return u[:, :rank], s


def fit_hc3(y, exposure, props=None, covariates=None, rule="adaptive"):
    y = np.asarray(y, float)
    x = np.asarray(exposure, float)
    n, g = y.shape
    out = dict(beta=np.full(g, np.nan), se=np.full(g, np.nan),
               good=np.zeros(g, bool), basic=False, reason="", df=np.nan,
               leverage=np.nan, condition=np.nan, composition_rank=0,
               nuisance_rank=0, exposure_r2=np.nan, exposure_residual_ss=np.nan,
               exposure_estimable=False, weights=np.full(n, np.nan),
               residualized_exposure=np.full(n, np.nan))
    if not np.isfinite(y).all() or not np.isfinite(x).all():
        out["reason"] = "nonfinite_outcome_or_exposure"
        return out
    pieces = [np.ones((n, 1))]
    if covariates is not None:
        c = np.asarray(covariates, float)
        if not np.isfinite(c).all():
            out["reason"] = "nonfinite_covariates"
            return out
        pieces.append(c)
    if props is not None:
        if not np.isfinite(props).all():
            out["reason"] = "proportion_estimation_failed"
            return out
        pc, _ = column_space(props - np.mean(props, 0))
        out["composition_rank"] = pc.shape[1]
        pieces.append(pc)
    u, _ = column_space(np.column_stack(pieces))
    out["nuisance_rank"] = u.shape[1]
    z = x - u @ (u.T @ x)
    zz = float(z @ z)
    centered = float(np.sum((x - x.mean()) ** 2))
    out.update(exposure_residual_ss=zz, residualized_exposure=z,
               exposure_r2=1 - zz / centered if centered > 0 else np.nan)
    out["exposure_estimable"] = centered > 0 and zz > max(1e-12, centered * 1e-10)
    if props is not None and rule == "full_six" and pc.shape[1] != 6:
        out["reason"] = "composition_rank_below_six"
        return out
    if centered == 0 or zz <= max(1e-12, centered * 1e-10):
        out["reason"] = "exposure_not_estimable"
        return out
    out["exposure_estimable"] = True
    df = n - u.shape[1] - 1
    hat = np.sum(u * u, 1) + z * z / zz
    out.update(df=df, leverage=float(hat.max()),
               condition=float(np.linalg.cond(np.column_stack([u, x]))))
    if df < 2 or hat.max() >= 1 - 1e-9:
        out["reason"] = "HC3_not_estimable"
        return out
    w = z / zz
    beta = w @ y
    residual = y - u @ (u.T @ y) - z[:, None] * beta
    se = np.sqrt(np.sum((w[:, None] * residual / (1 - hat[:, None])) ** 2, 0))
    good = np.isfinite(beta) & np.isfinite(se) & (se > 1e-12)
    out.update(beta=beta, se=se, good=good, basic=True, weights=w,
               reason="" if good.all() else "gene_SE_unavailable")
    return out


def evaluate(data, props, rule="full_six", reason=""):
    fit = fit_hc3(data["y"], data["group"], props, rule=rule)
    valid = fit["good"] & fit["basic"]
    truth = data["truth"]
    null = truth == 0
    p = np.ones(len(truth))
    p[valid] = 2 * stats.t.sf(abs(fit["beta"][valid] / fit["se"][valid]), fit["df"])
    q = M.bh(p)
    hit = valid & (q <= .05)
    covered = valid & (abs(fit["beta"] - truth) <= stats.t.ppf(.975, fit["df"]) * fit["se"])
    leakage = fit["weights"] @ data["nuisance"]
    noise = fit["weights"] @ data["noise"]
    error = float(np.max(abs(fit["beta"][valid] - truth[valid] - leakage[valid] - noise[valid]))) if valid.any() else np.nan
    if valid.any() and error > 1e-7 * (1 + abs(fit["beta"][valid]).max()):
        raise ValueError("FWL identity failed")
    avg = lambda v, mask: float(np.mean(v[mask])) if mask.any() else np.nan
    row = dict(availability=float(valid.mean()), fdp=float(np.sum(hit & null) / max(1, hit.sum())),
               any_false=float(np.any(hit & null)), power=float(np.sum(hit & ~null) / np.sum(~null)) if (~null).any() else np.nan,
               coverage=avg(covered, valid), delivered_coverage=float(covered.mean()),
               null_coverage=avg(covered, valid & null), null_delivered_coverage=avg(covered, null),
               bias=avg(fit["beta"] - truth, valid), abs_leakage=avg(abs(leakage), valid),
               proportion_rmse=float(np.sqrt(np.mean((props - data["p"]) ** 2))) if props is not None and np.isfinite(props).all() else np.nan,
               mean_se=avg(fit["se"], valid), identity_error=error,
               reason=";".join(t for t in (reason, fit["reason"]) if t),
               design_condition_number=fit["condition"],
               **{k: fit[k] for k in ("composition_rank", "nuisance_rank", "exposure_r2", "exposure_residual_ss", "exposure_estimable", "df", "leverage")})
    arrays = dict(truth=truth, beta=np.where(valid, fit["beta"], np.nan), se=np.where(valid, fit["se"], np.nan),
                  leakage=np.where(valid, leakage, np.nan), noise_projection=np.where(valid, noise, np.nan),
                  available=valid.astype(float), covered=covered.astype(float), p=p, q=q)
    return row, arrays, fit


def mechanism_checks():
    rng = np.random.default_rng(2026101301)
    x = rng.uniform(1, 100, (23, 11))
    d = rng.uniform(2000, 4000, (23, 1))
    delta = rng.uniform(100, 800, (23, 1))
    a, c = 1., 1e6
    transform = lambda denominator: np.log2(c * x / denominator + a)
    observed = transform(d + delta) - transform(d)
    exact = np.log2((d / (d + delta)) * (c * x + a * (d + delta)) / (c * x + a * d))
    eps = 1e-5
    finite_diff = (transform(d * np.exp(eps)) - transform(d * np.exp(-eps))) / (2 * eps)
    derivative = -(c * x / (c * x + a * d)) / np.log(2)
    pure = np.log2(c * x / (d + delta)) - np.log2(c * x / d)
    marker = x / x.sum(1, keepdims=True)
    checks = dict(exact_finite_change_error=float(abs(exact - observed).max()),
                  derivative_error=float(abs(finite_diff - derivative).max()),
                  pure_log_row_shift_error=float(np.ptp(pure, axis=1).max()),
                  matched_a_c_error=float(abs(transform(d) - (np.log2(c * .01 * x / d + a * .01) - np.log2(.01))).max()),
                  marker_sum_error=float(abs(marker.sum(1) - 1).max()))
    if max(checks.values()) > 1e-8:
        raise ValueError("Analytic mechanism self-test failed")
    return checks
